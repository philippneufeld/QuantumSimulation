// Philipp Neufeld, 2021-2022

#include <algorithm>
#include <cstdlib>
#include <sstream>

#include <fcntl.h>

#include <unistd.h> // close
#include <sys/socket.h> // socket functions
#include <sys/select.h> // select
#include <netinet/in.h> // sockaddr_in
#include <arpa/inet.h> // inet_ntoa
#include <netdb.h> // gethostbyname

#include "TCPIP.h"

namespace QSim
{

    //
    // IOEvent
    //

    IOEvent::IOEvent()
    {
        pipe(m_fd);
        fcntl(m_fd[0], F_SETFL, O_NONBLOCK);
    }

    IOEvent::~IOEvent()
    {
        for (int i=0; i<2; i++)
        {
            if (m_fd[i] >= 0)
                close(m_fd[i]);
            m_fd[i] = -1;
        }
    }

    void IOEvent::Set()
    {
        write(m_fd[1], "x", 1);
    }

    void IOEvent::Reset()
    {
        std::uint8_t buff = 0;
        while (read(m_fd[0], &buff, 1) == 1);
    }

    //
    // TCPIPSocket
    //
     
    TCPIPSocket::TCPIPSocket()
        : TCPIPSocket(socket(AF_INET, SOCK_STREAM, 0)) {}

    TCPIPSocket::TCPIPSocket(int handle)
        : m_handle(handle) {}

    TCPIPSocket::~TCPIPSocket()
    {
        Close();
    }
    
    TCPIPSocket::TCPIPSocket(TCPIPSocket&& rhs)
        : TCPIPSocket()
    {
        std::swap(m_handle, rhs.m_handle);
    }

    TCPIPSocket& TCPIPSocket::operator=(TCPIPSocket&& rhs)
    {
        std::swap(m_handle, rhs.m_handle);
        return *this;
    }

    void TCPIPSocket::Close()
    {
        if (IsValid())
            close(m_handle);
        m_handle = -1;
    }


    //
    // TCPIPConnection
    //

    TCPIPConnection::TCPIPConnection() 
        : TCPIPSocket() {}

    TCPIPConnection::TCPIPConnection(int handle)
        : TCPIPSocket(handle) {}

    std::make_signed_t<std::size_t> TCPIPConnection::Send(const void* data, std::size_t n)
    {
        if (!IsValid())
            return 0;

        return send(GetHandle(), data, n, 0);
    }
    
    std::make_signed_t<std::size_t> TCPIPConnection::Recv(void* buffer, std::size_t buffSize)
    {
        if (!IsValid())
            return 0;

        return recv(GetHandle(), buffer, buffSize, 0);
    }

    std::make_signed_t<std::size_t> TCPIPConnection::SendNonBlock(const void* data, std::size_t n)
    {
        if (!IsValid())
            return 0;

        return send(GetHandle(), data, n, MSG_DONTWAIT);
    }
    
    std::make_signed_t<std::size_t> TCPIPConnection::RecvNonBlock(void* buffer, std::size_t buffSize)
    {
        if (!IsValid())
            return 0;

        return recv(GetHandle(), buffer, buffSize, MSG_DONTWAIT);
    }

    //
    // TCPIPServerSocket
    //

    TCPIPServerSocket::TCPIPServerSocket()
        : TCPIPSocket() {}

    bool TCPIPServerSocket::Bind(short port)
    {
        if (!IsValid())
            return false;

        // lets server reuse previously bound address
        const int enable = 1;
        if (setsockopt(GetHandle(), SOL_SOCKET, SO_REUSEADDR, &enable, sizeof(enable)) < 0)
            return false;

        sockaddr_in addr = { 0 };
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY; // bind to all interfaces
        addr.sin_port = htons(port);

        if (bind(GetHandle(), reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
            return false;

        return true;
    }

    bool TCPIPServerSocket::Listen(int queue_length)
    {
        if (!IsValid())
            return false;

        if (listen(GetHandle(), queue_length) < 0)
            return false;

        return true;
    }

    std::tuple<int, std::string> TCPIPServerSocket::Accept()
    {
        sockaddr_in addr = { 0 };
        addr.sin_family = AF_INET;
        socklen_t len = sizeof(addr);

        int clientHandle = accept(GetHandle(), reinterpret_cast<sockaddr*>(&addr), &len);
        std::string ip = inet_ntoa(addr.sin_addr);

        return std::make_tuple(clientHandle, ip);
    }

    std::tuple<TCPIPConnection, std::string> TCPIPServerSocket::AcceptC()
    {
        auto [cid, ip] = Accept();
        return std::make_tuple(TCPIPConnection(cid), ip);
    }

    std::tuple<
        std::vector<TCPIPServerSocket*>, std::vector<TCPIPConnection*>, 
        std::vector<TCPIPConnection*>, std::vector<TCPIPConnection*>
    >
        TCPIPServerSocket::Select(
            const std::vector<TCPIPServerSocket*>& acceptable, 
            const std::vector<TCPIPConnection*>& readable, 
            const std::vector<TCPIPConnection*>& writable, 
            const std::vector<TCPIPConnection*>& exceptional,
            IOEvent& wakeupSignal)
    {
        // initialize file descriptor sets
        fd_set readFds, writeFds, excFds;

        // clear sets
        FD_ZERO(&readFds);
        FD_ZERO(&writeFds);
        FD_ZERO(&excFds);

        bool bA = !acceptable.empty();
        bool bR = !readable.empty();
        bool bW = !writable.empty();
        bool bE = !exceptional.empty();

        // populate sets
        FD_SET(wakeupSignal.GetFileDescriptor(), &readFds);
        for (const TCPIPServerSocket* pServer : acceptable) 
            FD_SET(pServer->GetHandle(), &readFds);
        for (const TCPIPConnection* pConn : readable) 
            FD_SET(pConn->GetHandle(), &readFds);
        for (const TCPIPConnection* pConn : writable) 
            FD_SET(pConn->GetHandle(), &writeFds);
        for (const TCPIPConnection* pConn : exceptional) 
            FD_SET(pConn->GetHandle(), &excFds);
    
        // Maybe error handling here?
        int cnt = select(FD_SETSIZE, 
            &readFds, 
            !writable.empty() ? &writeFds : nullptr, 
            !exceptional.empty() ? &excFds : nullptr, nullptr);

        
        // check result for available file descriptors
        std::vector<TCPIPServerSocket*> avA;
        std::vector<TCPIPConnection*> avR;
        std::vector<TCPIPConnection*> avW;
        std::vector<TCPIPConnection*> avE;

        // convenience function to reduce duplicate code
        auto populate = [&](auto& av, const auto& src, fd_set* pSet)
        {
            if (!src.empty())
            {
                av.reserve(cnt);
                for (auto* pEl : src)
                {
                    if (FD_ISSET(pEl->GetHandle(), pSet))
                        av.push_back(pEl);
                }
            }
        };

        if (cnt > 0)
        {
            populate(avA, acceptable, &readFds);
            populate(avR, readable, &readFds);
            populate(avW, writable, &writeFds);
            populate(avE, exceptional, &excFds);
        }
        
        // return list of available servers and connections
        return std::make_tuple(avA, avR, avW, avE);
    }


    //
    // TCPIPClientSocket
    //

    TCPIPClientSocket::TCPIPClientSocket()
        : TCPIPConnection() {}

    bool TCPIPClientSocket::Connect(const std::string& ip, short port)
    {
        sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);

        std::istringstream iss(ip);
        for (int i=0; i<4; i++)
        {
            std::string tmp;
            std::getline(iss, tmp, '.');
            reinterpret_cast<std::uint8_t*>(&addr.sin_addr.s_addr)[i] = std::atoi(tmp.c_str());
        }

        if (connect(GetHandle(), reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
            return false;

        return true;
    }

    std::string TCPIPClientSocket::GetHostByName(const std::string& hostname)
    {
        constexpr int bufferLen = 256;
        char buffer[bufferLen] = "";
        gethostname(buffer, bufferLen - 1);

        hostent* pRes = gethostbyname(hostname.c_str());
        if (!pRes || pRes->h_length < 4)
            return "";
        char** addresses = pRes->h_addr_list;
        
        std::string ip = 
            std::to_string((int)(unsigned char)addresses[0][0]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][1]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][2]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][3]);

        return ip;
    }

}
