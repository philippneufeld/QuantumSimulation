// Philipp Neufeld, 2021-2022

#include <algorithm>
#include <cstdlib>
#include <sstream>

#include <unistd.h> // close
#include <sys/socket.h> // socket functions
#include <sys/select.h> // select
#include <netinet/in.h> // sockaddr_in
#include <arpa/inet.h> // inet_ntoa
#include <netdb.h> // gethostbyname

#include "TCPIP.h"

#include <iostream>

namespace QSim
{
    
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

    std::size_t TCPIPConnection::Send(const void* data, std::size_t n)
    {
        if (!IsValid())
            return 0;

        return send(GetHandle(), data, n, 0);
    }
    
    std::size_t TCPIPConnection::Recv(void* buffer, std::size_t buffSize)
    {
        if (!IsValid())
            return 0;

        return recv(GetHandle(), buffer, buffSize, 0);
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
            const std::list<TCPIPServerSocket*>& acceptable, 
            const std::list<TCPIPConnection*>& readable, 
            const std::list<TCPIPConnection*>& writable, 
            const std::list<TCPIPConnection*>& exceptional)
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
            !acceptable.empty() || !readable.empty() ? &readFds : nullptr, 
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

    bool TCPIPClientSocket::ConnectHostname(const std::string& hostname, short port)
    {
        constexpr int bufferLen = 256;
        char buffer[bufferLen] = "";
        gethostname(buffer, bufferLen - 1);

        hostent* pRes = gethostbyname(hostname.c_str());
        char** addresses = pRes->h_addr_list;
        if (pRes->h_length < 4)
            return false;

        std::string ip = 
            std::to_string((int)(unsigned char)addresses[0][0]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][1]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][2]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][3]);
        
        return Connect(ip, port);
    }


    //
    // SocketDataPackage
    //

    SocketDataPackage::SocketDataPackage() 
        : m_size(0), m_pData(nullptr) {}
    
    SocketDataPackage::~SocketDataPackage() 
    { 
        Allocate(0);
    }

    SocketDataPackage::SocketDataPackage(const SocketDataPackage& rhs)
        : SocketDataPackage()
    {
        if (rhs.m_size > 0)
        {
            Allocate(rhs.m_size);
            std::copy(rhs.m_pData, rhs.m_pData+rhs.m_size, m_pData);
        }
    }
    
    SocketDataPackage::SocketDataPackage(SocketDataPackage&& rhs)
        : m_pData(rhs.m_pData), m_size(rhs.m_size)
    {
        rhs.m_pData = nullptr;
        rhs.m_size = 0;
    }

    SocketDataPackage& SocketDataPackage::operator=(const SocketDataPackage& rhs)
    {
        SocketDataPackage tmp(rhs);
        std::swap(*this, tmp);
        return *this;
    }

    SocketDataPackage& SocketDataPackage::operator=(SocketDataPackage&& rhs)
    {
        std::swap(m_pData, rhs.m_pData);
        std::swap(m_size, rhs.m_size);
        return *this;
    }
    
    bool SocketDataPackage::Allocate(std::uint64_t size)
    {
        if (size == m_size)
            return true;

        if (m_pData)
            delete[] m_pData;
        m_pData = nullptr;
        m_size = 0;      
        
        if (size > 0)
        {
            m_pData = new std::uint8_t[size];
            m_size = size;
            if (!m_pData)
            {
                m_pData = nullptr;
                m_size = 0;
            }
        }

        return (size == m_size);
    }


    //
    // TCPIPServer
    //

    TCPIPServer::TCPIPServer()
        : m_bRunning(false) { }
    
    TCPIPServer::~TCPIPServer()
    {
        Purge();
    }

    bool TCPIPServer::Run(short port)
    {
        // clear all buffers
        Purge();
        
        // start socket server
        TCPIPServerSocket server;
        if (!server.Bind(8000))
            return false;
        if (!server.Listen(5))
            return false;

        while (true)
        {
            auto [ac, re, wr, ex] = TCPIPServerSocket::Select({&server}, m_pClients, {}, {});
            
            // handle new connections
            if (std::find(ac.begin(), ac.end(), &server) != ac.end())
            {
                auto [cid, ip] = server.Accept();
                auto pClient = new TCPIPConnection(cid);
                if (pClient->IsValid())
                {
                    if (this->OnClientConnected(ip))
                        m_pClients.push_back(pClient);
                }
            }

            // handle receiving messages
            for (TCPIPConnection* pConn: re)
            {
                constexpr int buffsize = 4096;
                std::uint8_t buff[buffsize] = "";
                int cnt = pConn->Recv(buff, buffsize);

                if (cnt <= 0)
                {
                    CloseClientConnection(pConn);
                    OnClientDisconnected();
                }
                else
                {
                    SocketDataPackage data;
                    data.Allocate(cnt);
                    std::copy_n(buff, cnt, data.GetData());
                    OnMessageReceived(std::move(data));
                }
            }
            
        }

        return 0;
    }

    void TCPIPServer::Purge()
    {
        for (TCPIPConnection* pClient: m_pClients)
            delete pClient;
        m_pClients.clear();
    }

    void TCPIPServer::CloseClientConnection(TCPIPConnection* pClient)
    {
        m_pClients.remove(pClient);
        delete pClient;
    }


}
