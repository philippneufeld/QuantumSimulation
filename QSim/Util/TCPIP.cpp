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

// https://stackoverflow.com/questions/3022552/is-there-any-standard-htonl-like-function-for-64-bits-integers-in-c
#ifdef __BIG_ENDIAN__
#   define htonll(x) (x)
#   define ntohll(x) (x)
#else
#   define htonll(x) (((std::uint64_t)htonl((x) & 0xFFFFFFFF) << 32) | htonl((x) >> 32))
#   define ntohll(x) (((std::uint64_t)ntohl((x) & 0xFFFFFFFF) << 32) | ntohl((x) >> 32))
#endif


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

    std::int64_t TCPIPConnection::Send(const void* data, std::size_t n)
    {
        if (!IsValid())
            return 0;

        return send(GetHandle(), data, n, 0);
    }
    
    std::int64_t TCPIPConnection::Recv(void* buffer, std::size_t buffSize)
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
        if (!pRes || pRes->h_length < 4)
            return false;
        char** addresses = pRes->h_addr_list;
        
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
        : m_size(0), m_pData(nullptr), m_status(0) {}

    SocketDataPackage::SocketDataPackage(std::uint64_t size) 
        : SocketDataPackage(size, 0) {}
    
    SocketDataPackage::SocketDataPackage(std::uint64_t size, std::uint8_t status) 
        : SocketDataPackage() 
    {
        m_status = status;
        Allocate(size);
    }

    SocketDataPackage::SocketDataPackage(const std::array<std::uint8_t, 16>& header)
        : SocketDataPackage()
    {
        if (IsValidHeader(header))
        {
            m_status = GetStatusFromHeader(header);
            Allocate(GetSizeFromHeader(header));
        }
    }
    
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

    bool SocketDataPackage::IsValidHeader(const std::array<std::uint8_t, 16>& header)
    {
        std::uint64_t npid = 0;
        std::copy_n(header.data(), 7, reinterpret_cast<std::uint8_t*>(&npid) + 1);
        std::uint64_t pid = ntohll(npid);
        return (pid == s_protocolId);
    }

    std::uint8_t SocketDataPackage::GetStatusFromHeader(const std::array<std::uint8_t, 16>& header)
    {
        return header[7];
    }
    
    std::uint64_t SocketDataPackage::GetSizeFromHeader(const std::array<std::uint8_t, 16>& header)
    {
        std::uint64_t nsize = 0;
        std::copy_n(header.data() + 8, 8, reinterpret_cast<std::uint8_t*>(&nsize));
        return ntohll(nsize);
    }

    SocketDataPackage::Header_t SocketDataPackage::GenerateHeader(std::uint64_t size, std::uint8_t status)
    {
        Header_t header;

        std::uint64_t npid = htonll(s_protocolId);
        std::uint64_t nsize = htonll(size);
        
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&npid) + 1, 7, header.data());
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&status), 1, header.data() + 7);
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&nsize), 8, header.data() + 8);

        return header;
    }

    SocketDataPackage::Header_t SocketDataPackage::GetHeader() const
    {
        return GenerateHeader(m_size, m_status);
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

        std::set<TCPIPConnection*> killed;
        while (true)
        {
            killed.clear();
            auto [ac, re, wr, ex] = TCPIPServerSocket::Select({&server}, m_pClients, m_pClients, {});
            
            // handle new connections
            if (std::find(ac.begin(), ac.end(), &server) != ac.end())
            {
                auto [cid, ip] = server.Accept();
                auto pConn = new TCPIPConnection(cid);
                if (pConn->IsValid())
                {
                    if (OnClientConnected(GetIdFromConnection(pConn), ip))
                    {
                        m_pClients.push_back(pConn);
                    }
                }
            }

            // handle receiving messages
            for (TCPIPConnection* pConn: re)
            {
                if (killed.count(pConn) != 0)
                    continue;
                
                auto id = GetIdFromConnection(pConn);

                // check if data in read buffer
                if (m_readBuffer.count(pConn) == 0)
                {
                    SocketDataPackage::Header_t header;
                    if (pConn->Recv(header.data(), sizeof(header)) != sizeof(header))
                    {
                        killed.insert(pConn);
                        continue;
                    }

                    if (!SocketDataPackage::IsValidHeader(header))
                    {
                        killed.insert(pConn);
                        continue;
                    }
                    
                    m_readBuffer[pConn] = std::make_tuple(SocketDataPackage(header), std::uint64_t(0));                    
                }
                else
                {
                    auto& [buffer, alreadyRead] = m_readBuffer[pConn];
                    std::uint64_t remaining = buffer.GetSize() - alreadyRead;

                    std::size_t cnt = pConn->Recv(buffer.GetData() + alreadyRead, remaining);
                    if (cnt <= 0 || cnt > remaining)
                    {
                        killed.insert(pConn);
                        continue;
                    }
                    else
                    {
                        alreadyRead += cnt; // update read buffer (reference)
                        if (alreadyRead == buffer.GetSize())
                        {
                            SocketDataPackage resp = OnMessageReceived(id, std::move(buffer));
                            m_writeBuffer[pConn].push_back(std::move(resp));
                            m_readBuffer.erase(m_readBuffer.find(pConn));
                        }
                    }
                }
            }

            // handle sending messages
            for (TCPIPConnection* pConn: wr)
            {
                if (killed.count(pConn) != 0)
                    continue;

                auto& outputBuff = m_writeBuffer[pConn];
                for (auto& package: outputBuff)
                {
                    
                    // send package header
                    SocketDataPackage::Header_t header = package.GetHeader();
                    if (pConn->Send(header.data(), sizeof(header)) != sizeof(header))
                    {
                        killed.insert(pConn);
                        break;
                    }

                    // send package data
                    std::uint64_t size = package.GetSize();
                    if (size > 0)
                    {
                        if (pConn->Send(package.GetData(), size) != size)
                        {
                            killed.insert(pConn);
                            break;
                        }
                    }
                }
                outputBuff.clear();
            }

            // handle killed connections
            for (TCPIPConnection* pConn: killed)
            {
                auto id = GetIdFromConnection(pConn);
                CloseClientConnection(pConn);
                OnClientDisconnected(id);
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

        auto it1 = m_readBuffer.find(pClient);
        if (it1 != m_readBuffer.end())
            m_readBuffer.erase(it1);

        auto it2 = m_writeBuffer.find(pClient);
        if (it2 != m_writeBuffer.end())
            m_writeBuffer.erase(it2);

        delete pClient;
    }

    std::size_t TCPIPServer::GetIdFromConnection(TCPIPConnection* pClient) const
    {
        static_assert(sizeof(std::size_t) == sizeof(TCPIPConnection*), "id type not large enough");
        return reinterpret_cast<std::size_t>(pClient);
    }

    std::optional<SocketDataPackage> TCPIPClient::Query(const void* data, std::uint64_t n)
    {
        SocketDataPackage::Header_t header;
        
        // send package header
        header = SocketDataPackage::GenerateHeader(n, 0);
        if (TCPIPClientSocket::Send(header.data(), sizeof(header)) != sizeof(header))
        {
            Close();
            return std::nullopt;
        }

        // send package data
        if (n > 0)
        {
            if (TCPIPClientSocket::Send(data, n) != n)
            {
                Close();
                return std::nullopt;
            }
        }

        // receiving package header of response
        std::uint64_t cnt = TCPIPClientSocket::Recv(header.data(), sizeof(header));
        if (cnt != sizeof(header) || !SocketDataPackage::IsValidHeader(header))
        {
            Close();
            return std::nullopt;
        }
        
        SocketDataPackage package(header);
        for (std::size_t read = 0; read < package.GetSize();)
        {
            cnt = TCPIPClientSocket::Recv(package.GetData() + read, package.GetSize()-read);
            if (cnt <= 0)
            {
                Close();
                return std::nullopt;
            }
            read += cnt;
        }

        return package;
    }
}
