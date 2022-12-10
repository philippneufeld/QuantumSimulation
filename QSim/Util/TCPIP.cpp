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
            const std::vector<TCPIPServerSocket*>& acceptable, 
            const std::vector<TCPIPConnection*>& readable, 
            const std::vector<TCPIPConnection*>& writable, 
            const std::vector<TCPIPConnection*>& exceptional)
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


    //
    // SocketDataPackage
    //

    SocketDataPackage::SocketDataPackage() 
        : m_size(0), m_pData(nullptr), m_status(SocketDataPackageStatus_OK), m_messageId(0) {}

    SocketDataPackage::SocketDataPackage(std::uint64_t size) 
        : SocketDataPackage() 
    {
        Allocate(size);
    }

    SocketDataPackage::SocketDataPackage(const Header_t& header)
        : SocketDataPackage()
    {
        if (IsValidHeader(header))
        {
            m_status = GetStatusFromHeader(header);
            m_messageId = GetMessageIdFromHeader(header);
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
        m_status = rhs.m_status;
        m_messageId = rhs.m_messageId;
        if (rhs.m_size > 0)
        {
            Allocate(rhs.m_size);
            std::copy(rhs.m_pData, rhs.m_pData+rhs.m_size, m_pData);
        }
    }
    
    SocketDataPackage::SocketDataPackage(SocketDataPackage&& rhs)
        : m_pData(rhs.m_pData), m_size(rhs.m_size), m_status(rhs.m_status), m_messageId(rhs.m_messageId)
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

    SocketDataPackage::operator bool() const
    {
        return m_status == SocketDataPackageStatus_OK;
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

    bool SocketDataPackage::IsValidHeader(const Header_t& header)
    {
        std::uint64_t npid = 0;
        std::copy_n(header.data(), 7, reinterpret_cast<std::uint8_t*>(&npid) + 1);
        std::uint64_t pid = ntohll(npid);
        return (pid == s_protocolId);
    }

    std::uint8_t SocketDataPackage::GetStatusFromHeader(const Header_t& header)
    {
        return header[7];
    }
    
    std::uint32_t SocketDataPackage::GetMessageIdFromHeader(const Header_t& header)
    {
        std::uint32_t nmid = 0;
        std::copy_n(header.data() + 16, sizeof(nmid), reinterpret_cast<std::uint8_t*>(&nmid));
        return ntohl(nmid);
    }
    
    std::uint64_t SocketDataPackage::GetSizeFromHeader(const Header_t& header)
    {
        std::uint64_t nsize = 0;
        std::copy_n(header.data() + 8, sizeof(nsize), reinterpret_cast<std::uint8_t*>(&nsize));
        return ntohll(nsize);
    }

    SocketDataPackage::Header_t SocketDataPackage::GenerateHeader(std::uint64_t size, std::uint8_t status, std::uint32_t msgId)
    {
        Header_t header;

        std::uint64_t npid = htonll(s_protocolId);
        std::uint64_t nsize = htonll(size);
        std::uint32_t nmid = htonl(msgId);
        
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&npid) + 1, 7, header.data());
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&status), 1, header.data() + 7);
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&nsize), 8, header.data() + 8);
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&nmid), 4, header.data() + 16);

        return header;
    }

    SocketDataPackage::Header_t SocketDataPackage::GetHeader() const
    {
        return GenerateHeader(m_size, m_status, m_messageId);
    }
    
    SocketDataPackage SocketDataPackage::CreateError(std::uint8_t status)
    {
        SocketDataPackage package;
        package.SetStatus(status);
        return package;
    }


    //
    //
    //

    void TCPIPConnectionHandler::AddConnection(TCPIPConnection* pConn)
    {
        auto it = std::find(m_connections.begin(), m_connections.end(), pConn);
        if (it == m_connections.end())
            m_connections.push_back(pConn);
    }

    void TCPIPConnectionHandler::RemoveConnection(TCPIPConnection* pConn)
    {
        auto it = std::find(m_connections.begin(), m_connections.end(), pConn);
        if (it != m_connections.end())
            m_connections.erase(it);
    }

    void TCPIPConnectionHandler::AddServer(TCPIPServerSocket* pServer)
    {
        auto it = std::find(m_servers.begin(), m_servers.end(), pServer);
        if (it == m_servers.end())
            m_servers.push_back(pServer);
    }

    void TCPIPConnectionHandler::RemoveServer(TCPIPServerSocket* pServer)
    {
        auto it = std::find(m_servers.begin(), m_servers.end(), pServer);
        if (it != m_servers.end())
            m_servers.erase(it);
    }

    void TCPIPConnectionHandler::RunHandler()
    {
        DoRun();
    }

    void TCPIPConnectionHandler::StopHandler()
    {
        m_keepRunning = true;
    }

    void TCPIPConnectionHandler::DoRunHandler()
    {
        std::map<TCPIPConnection*, std::tuple<SocketDataPackage, std::uint64_t>> readBuffer;
        std::set<TCPIPConnection*> killed;

        std::vector<TCPIPServerSocket*> checkAcceptable;
        std::vector<TCPIPConnection*> checkReadable;
        std::vector<TCPIPConnection*> checkWritable;

        m_keepRunning = true;
        while (m_keepRunning)
        { 
            checkAcceptable.clear();
            checkReadable.clear();
            checkWritable.clear();
            killed.clear();
 
            // prepare lists of connections to be checked
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                checkAcceptable.reserve(m_servers.size());
                checkReadable.reserve(m_connections.size());
                checkWritable.reserve(m_writeBuffer.size());
                std::copy(m_servers.begin(), m_servers.end(), std::back_inserter(checkAcceptable));
                std::copy(m_connections.begin(), m_connections.end(), std::back_inserter(checkReadable));
                std::transform(m_writeBuffer.begin(), m_writeBuffer.end(), 
                    std::back_inserter(checkWritable), [](auto el) { return el.first; });
            }

            // block until either accept, read or write can be performed
            // TODO: wakup with pipe
            auto [ac, re, wr, ex] = TCPIPServerSocket::Select(checkAcceptable, checkReadable, checkWritable, {});
            
            // handle new connections
            for (TCPIPServerSocket* pServer: ac)
                this->OnConnectionAcceptale(pServer);

            // handle receiving messages
            for (TCPIPConnection* pConn: re)
            {
                if (killed.count(pConn) != 0)
                    continue;
                
                // check if data in read buffer
                if (readBuffer.count(pConn) == 0)
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
                    
                    readBuffer[pConn] = std::make_tuple(SocketDataPackage(header), std::uint64_t(0));                    
                }
                else
                {
                    auto& [buffer, alreadyRead] = readBuffer[pConn];
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
                            SocketDataPackage resp = OnMessageReceived(pConn, std::move(buffer));
                            m_writeBuffer[pConn].push_back(std::move(resp));
                            readBuffer.erase(readBuffer.find(pConn));
                        }
                    }
                }
            }

            // handle sending messages
            for (TCPIPConnection* pConn: wr)
            {
                // access to write buffer is thread-safe without lock in this case,
                // since all other threads may only append to the list 
                // (this does not invalidate the iterators to the existing items)
                if (killed.count(pConn) != 0)
                    continue;

                auto& outputQueue = m_writeBuffer[pConn];

                while (true)
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    if (outputQueue.empty())
                        break;
                    SocketDataPackage package = std::move(outputQueue.front());
                    outputQueue.pop();
                    lock.unlock();

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
            }
            
            // purge write buffer of empty queues
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                for (auto it = m_writeBuffer.begin(); it != m_writeBuffer.end();)
                {
                    if (it->second.empty())
                        it = m_writeBuffer.erase(it);
                    else
                        it = std::next(it);
                }
                
            }

            // handle killed connections
            for (TCPIPConnection* pConn: killed)
            {
                this->RemoveConnection(pConn);
                pConn->Close();
                this->OnConnectionClosed(pConn);

                // clear read buffer
                auto it1 = readBuffer.find(pConn);
                if (it1 != readBuffer.end())
                    readBuffer.erase(it1);

                // clear write buffer
                std::unique_lock<std::mutex> lock(m_mutex);
                auto it2 = m_writeBuffer.find(pConn);
                if (it2 != m_writeBuffer.end())
                    m_writeBuffer.erase(it2);
            }
        }
    }

    //
    // TCPIPServer
    //

    TCPIPServer::TCPIPServer()
        : m_bRunning(false) { }
    
    TCPIPServer::~TCPIPServer()
    {
    }

    bool TCPIPServer::Run(short port)
    {
        // start socket server
        TCPIPServerSocket server;
        if (!server.Bind(8000))
            return false;
        if (!server.Listen(5))
            return false;

        this->AddServer(&server);
        this->RunHandler();
    }


    //
    // TCPIPClient
    //

    SocketDataPackage TCPIPClient::Query(const void* data, std::uint64_t n, std::uint32_t msgId)
    {
        SocketDataPackage::Header_t header;
        
        // send package header
        header = SocketDataPackage::GenerateHeader(n, SocketDataPackageStatus_OK, msgId);
        if (TCPIPClientSocket::Send(header.data(), sizeof(header)) != sizeof(header))
        {
            Close();
            return SocketDataPackage::CreateError(SocketDataPackageStatus_IO_ERROR);
        }

        // send package data
        if (n > 0)
        {
            if (TCPIPClientSocket::Send(data, n) != n)
            {
                Close();
                return SocketDataPackage::CreateError(SocketDataPackageStatus_IO_ERROR);
            }
        }

        // receiving package header of response
        std::uint64_t cnt = TCPIPClientSocket::Recv(header.data(), sizeof(header));
        if (cnt != sizeof(header) || !SocketDataPackage::IsValidHeader(header))
        {
            Close();
            return SocketDataPackage::CreateError(SocketDataPackageStatus_INVALID_PACKAGE);
        }
        
        SocketDataPackage package(header);
        for (std::size_t read = 0; read < package.GetSize();)
        {
            cnt = TCPIPClientSocket::Recv(package.GetData() + read, package.GetSize()-read);
            if (cnt <= 0)
            {
                Close();
                return SocketDataPackage::CreateError(SocketDataPackageStatus_IO_ERROR);
            }
            read += cnt;
        }

        return package;
    }

    //
    // TCPIPMultiClient
    //

    TCPIPMultiClient::TCPIPMultiClient() { }

    std::size_t TCPIPMultiClient::Connect(const std::string& ip, short port)
    {
        auto pConn = new TCPIPClientSocket();
        if (!pConn->Connect(ip, port))
            return 0;
        
        m_connections.push_back(pConn);
        return GetIdFromConnection(pConn);
    }
    
    std::size_t TCPIPMultiClient::ConnectHostname(const std::string& hostname, short port)
    {
        return Connect(TCPIPClientSocket::GetHostByName(hostname), port);
    }

    std::size_t TCPIPMultiClient::GetIdFromConnection(TCPIPConnection* pConnection) const
    {
        static_assert(sizeof(std::size_t) == sizeof(pConnection), "id type not large enough");
        return reinterpret_cast<std::size_t>(pConnection);
    }
}
