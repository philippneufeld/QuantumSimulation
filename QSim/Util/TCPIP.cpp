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
    // TCPIPConnectionHandler
    //

    TCPIPConnectionHandler::TCPIPConnectionHandler() { }
    
    TCPIPConnectionHandler::~TCPIPConnectionHandler() 
    {
        while (!m_connections.empty())
            RemoveConnection(*m_connections.begin());

        while (!m_servers.empty())
            RemoveServer(m_servers.front());
    }

    void TCPIPConnectionHandler::AddConnection(TCPIPConnection* pConn)
    {
        if (!pConn)
            return;
        
        std::unique_lock<std::mutex> lock(m_mutex);
        auto it = m_connections.find(pConn);
        if (it == m_connections.end())
        {
            m_connections.insert(pConn);
            m_wakeupSignal.Set();
        }
    }

    void TCPIPConnectionHandler::RemoveConnection(TCPIPConnection* pConn)
    {   
        std::unique_lock<std::mutex> lock(m_mutex);
        auto it = m_connections.find(pConn);
        if (it != m_connections.end())
        {
            m_connections.erase(it);

            // clear write buffer from pConn
            auto it2 = m_writeBuffer.find(pConn);
            if (it2 != m_writeBuffer.end())
                m_writeBuffer.erase(it2);

            m_wakeupSignal.Set();

            lock.unlock();
            OnConnectionRemoved(pConn);
        }
    }

    void TCPIPConnectionHandler::AddServer(TCPIPServerSocket* pServer)
    {
        if (!pServer)
            return;
        
        std::unique_lock<std::mutex> lock(m_mutex);
        auto it = std::find(m_servers.begin(), m_servers.end(), pServer);
        if (it == m_servers.end())
        {
            m_servers.push_back(pServer);
            m_wakeupSignal.Set();
        }
    }

    void TCPIPConnectionHandler::RemoveServer(TCPIPServerSocket* pServer)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        auto it = std::find(m_servers.begin(), m_servers.end(), pServer);
        if (it != m_servers.end())
        {
            m_servers.erase(it);
            m_wakeupSignal.Set();
        }
    }

    void TCPIPConnectionHandler::WriteTo(TCPIPConnection* pConn, SocketDataPackage msg)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (m_connections.count(pConn) != 0)
        {
            m_writeBuffer[pConn].push(std::move(msg));
            m_wakeupSignal.Set();
        }
    }

    void TCPIPConnectionHandler::RunHandler()
    {
        std::map<TCPIPConnection*, std::tuple<SocketDataPackage, std::uint64_t>> readBuffer;
        std::set<TCPIPConnection*> killed;

        std::vector<TCPIPServerSocket*> checkAcceptable;
        std::vector<TCPIPConnection*> checkReadable;
        std::vector<TCPIPConnection*> checkWritable;

        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_keepRunning = true;
        }

        while (true)
        { 
            checkAcceptable.clear();
            checkReadable.clear();
            checkWritable.clear();
            killed.clear();
 
            // prepare lists of connections to be checked
            {
                std::unique_lock<std::mutex> lock(m_mutex);

                if (!m_keepRunning)
                    break;
                m_wakeupSignal.Reset();

                checkAcceptable.reserve(m_servers.size());
                checkReadable.reserve(m_connections.size());
                checkWritable.reserve(m_writeBuffer.size());
                std::copy(m_servers.begin(), m_servers.end(), std::back_inserter(checkAcceptable));
                std::copy(m_connections.begin(), m_connections.end(), std::back_inserter(checkReadable));
                std::transform(m_writeBuffer.begin(), m_writeBuffer.end(), 
                    std::back_inserter(checkWritable), [](auto el) { return el.first; });
            }

            // block until either accept, read or write can be performed
            auto [ac, re, wr, ex] = TCPIPServerSocket::Select(checkAcceptable, checkReadable, checkWritable, {}, m_wakeupSignal);
            
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if (!m_keepRunning)
                    break;
            }

            // handle new connections
            for (TCPIPServerSocket* pServer: ac)
                this->OnConnectionAcceptale(pServer);

            // validate read buffer
            for (auto it = readBuffer.begin(); it != readBuffer.end();)
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if (m_connections.count(it->first) == 0)
                    it = readBuffer.erase(it);
                else
                    it = std::next(it);
            }

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
                    
                    alreadyRead += cnt; // update read buffer (reference)

                    if (alreadyRead == buffer.GetSize())
                    {
                        SocketDataPackage msg = std::move(buffer);
                        readBuffer.erase(readBuffer.find(pConn));
                        OnConnectionMessage(pConn, std::move(msg)); 
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
                pConn->Close();
                this->OnConnectionClosed(pConn);
                this->RemoveConnection(pConn);
                
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

        // cleanup
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            while (!m_connections.empty())
            {
                TCPIPConnection* pConn = *m_connections.begin();
                pConn->Close();

                lock.unlock();
                this->RemoveConnection(pConn);
                lock.lock();
            }

            m_writeBuffer.clear();
        }
    }

    void TCPIPConnectionHandler::StopHandler()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_keepRunning = false;
        m_wakeupSignal.Set();
    }


    //
    // TCPIPServer
    //

    TCPIPServer::TCPIPServer() : m_pServer(nullptr) { }
    
    TCPIPServer::~TCPIPServer() { }

    bool TCPIPServer::Run(short port)
    {
        if (m_pServer)
            return false;

        // start socket server
        m_pServer = new TCPIPServerSocket();
        if (!m_pServer->Bind(8000))
            return false;
        if (!m_pServer->Listen(5))
            return false;

        this->AddServer(m_pServer);
        this->RunHandler();
        this->RemoveServer(m_pServer);

        delete m_pServer;
        m_pServer = nullptr;

        return true;
    }

    void TCPIPServer::Stop()
    {
        StopHandler();
    }

    void TCPIPServer::WriteTo(std::size_t id, SocketDataPackage msg)
    {
        TCPIPConnectionHandler::WriteTo(reinterpret_cast<TCPIPConnection*>(id), std::move(msg));
    }

    void TCPIPServer::OnConnectionAcceptale(TCPIPServerSocket* pServer) 
    {
        auto [fd, ip] = pServer->Accept();
        TCPIPConnection* pConn = new TCPIPConnection(fd);
        if (pConn->IsValid())
        {
            AddConnection(pConn);
            OnClientConnected(reinterpret_cast<std::size_t>(pConn), ip);
        }
    }

    void TCPIPServer::OnConnectionClosed(TCPIPConnection* pConn) 
    {
        OnClientDisconnected(reinterpret_cast<std::size_t>(pConn));
    }

    void TCPIPServer::OnConnectionRemoved(TCPIPConnection* pConn) 
    {
        if (pConn) delete pConn;
    }

    void TCPIPServer::OnConnectionMessage(TCPIPConnection* pConn, SocketDataPackage msg) 
    {
        OnMessageReceived(reinterpret_cast<std::size_t>(pConn), std::move(msg));
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
        
        AddConnection(pConn);
        return reinterpret_cast<std::size_t>(pConn);
    }
    
    std::size_t TCPIPMultiClient::ConnectHostname(const std::string& hostname, short port)
    {
        return Connect(TCPIPClientSocket::GetHostByName(hostname), port);
    }

    void TCPIPMultiClient::Run()
    {
        RunHandler();
    }

    void TCPIPMultiClient::Stop()
    {
        StopHandler();
    }

    void TCPIPMultiClient::WriteTo(std::size_t id, SocketDataPackage data)
    {
        TCPIPConnectionHandler::WriteTo(reinterpret_cast<TCPIPConnection*>(id), std::move(data));
    }

    void TCPIPMultiClient::OnConnectionAcceptale(TCPIPServerSocket* pServer) { }
    
    void TCPIPMultiClient::OnConnectionClosed(TCPIPConnection* pConn) 
    {
        OnClientDisconnected(reinterpret_cast<std::size_t>(pConn));
    }

    void TCPIPMultiClient::OnConnectionRemoved(TCPIPConnection* pConn) 
    { 
        if (pConn) delete pConn;
    }
    
    void TCPIPMultiClient::OnConnectionMessage(TCPIPConnection* pConn, SocketDataPackage msg) 
    {
        OnMessageReceived(reinterpret_cast<std::size_t>(pConn), std::move(msg));
    }

}
