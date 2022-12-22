// Philipp Neufeld, 2021-2022

#include <algorithm>
#include <cstdlib>
#include <sstream>

#include <arpa/inet.h> // htonl, ntohl

#include "PackageServer.h"

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
    // NetworkDataPackage
    //

    NetworkDataPackage::NetworkDataPackage() 
        : m_size(0), m_pData(nullptr), m_status(NetworkDataPackageStatus_OK), m_messageId(0) {}

    NetworkDataPackage::NetworkDataPackage(std::uint64_t size) 
        : NetworkDataPackage() 
    {
        Allocate(size);
    }

    NetworkDataPackage::NetworkDataPackage(const Header_t& header)
        : NetworkDataPackage()
    {
        if (IsValidHeader(header))
        {
            m_status = GetStatusFromHeader(header);
            m_messageId = GetMessageIdFromHeader(header);
            Allocate(GetSizeFromHeader(header));
        }
    }
    
    NetworkDataPackage::~NetworkDataPackage() 
    { 
        Allocate(0);
    }

    NetworkDataPackage::NetworkDataPackage(const NetworkDataPackage& rhs)
        : NetworkDataPackage()
    {
        m_status = rhs.m_status;
        m_messageId = rhs.m_messageId;
        if (rhs.m_size > 0)
        {
            Allocate(rhs.m_size);
            std::copy(rhs.m_pData, rhs.m_pData+rhs.m_size, m_pData);
        }
    }
    
    NetworkDataPackage::NetworkDataPackage(NetworkDataPackage&& rhs)
        : m_pData(rhs.m_pData), m_size(rhs.m_size), m_status(rhs.m_status), m_messageId(rhs.m_messageId)
    {
        rhs.m_pData = nullptr;
        rhs.m_size = 0;
    }

    NetworkDataPackage& NetworkDataPackage::operator=(const NetworkDataPackage& rhs)
    {
        NetworkDataPackage tmp(rhs);
        std::swap(*this, tmp);
        return *this;
    }

    NetworkDataPackage& NetworkDataPackage::operator=(NetworkDataPackage&& rhs)
    {
        std::swap(m_pData, rhs.m_pData);
        std::swap(m_size, rhs.m_size);
        return *this;
    }

    NetworkDataPackage::operator bool() const
    {
        return m_status == NetworkDataPackageStatus_OK;
    }
    
    bool NetworkDataPackage::Allocate(std::uint64_t size)
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

    bool NetworkDataPackage::IsValidHeader(const Header_t& header)
    {
        std::uint64_t npid = 0;
        std::copy_n(header.data(), 7, reinterpret_cast<std::uint8_t*>(&npid) + 1);
        std::uint64_t pid = ntohll(npid);
        return (pid == s_protocolId);
    }

    std::uint8_t NetworkDataPackage::GetStatusFromHeader(const Header_t& header)
    {
        return header[7];
    }
    
    std::uint32_t NetworkDataPackage::GetMessageIdFromHeader(const Header_t& header)
    {
        std::uint32_t nmid = 0;
        std::copy_n(header.data() + 16, sizeof(nmid), reinterpret_cast<std::uint8_t*>(&nmid));
        return ntohl(nmid);
    }
    
    std::uint64_t NetworkDataPackage::GetSizeFromHeader(const Header_t& header)
    {
        std::uint64_t nsize = 0;
        std::copy_n(header.data() + 8, sizeof(nsize), reinterpret_cast<std::uint8_t*>(&nsize));
        return ntohll(nsize);
    }

    NetworkDataPackage::Header_t NetworkDataPackage::GenerateHeader(std::uint64_t size, std::uint8_t status, std::uint32_t msgId)
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

    NetworkDataPackage::Header_t NetworkDataPackage::GetHeader() const
    {
        return GenerateHeader(m_size, m_status, m_messageId);
    }
    
    NetworkDataPackage NetworkDataPackage::CreateError(std::uint8_t status)
    {
        NetworkDataPackage package;
        package.SetStatus(status);
        return package;
    }


    //
    // PackageConnectionHandler
    //

    PackageConnectionHandler::PackageConnectionHandler() { }
    
    PackageConnectionHandler::~PackageConnectionHandler() 
    {
        while (!m_connections.empty())
            RemoveConnection(*m_connections.begin());

        while (!m_servers.empty())
            RemoveServer(m_servers.front());
    }

    void PackageConnectionHandler::AddConnection(TCPIPConnection* pConn)
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

    void PackageConnectionHandler::RemoveConnection(TCPIPConnection* pConn)
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

    void PackageConnectionHandler::AddServer(TCPIPServerSocket* pServer)
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

    void PackageConnectionHandler::RemoveServer(TCPIPServerSocket* pServer)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        auto it = std::find(m_servers.begin(), m_servers.end(), pServer);
        if (it != m_servers.end())
        {
            m_servers.erase(it);
            m_wakeupSignal.Set();
        }
    }

    void PackageConnectionHandler::WriteTo(TCPIPConnection* pConn, NetworkDataPackage msg)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        if (m_connections.count(pConn) != 0)
        {
            m_writeBuffer[pConn].push(std::move(msg));
            m_wakeupSignal.Set();
        }
    }

    void PackageConnectionHandler::RunHandler()
    {
        std::map<TCPIPConnection*, std::tuple<NetworkDataPackage, std::uint64_t>> readBuffer;
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
                    NetworkDataPackage::Header_t header;
                    if (pConn->Recv(header.data(), sizeof(header)) != sizeof(header))
                    {
                        killed.insert(pConn);
                        continue;
                    }

                    if (!NetworkDataPackage::IsValidHeader(header))
                    {
                        killed.insert(pConn);
                        continue;
                    }
                    
                    readBuffer[pConn] = std::make_tuple(NetworkDataPackage(header), std::uint64_t(0));                    
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
                        NetworkDataPackage msg = std::move(buffer);
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
                    NetworkDataPackage package = std::move(outputQueue.front());
                    outputQueue.pop();
                    lock.unlock();

                    // send package header
                    NetworkDataPackage::Header_t header = package.GetHeader();
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

    void PackageConnectionHandler::StopHandler()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_keepRunning = false;
        m_wakeupSignal.Set();
    }


    //
    // PackageServer
    //

    PackageServer::PackageServer() : m_pServer(nullptr) { }
    
    PackageServer::~PackageServer() { }

    bool PackageServer::Run(short port)
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

    void PackageServer::Stop()
    {
        StopHandler();
    }

    void PackageServer::WriteTo(std::size_t id, NetworkDataPackage msg)
    {
        PackageConnectionHandler::WriteTo(reinterpret_cast<TCPIPConnection*>(id), std::move(msg));
    }

    void PackageServer::OnConnectionAcceptale(TCPIPServerSocket* pServer) 
    {
        auto [fd, ip] = pServer->Accept();
        TCPIPConnection* pConn = new TCPIPConnection(fd);
        if (pConn->IsValid())
        {
            AddConnection(pConn);
            OnClientConnected(reinterpret_cast<std::size_t>(pConn), ip);
        }
    }

    void PackageServer::OnConnectionClosed(TCPIPConnection* pConn) 
    {
        OnClientDisconnected(reinterpret_cast<std::size_t>(pConn));
    }

    void PackageServer::OnConnectionRemoved(TCPIPConnection* pConn) 
    {
        if (pConn) delete pConn;
    }

    void PackageServer::OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) 
    {
        OnMessageReceived(reinterpret_cast<std::size_t>(pConn), std::move(msg));
    }

    //
    // PackageMultiClient
    //

    PackageMultiClient::PackageMultiClient() { }

    std::size_t PackageMultiClient::Connect(const std::string& ip, short port)
    {
        auto pConn = new TCPIPClientSocket();
        if (!pConn->Connect(ip, port))
            return 0;
        
        AddConnection(pConn);
        return reinterpret_cast<std::size_t>(pConn);
    }
    
    std::size_t PackageMultiClient::ConnectHostname(const std::string& hostname, short port)
    {
        return Connect(TCPIPClientSocket::GetHostByName(hostname), port);
    }

    void PackageMultiClient::Run()
    {
        RunHandler();
    }

    void PackageMultiClient::Stop()
    {
        StopHandler();
    }

    void PackageMultiClient::WriteTo(std::size_t id, NetworkDataPackage data)
    {
        PackageConnectionHandler::WriteTo(reinterpret_cast<TCPIPConnection*>(id), std::move(data));
    }

    void PackageMultiClient::OnConnectionAcceptale(TCPIPServerSocket* pServer) { }
    
    void PackageMultiClient::OnConnectionClosed(TCPIPConnection* pConn) 
    {
        OnClientDisconnected(reinterpret_cast<std::size_t>(pConn));
    }

    void PackageMultiClient::OnConnectionRemoved(TCPIPConnection* pConn) 
    { 
        if (pConn) delete pConn;
    }
    
    void PackageMultiClient::OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) 
    {
        OnMessageReceived(reinterpret_cast<std::size_t>(pConn), std::move(msg));
    }

}
