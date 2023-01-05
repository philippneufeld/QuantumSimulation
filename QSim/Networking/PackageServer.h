// Philipp Neufeld, 2021-2022

#ifndef QSim_Networking_PackageServer_H_
#define QSim_Networking_PackageServer_H_

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <array>
#include <mutex>
#include <queue>

#include "TCPIP.h"
#include "../Util/UUID.h"

namespace QSim
{

    class DataPackagePayload
    {
    public:
        DataPackagePayload();
        DataPackagePayload(std::uint64_t size);
        virtual ~DataPackagePayload();

        DataPackagePayload(const DataPackagePayload& rhs);
        DataPackagePayload(DataPackagePayload&& rhs);

        DataPackagePayload& operator=(const DataPackagePayload& rhs);
        DataPackagePayload& operator=(DataPackagePayload&& rhs);

        bool Allocate(std::uint64_t size);

        // getter
        std::uint64_t GetSize() const { return m_size; };
        std::uint8_t* GetData() const { return m_pData; };

    private:
        std::uint64_t m_size;
        std::uint8_t* m_pData;
    };

    class NetworkDataPackage : public DataPackagePayload
    {
        constexpr static std::uint64_t s_protocolId = 0x13B5F6A39464CF1A;

        constexpr static std::size_t s_pidOff = 0;
        constexpr static std::size_t s_sizeOff = 8;
        constexpr static std::size_t s_msgOff = 16;
        constexpr static std::size_t s_topicOff = 20;
    public:
        using Header_t = std::array<std::uint8_t, 36>;

        NetworkDataPackage();
        NetworkDataPackage(std::uint64_t size);
        NetworkDataPackage(const Header_t& header);

        NetworkDataPackage(const NetworkDataPackage& rhs);
        NetworkDataPackage(NetworkDataPackage&& rhs);
        NetworkDataPackage(const DataPackagePayload& rhs);
        NetworkDataPackage(DataPackagePayload&& rhs);

        NetworkDataPackage& operator=(const NetworkDataPackage& rhs);
        NetworkDataPackage& operator=(NetworkDataPackage&& rhs);
        NetworkDataPackage& operator=(const DataPackagePayload& rhs);
        NetworkDataPackage& operator=(DataPackagePayload&& rhs);

        // function to change package state
        void SetMessageId(std::uint32_t msgId) { m_messageId = msgId; }
        void SetTopic(const UUIDv4& topic) { m_topic = topic; }

        // getter
        std::uint8_t GetMessageId() const { return m_messageId; }
        UUIDv4 GetTopic() const { return m_topic; }
        
        // header functions
        static bool IsValidHeader(const Header_t& header);
        static std::uint64_t GetSizeFromHeader(const Header_t& header);
        static std::uint32_t GetMessageIdFromHeader(const Header_t& header);
        static UUIDv4 GetTopicFromHeader(const Header_t& header);
        static Header_t GenerateHeader(std::uint64_t size, std::uint32_t msgId, const UUIDv4& topic);
        Header_t GetHeader() const;

        static NetworkDataPackage CreateEmptyPackage(std::uint32_t msgId);

    private:
        std::uint32_t m_messageId = 0;
        UUIDv4 m_topic;
    };

    class PackageConnectionHandler
    {
    public:
        PackageConnectionHandler();
        virtual ~PackageConnectionHandler();

        PackageConnectionHandler(const PackageConnectionHandler&) = delete;
        PackageConnectionHandler(PackageConnectionHandler&& rhs) = delete;

        PackageConnectionHandler& operator=(const PackageConnectionHandler&) = delete;
        PackageConnectionHandler& operator=(PackageConnectionHandler&&) = delete;

        void AddConnection(TCPIPConnection* pConn);
        void RemoveConnection(TCPIPConnection* pConn);
        void AddServer(TCPIPServerSocket* pServer);
        void RemoveServer(TCPIPServerSocket* pServer);

        void Broadcast(NetworkDataPackage msg);
        void WriteTo(TCPIPConnection* pConn, NetworkDataPackage msg);

        void RunHandler();
        void StopHandler();

        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) = 0;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) = 0;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) = 0;
        virtual void OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) = 0;

    private:
        std::mutex m_mutex;
        bool m_keepRunning;
        bool m_cleanup;

        IOEvent m_wakeupSignal;
        std::list<TCPIPServerSocket*> m_servers;
        std::set<TCPIPConnection*> m_connections;
        std::map<TCPIPConnection*, std::queue<std::tuple<NetworkDataPackage, std::uint64_t>>> m_writeBuffer;
    };

    class ConnectionUUIDMap
    {
    public:
        UUIDv4 GetIdFromConnection(const TCPIPConnection* pConn) const;
        TCPIPConnection* GetConnectionFromId(UUIDv4 id) const;

        UUIDv4 Insert(TCPIPConnection* pConn);
        void Remove(UUIDv4 id);
        void Remove(const TCPIPConnection* pConn);

    private:
        std::map<const TCPIPConnection*, UUIDv4> m_fromConn;
        std::map<UUIDv4, TCPIPConnection*> m_fromUuid;
    };

    class PackageServer : private PackageConnectionHandler
    {
    public:
        PackageServer();
        virtual ~PackageServer();

        bool Run(short port);
        void Stop();

        void Broadcast(NetworkDataPackage msg);
        void WriteTo(UUIDv4 id, NetworkDataPackage msg);

        void Disconnect(UUIDv4 id);

        // callbacks
        virtual void OnClientConnected(UUIDv4 id, const std::string& ip) { }
        virtual void OnClientDisconnected(UUIDv4 id) {}
        virtual void OnMessageReceived(UUIDv4 id, NetworkDataPackage data) { }

    private:
        // connection handler callbacks
        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) override;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) override;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) override;
        virtual void OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) override;

    private:
        TCPIPServerSocket* m_pServer;
        ConnectionUUIDMap m_ids;
    };

    class PackageClient : private PackageConnectionHandler
    {
    public:
        PackageClient();

        UUIDv4 Connect(const std::string& ip, short port);
        UUIDv4 ConnectHostname(const std::string& hostname, short port);

        virtual void OnClientDisconnected(UUIDv4 id) {}
        virtual void OnMessageReceived(UUIDv4 id, NetworkDataPackage data) { }

        void Run();
        void Stop();

        void WriteTo(UUIDv4 id, NetworkDataPackage data);
        void Disconnect(UUIDv4 id);

    private:
        // connection handler callbacks
        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) override;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) override;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) override;
        virtual void OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) override;

    private:
        ConnectionUUIDMap m_ids;
    };

}

#endif
