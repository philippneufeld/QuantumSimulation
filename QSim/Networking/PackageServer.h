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

namespace QSim
{

    enum NetworkDataPackageStatus
    {
        NetworkDataPackageStatus_OK = 0x00,
        NetworkDataPackageStatus_FAIL = 0x01,
        NetworkDataPackageStatus_INVALID_PACKAGE = 0x02,
        NetworkDataPackageStatus_IO_ERROR = 0x03,
    };

    class NetworkDataPackage
    {
        constexpr static std::uint64_t s_protocolId = 0x00A3F6A39464CF1A;
    public:
        using Header_t = std::array<std::uint8_t, 20>;

        NetworkDataPackage();
        NetworkDataPackage(std::uint64_t size);
        NetworkDataPackage(const Header_t& header);
        ~NetworkDataPackage();

        NetworkDataPackage(const NetworkDataPackage& rhs);
        NetworkDataPackage(NetworkDataPackage&& rhs);

        NetworkDataPackage& operator=(const NetworkDataPackage& rhs);
        NetworkDataPackage& operator=(NetworkDataPackage&& rhs);

        operator bool() const;

        // function to change package state
        bool Allocate(std::uint64_t size);
        void SetStatus(std::uint8_t status) { m_status = status; }
        void SetMessageId(std::uint8_t msgId) { m_messageId = msgId; }

        // getter
        std::uint8_t GetStatus() const { return m_status; }
        std::uint8_t GetMessageId() const { return m_messageId; }
        std::uint64_t GetSize() const { return m_size; };
        std::uint8_t* GetData() const { return m_pData; };

        // header functions
        static bool IsValidHeader(const Header_t& header);
        static std::uint8_t GetStatusFromHeader(const Header_t& header);
        static std::uint32_t GetMessageIdFromHeader(const Header_t& header);
        static std::uint64_t GetSizeFromHeader(const Header_t& header);
        static Header_t GenerateHeader(std::uint64_t size, std::uint8_t status, std::uint32_t msgId);
        Header_t GetHeader() const;

        // creation function
        static NetworkDataPackage CreateError(std::uint8_t status);

    private:
        std::uint8_t m_status;
        std::uint64_t m_size;
        std::uint32_t m_messageId;
        std::uint8_t* m_pData;
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
        std::map<TCPIPConnection*, std::queue<NetworkDataPackage>> m_writeBuffer;
    };

    class PackageServer : private PackageConnectionHandler
    {
    public:
        PackageServer();
        virtual ~PackageServer();

        bool Run(short port);
        void Stop();

        void WriteTo(std::size_t id, NetworkDataPackage msg);

        // callbacks
        virtual bool OnClientConnected(std::size_t id, const std::string& ip) { return true; }
        virtual void OnClientDisconnected(std::size_t id) {}
        virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) { }

    private:
        // connection handler callbacks
        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) override;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) override;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) override;
        virtual void OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) override;

    private:
        TCPIPServerSocket* m_pServer;
    };

    class PackageMultiClient : private PackageConnectionHandler
    {
    public:
        PackageMultiClient();

        std::size_t Connect(const std::string& ip, short port);
        std::size_t ConnectHostname(const std::string& hostname, short port);

        virtual void OnClientDisconnected(std::size_t id) {}
        virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) { }

        void Run();
        void Stop();

        void WriteTo(std::size_t id, NetworkDataPackage data);

    private:
        // connection handler callbacks
        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) override;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) override;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) override;
        virtual void OnConnectionMessage(TCPIPConnection* pConn, NetworkDataPackage msg) override;
    };

}

#endif
