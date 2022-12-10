// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_TCPIP_H_
#define QSim_Util_TCPIP_H_

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <array>
#include <functional>
#include <atomic>
#include <mutex>
#include <queue>

namespace QSim
{

    class TCPIPSocket
    {
    public:
        TCPIPSocket();
        TCPIPSocket(int handle);
        virtual ~TCPIPSocket();

        TCPIPSocket(const TCPIPSocket&) = delete;
        TCPIPSocket(TCPIPSocket&& rhs);

        TCPIPSocket& operator=(const TCPIPSocket&) = delete;
        TCPIPSocket& operator=(TCPIPSocket&& rhs);

        bool IsValid() const { return (m_handle > 0);}
        operator bool() const { return IsValid(); }
        int GetHandle() const { return m_handle; }

        void Close();

    private:
        int m_handle;
    };

    class TCPIPConnection : public TCPIPSocket
    {
    public:
        TCPIPConnection();
        TCPIPConnection(int handle);

        std::int64_t Send(const void* data, std::size_t n);
        std::int64_t Recv(void* buffer, std::size_t buffSize);
    };

    class TCPIPServerSocket : public TCPIPSocket
    {
    public:
        TCPIPServerSocket();
        
        bool Bind(short port);
        bool Listen(int queue_length);
        std::tuple<int, std::string> Accept();
        std::tuple<TCPIPConnection, std::string> AcceptC();

        static std::tuple<
            std::vector<TCPIPServerSocket*>, std::vector<TCPIPConnection*>, 
            std::vector<TCPIPConnection*>, std::vector<TCPIPConnection*>
        >
            Select(const std::vector<TCPIPServerSocket*>& acceptable,
                   const std::vector<TCPIPConnection*>& readable, 
                   const std::vector<TCPIPConnection*>& writable, 
                   const std::vector<TCPIPConnection*>& exceptional);
    };

    class TCPIPClientSocket : public TCPIPConnection
    {
    public:
        TCPIPClientSocket();
        bool Connect(const std::string& ip, short port);
        static std::string GetHostByName(const std::string& hostname);
    };


    enum SocketDataPackageStatus
    {
        SocketDataPackageStatus_OK = 0x00,
        SocketDataPackageStatus_FAIL = 0x01,
        SocketDataPackageStatus_INVALID_PACKAGE = 0x02,
        SocketDataPackageStatus_IO_ERROR = 0x03,
    };

    class SocketDataPackage
    {
        constexpr static std::uint64_t s_protocolId = 0x00A3F6A39464CF1A;
    public:
        using Header_t = std::array<std::uint8_t, 20>;

        SocketDataPackage();
        SocketDataPackage(std::uint64_t size);
        SocketDataPackage(const Header_t& header);
        ~SocketDataPackage();

        SocketDataPackage(const SocketDataPackage& rhs);
        SocketDataPackage(SocketDataPackage&& rhs);

        SocketDataPackage& operator=(const SocketDataPackage& rhs);
        SocketDataPackage& operator=(SocketDataPackage&& rhs);

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
        static SocketDataPackage CreateError(std::uint8_t status);

    private:
        std::uint8_t m_status;
        std::uint64_t m_size;
        std::uint32_t m_messageId;
        std::uint8_t* m_pData;
    };

    using UUID = std::size_t;

    class TCPIPConnectionHandler
    {
    public:
        TCPIPConnectionHandler();
        virtual ~TCPIPConnectionHandler();

        TCPIPConnectionHandler(const TCPIPConnectionHandler&) = delete;
        TCPIPConnectionHandler(TCPIPConnectionHandler&& rhs) = delete;

        TCPIPConnectionHandler& operator=(const TCPIPConnectionHandler&) = delete;
        TCPIPConnectionHandler& operator=(TCPIPConnectionHandler&&) = delete;

        void AddConnection(TCPIPConnection* pConn);
        void RemoveConnection(TCPIPConnection* pConn);
        void AddServer(TCPIPServerSocket* pServer);
        void RemoveServer(TCPIPServerSocket* pServer);

        void RunHandler();
        void StopHandler();

        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) = 0;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) = 0;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) = 0;
        virtual SocketDataPackage OnConnectionMessage(TCPIPConnection* pConn, SocketDataPackage msg) = 0;

    private:
        void DoRunHandler();

    private:
        std::mutex m_mutex;
        bool m_keepRunning;

        std::list<TCPIPServerSocket*> m_servers;
        std::list<TCPIPConnection*> m_connections;
        std::map<TCPIPConnection*, std::queue<SocketDataPackage>> m_writeBuffer;
    };

    class TCPIPServer : private TCPIPConnectionHandler
    {
    public:
        TCPIPServer();
        virtual ~TCPIPServer();

        bool Run(short port);

        // callbacks
        virtual bool OnClientConnected(std::size_t id, const std::string& ip) { return true; }
        virtual void OnClientDisconnected(std::size_t id) {}
        virtual SocketDataPackage OnMessageReceived(std::size_t id, SocketDataPackage data) { return SocketDataPackage(); }

    private:
        // connection handler callbacks
        virtual void OnConnectionAcceptale(TCPIPServerSocket* pServer) override;
        virtual void OnConnectionClosed(TCPIPConnection* pConn) override;
        virtual void OnConnectionRemoved(TCPIPConnection* pConn) override;
        virtual SocketDataPackage OnConnectionMessage(TCPIPConnection* pConn, SocketDataPackage msg) override;

        std::size_t GetIdFromConnection(TCPIPConnection* pConnection) const;

    private:
        TCPIPServerSocket* m_pServer;
    };

    class TCPIPClient : public TCPIPClientSocket
    {
        // disallow the use of base class send and receive
    protected:
        using TCPIPClientSocket::Recv;
        using TCPIPClientSocket::Send;

    public:
        SocketDataPackage Query(const void* data, std::uint64_t n, std::uint32_t msgId = 0);
    };

    class TCPIPMultiClient
    {
    public:
        TCPIPMultiClient();

        std::size_t Connect(const std::string& ip, short port);
        std::size_t ConnectHostname(const std::string& hostname, short port);

    private:
        std::size_t GetIdFromConnection(TCPIPConnection* pConnection) const;

    private:
        std::list<TCPIPConnection*> m_connections;
    };

}

#endif
