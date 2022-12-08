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
            Select(const std::list<TCPIPServerSocket*>& acceptable,
                   const std::list<TCPIPConnection*>& readable, 
                   const std::list<TCPIPConnection*>& writable, 
                   const std::list<TCPIPConnection*>& exceptional);
    };

    class TCPIPClientSocket : public TCPIPConnection
    {
    public:
        TCPIPClientSocket();
        bool Connect(const std::string& ip, short port);
        bool ConnectHostname(const std::string& hostname, short port);
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

    class TCPIPServer
    {
    public:
        TCPIPServer();
        virtual ~TCPIPServer();

        bool Run(short port);

        // callbacks
        virtual bool OnClientConnected(std::size_t id, const std::string& ip) { return true; }
        virtual void OnClientDisconnected(std::size_t id) {}
        virtual SocketDataPackage OnMessageReceived(std::size_t id, SocketDataPackage data) { return SocketDataPackage(); }

    protected:
        void Purge();
        void CloseClientConnection(TCPIPConnection* pClient);
        std::size_t GetIdFromConnection(TCPIPConnection* pClient) const;
        
    private:
        bool m_bRunning;
        std::list<TCPIPConnection*> m_pClients;
        
        // buffers
        std::map<TCPIPConnection*, std::tuple<SocketDataPackage, std::uint64_t>> m_readBuffer;
        std::map<TCPIPConnection*, std::list<SocketDataPackage>> m_writeBuffer;
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

}

#endif
