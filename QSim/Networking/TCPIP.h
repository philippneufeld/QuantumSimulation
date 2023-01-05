// Philipp Neufeld, 2021-2022

#ifndef QSim_Networking_TCPIP_H_
#define QSim_Networking_TCPIP_H_

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace QSim
{

    class IOEvent
    {
    public:
        IOEvent();
        virtual ~IOEvent();

        IOEvent(const IOEvent&) = delete;
        IOEvent(IOEvent&&) = default;

        IOEvent& operator=(const IOEvent&) = delete;
        IOEvent& operator=(IOEvent&&) = default;

        void Set();
        void Reset();
        int GetFileDescriptor() const { return m_fd[0]; }

    private:
        int m_fd[2];
    };

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

        std::make_signed_t<std::size_t> Send(const void* data, std::size_t n);
        std::make_signed_t<std::size_t> Recv(void* buffer, std::size_t buffSize);

        std::make_signed_t<std::size_t> SendNonBlock(const void* data, std::size_t n);
        std::make_signed_t<std::size_t> RecvNonBlock(void* buffer, std::size_t buffSize);
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
                   const std::vector<TCPIPConnection*>& exceptional,
                   IOEvent& wakeupSignal);
    };

    class TCPIPClientSocket : public TCPIPConnection
    {
    public:
        TCPIPClientSocket();
        bool Connect(const std::string& ip, short port);
        static std::string GetHostByName(const std::string& hostname);
    };

}

#endif
