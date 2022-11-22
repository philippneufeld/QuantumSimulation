// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>
#include <string>
#include <tuple>

#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

    class SocketDataPackageBin
    {
    public:
        SocketDataPackageBin();
        ~SocketDataPackageBin();

        SocketDataPackageBin(const SocketDataPackageBin& rhs);
        SocketDataPackageBin(SocketDataPackageBin&& rhs);

        SocketDataPackageBin& operator=(const SocketDataPackageBin& rhs);
        SocketDataPackageBin& operator=(SocketDataPackageBin&& rhs);

        bool Allocate(std::uint64_t size);

        std::uint64_t GetSize() const { return m_size; };
        std::uint8_t* GetData() const { return m_pData; };
        
    private:
        std::uint64_t m_size;
        std::uint8_t* m_pData;
    };

    class SocketDataPackage
    {
    public:
        virtual SocketDataPackageBin Serialize() = 0;
        virtual void Deserialize(const SocketDataPackageBin& data) = 0;
    };


    class TCPIPSocket
    {
    protected:
        TCPIPSocket(int handle);
    public:
        TCPIPSocket();
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
        friend class TCPIPSocketServer;
    protected:
        TCPIPConnection(int handle);
    public:
        TCPIPConnection();

        std::size_t Send(const void* data, std::size_t n);
    };

    class TCPIPSocketServer : public TCPIPSocket
    {
    protected:
        TCPIPSocketServer(int handle);
    public:
        TCPIPSocketServer();
        
        bool Bind(short port);
        bool Listen(int queue_length);
        std::tuple<TCPIPConnection, std::string> Accept();
    };

    class TCPIPSocketClient : public TCPIPConnection
    {
    public:
        TCPIPSocketClient();
        bool Connect(const std::string& ip, short port);
        bool ConnectHostname(const std::string& hostname, short port);
    };




    class ServerPoolWorker
    {
    public:
        ServerPoolWorker();
        virtual ~ServerPoolWorker();

    private:
        ThreadPool m_pool;
    };

    class ServerPoolMaster
    {
        ServerPoolMaster();

    public:


    };

}


#endif
