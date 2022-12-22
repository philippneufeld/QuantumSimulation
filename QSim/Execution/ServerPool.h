// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>
#include <string>
#include <tuple>
#include <set>

#include "../Util/TCPIP.h"
#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

    // 32-bit message descriptor
    enum ServerPoolQuery
    {
        ServerPoolQuery_Reserve = 0x01,
        ServerPoolQuery_Post = 0x02,
    };

    class ServerPoolWorker : public TCPIPServer
    {
    public:
        ServerPoolWorker();
        ServerPoolWorker(std::size_t threadCnt);

        virtual void OnMessageReceived(std::size_t id, SocketDataPackage data) override;

        virtual SocketDataPackage DoWork(SocketDataPackage data) = 0;

    private:
        std::uint32_t TryReserve(std::size_t id, std::uint32_t cnt);
        bool TryRun(std::size_t id, SocketDataPackage data);

    private:
        ThreadPool m_pool;
        std::multiset<std::size_t> m_tickets;
    };

    class ServerPoolMaster
    {
        ServerPoolMaster();

    public:


    };

}


#endif
