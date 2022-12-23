// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>
#include <string>
#include <tuple>
#include <set>

#include "../Networking/PackageServer.h"
#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

    // 32-bit message descriptor
    enum ServerPoolMessage
    {
        ServerPool_InvalidRequest,
        ServerPool_ReserveRequest,
        ServerPool_Reserved,
        ServerPool_PostRequest,
        ServerPool_NotPosted,
        ServerPool_Posted,
        ServerPool_TaskCompleted,
    };

    class ServerPoolWorker : public PackageServer
    {
    public:
        ServerPoolWorker();
        ServerPoolWorker(std::size_t threadCnt);

        virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) override;

        virtual DataPackagePayload DoWork(DataPackagePayload data) = 0;

    private:
        std::uint32_t TryReserve(std::size_t id, std::uint32_t cnt);
        bool TryRun(std::size_t id, NetworkDataPackage data);

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
