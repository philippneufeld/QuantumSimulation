// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>
#include <string>
#include <tuple>
#include <list>
#include <map>
#include <set>
#include <future>
#include <deque>
#include <mutex>
#include <condition_variable>

#include "../Networking/PackageServer.h"
#include "../Util/UUID.h"
#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

    // 32-bit message descriptor
    enum ServerPoolMessage
    {
        ServerPool_InvalidRequest,
        ServerPool_ReserveRequest,
        ServerPool_CancelReservationRequest,
        ServerPool_NotReserved,
        ServerPool_Reserved,
        ServerPool_PostRequest,
        ServerPool_NotPosted,
        ServerPool_Posted,
        ServerPool_TaskCompleted,
        ServerPool_CapacityAvailable,
    };

    class ServerPoolWorker : private PackageServer
    {
    public:
        ServerPoolWorker();
        ServerPoolWorker(std::size_t threadCnt);

        using PackageServer::Run;
        using PackageServer::Stop;

        virtual DataPackagePayload DoWork(DataPackagePayload data) = 0;

    private:
        virtual void OnClientConnected(std::size_t id, const std::string& ip) override;
        virtual void OnClientDisconnected(std::size_t id) override;
        virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) override;

        UUIDv4 Reserve(std::size_t id);
        void CancelReservation(std::size_t id, const UUIDv4& ticket);
        bool Run(std::size_t id, const UUIDv4& ticket, NetworkDataPackage data);
        void BroadcastAvailability();

    private:
        ThreadPool m_pool;
        std::map<std::size_t, std::set<UUIDv4>> m_tickets;
    };

    class ServerPool : private PackageClient
    {
    public:
        ServerPool() = default;

        using PackageClient::Run;
        using PackageClient::Stop;

        std::size_t ConnectWorker(const std::string& ip, short port);
        std::size_t ConnectWorkerHostname(const std::string& hostname, short port);
        std::size_t GetWorkerCount() const;

        std::future<DataPackagePayload> Submit(DataPackagePayload task);
        void WaitUntilFinished();

    private:
        virtual void OnClientDisconnected(std::size_t worker) override;
        virtual void OnMessageReceived(std::size_t worker, NetworkDataPackage data) override;
        
    private:
        mutable std::mutex m_mutex;
        std::condition_variable m_taskFinished;
        std::set<std::size_t> m_workers;
        std::deque<std::tuple<
            DataPackagePayload, 
            std::promise<DataPackagePayload>
        >> m_unscheduled;
        std::list<std::tuple<
            DataPackagePayload, 
            std::promise<DataPackagePayload>, 
            std::size_t, 
            UUIDv4
        >> m_executing;
    };

}


#endif
