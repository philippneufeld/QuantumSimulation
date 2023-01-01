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
#include <functional>
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
        virtual void OnClientConnected(UUIDv4 cid, const std::string& ip) override;
        virtual void OnClientDisconnected(UUIDv4 cid) override;
        virtual void OnMessageReceived(UUIDv4 cid, NetworkDataPackage data) override;

        UUIDv4 Reserve(UUIDv4 cid);
        std::size_t GetReservedCount() const;
        void CancelReservation(UUIDv4 cid, const UUIDv4& ticket);
        bool Run(UUIDv4 cid, const UUIDv4& ticket, NetworkDataPackage data);
        void BroadcastAvailability();

    private:
        ThreadPool m_pool;
        std::map<UUIDv4, std::set<UUIDv4>> m_tickets;
    };

    class ServerPool : private PackageClient
    {
    public:
        ServerPool() = default;

        using PackageClient::Run;
        using PackageClient::Stop;

        UUIDv4 ConnectWorker(const std::string& ip, short port);
        UUIDv4 ConnectWorkerHostname(const std::string& hostname, short port);
        std::size_t GetWorkerCount() const;

        std::pair<UUIDv4, std::future<DataPackagePayload>> Submit(DataPackagePayload task);
        std::pair<UUIDv4, std::future<DataPackagePayload>> Submit(std::function<DataPackagePayload()> task);
        virtual void OnTaskCompleted(UUIDv4 id) {};
        void WaitUntilFinished();

    private:
        virtual void OnClientDisconnected(UUIDv4 worker) override;
        virtual void OnMessageReceived(UUIDv4 worker, NetworkDataPackage data) override;
        
    private:
        mutable std::mutex m_mutex;
        std::condition_variable m_taskFinished;
        std::set<UUIDv4> m_workers;
        std::deque<std::tuple<
            UUIDv4, // task id
            std::function<DataPackagePayload()>, // task data generator
            std::promise<DataPackagePayload> // task result promise
        >> m_unscheduled;
        std::list<std::tuple<
            UUIDv4, // task id
            std::function<DataPackagePayload()>, // task data generator
            std::promise<DataPackagePayload>, // task result promise
            UUIDv4, // worker id
            UUIDv4 // ticket
        >> m_executing;
    };

}


#endif
