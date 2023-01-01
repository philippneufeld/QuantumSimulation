// Philipp Neufeld, 2021-2022

#include <algorithm>
#include <cstdlib>
#include <sstream>

#include <unistd.h> // close
#include <sys/socket.h> // socket functions
#include <netinet/in.h> // sockaddr_in
#include <arpa/inet.h> // inet_ntoa
#include <netdb.h> // gethostbyname

#include "ServerPool.h"

namespace QSim
{

    ServerPoolWorker::ServerPoolWorker() : m_pool() { }

    ServerPoolWorker::ServerPoolWorker(std::size_t threadCnt) : m_pool(threadCnt) { }

    void ServerPoolWorker::OnClientConnected(UUIDv4 cid, const std::string& ip)
    {
        std::size_t available = m_pool.GetAvailableThreads();
        std::size_t reserved = GetReservedCount();
        if (available > reserved)
        {
            NetworkDataPackage msg;
            msg.SetMessageId(ServerPool_CapacityAvailable);
            WriteTo(cid, std::move(msg));
        }
    }

    void ServerPoolWorker::OnClientDisconnected(UUIDv4 cid)
    {
        // remove outstanding tickets
        auto it = m_tickets.find(cid);
        for (; it!=m_tickets.end(); it = m_tickets.find(cid))
            m_tickets.erase(it);
    }

    void ServerPoolWorker::OnMessageReceived(UUIDv4 cid, NetworkDataPackage data)
    {
        NetworkDataPackage response;
        response.SetMessageId(ServerPool_InvalidRequest);
        response.SetTopic(data.GetTopic());

        switch (data.GetMessageId())
        {
        case ServerPool_ReserveRequest:
        {
            auto uuid = Reserve(cid);
            bool reserved = !uuid.IsNullUUID();

            response.SetMessageId(reserved ? ServerPool_Reserved : ServerPool_NotReserved);

            if (reserved)
            {
                response.Allocate(sizeof(uuid));
                uuid.StoreToBufferLE(response.GetData());
            }          
        } 
        break;
        
        case ServerPool_PostRequest:
        {
            auto ticket = data.GetTopic();
            if (Run(cid, ticket, std::move(data)))
                response.SetMessageId(ServerPool_Posted);
            else
                response.SetMessageId(ServerPool_NotPosted);
        } 
        break;

        case ServerPool_CancelReservationRequest:
        {
            auto ticket = data.GetTopic();
            CancelReservation(cid, ticket);
        }
        break;
        
        default:
        {
            response.SetMessageId(ServerPool_InvalidRequest);
        } 
        break;
        }

        WriteTo(cid, std::move(response));
    }

    UUIDv4 ServerPoolWorker::Reserve(UUIDv4 cid)
    {
        std::size_t available = m_pool.GetAvailableThreads();
        std::size_t reserved = GetReservedCount();
        if (available <= reserved)
            return UUIDv4::NullUUID();

        UUIDv4 uuid;
        m_tickets[cid].insert(uuid);
        return uuid;
    }

    std::size_t ServerPoolWorker::GetReservedCount() const
    {
        std::size_t reserved = 0;
        for (auto& tickets: m_tickets) 
            reserved += tickets.second.size();
        return reserved;
    }

    void ServerPoolWorker::CancelReservation(UUIDv4 cid, const UUIDv4& ticket)
    {
        auto it = m_tickets[cid].find(ticket);
        if (it != m_tickets[cid].end())
            m_tickets[cid].erase(it);

        BroadcastAvailability();
    }

    bool ServerPoolWorker::Run(UUIDv4 cid, const UUIDv4& ticket, NetworkDataPackage data)
    {
        // check for ticket
        auto it = m_tickets[cid].find(ticket);
        if (it == m_tickets[cid].end())
            return false;

        // remove ticket
        m_tickets[cid].erase(it);

        m_pool.Submit([this, cid, ticket, data=std::move(data)]()
        { 
            NetworkDataPackage result = DoWork(std::move(data));
            result.SetMessageId(ServerPool_TaskCompleted);
            result.SetTopic(ticket);
            WriteTo(cid, std::move(result));
            BroadcastAvailability();
        });

        return true;
    }

    void ServerPoolWorker::BroadcastAvailability()
    {
        std::size_t available = m_pool.GetAvailableThreads();
        std::size_t reserved = GetReservedCount();
        if (available > reserved)
        {
            NetworkDataPackage msg;
            msg.SetMessageId(ServerPool_CapacityAvailable);
            Broadcast(std::move(msg));
        }
    }


    //
    // ServerPool
    //

    UUIDv4 ServerPool::ConnectWorker(const std::string& ip, short port)
    {
        auto worker = PackageClient::Connect(ip, port);
        
        std::unique_lock<std::mutex> lock(m_mutex);
        if (worker) m_workers.insert(worker);
        lock.unlock();
        
        return worker;
    }

    UUIDv4 ServerPool::ConnectWorkerHostname(const std::string& hostname, short port)
    {
        auto worker = PackageClient::ConnectHostname(hostname, port);
        
        std::unique_lock<std::mutex> lock(m_mutex);
        if (worker) m_workers.insert(worker);
        lock.unlock();

        return worker;
    }

    std::size_t ServerPool::GetWorkerCount() const
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        return m_workers.size();
    }

    std::pair<UUIDv4, std::future<DataPackagePayload>> ServerPool::Submit(DataPackagePayload task)
    {
        return Submit([task=std::move(task)]{ return task; });
    }

    std::pair<UUIDv4, std::future<DataPackagePayload>> ServerPool::Submit(std::function<DataPackagePayload()> task)
    {
        std::promise<DataPackagePayload> promise;
        std::future<DataPackagePayload> future = promise.get_future();

        UUIDv4 id;

        std::unique_lock<std::mutex> lock(m_mutex);
        m_unscheduled.push_back(std::make_tuple(id, std::move(task), std::move(promise)));
        lock.unlock();

        // send reserve request to all workers
        for (auto worker: m_workers)
        {
            NetworkDataPackage request;
            request.SetMessageId(ServerPool_ReserveRequest);
            WriteTo(worker, request);
        }

        return std::make_pair(id, std::move(future));
    }

    void ServerPool::WaitUntilFinished()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_taskFinished.wait(lock, [this](){ return m_unscheduled.empty() && m_executing.empty(); });
    }

    void ServerPool::OnClientDisconnected(UUIDv4 worker)
    {
        // check if there were uncompleted tasks that were scheduled with the worker
        std::unique_lock<std::mutex> lock(m_mutex);
        m_workers.erase(worker);
        for (auto it=m_executing.begin(); it!=m_executing.end(); it++)
        {
            auto& [id, task, prom, w, t] = *it;
            if (w == worker)
            {
                m_unscheduled.push_front(std::make_tuple(id, std::move(task), std::move(prom)));
                it = m_executing.erase(it);
            }
        }
    }

    void ServerPool::OnMessageReceived(UUIDv4 worker, NetworkDataPackage data)
    {
        switch (data.GetMessageId())
        {
        case ServerPool_Reserved:
        {
            auto ticket = UUIDv4::LoadFromBufferLE(data.GetData());

            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_unscheduled.size() == 0)
            {
                NetworkDataPackage request;
                request.SetMessageId(ServerPool_CancelReservationRequest);
                request.SetTopic(ticket);
                WriteTo(worker, request);
            }
            else
            {
                auto [id, task, prom] = std::move(m_unscheduled.front());
                m_unscheduled.pop_front();
                
                NetworkDataPackage request = std::invoke(task);
                request.SetMessageId(ServerPool_PostRequest);
                request.SetTopic(ticket);

                m_executing.push_back(std::make_tuple(id, std::move(task), std::move(prom), worker, ticket));

                WriteTo(worker, std::move(request));
            }
        }
        break;

        case ServerPool_NotPosted:
        {
            auto ticket = data.GetTopic();
            std::unique_lock<std::mutex> lock(m_mutex);

            // search for task
            auto it=m_executing.begin();
            for (; it!=m_executing.end(); it++)
            {
                auto& [id, task, prom, w, t] = *it;
                if (w == worker && t == ticket)
                    break;
            }

            // enlist task back into the unscheduled list
            if (it != m_executing.end())
            {
                auto [id, task, prom, w, t] = std::move(*it);
                m_executing.erase(it);
                m_unscheduled.push_front(std::make_tuple(id, std::move(task), std::move(prom)));
            }                
        }
        break;

        case ServerPool_CapacityAvailable:
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            if (!m_unscheduled.empty())
            {
                NetworkDataPackage request;
                request.SetMessageId(ServerPool_ReserveRequest);
                WriteTo(worker, request);
            }
        }
        break;

        case ServerPool_TaskCompleted:
        {
            auto ticket = data.GetTopic();
            std::unique_lock<std::mutex> lock(m_mutex);

            // search for task
            auto it=m_executing.begin();
            for (; it!=m_executing.end(); it++)
            {
                auto& [id, task, prom, w, t] = *it;
                if (w == worker && t == ticket)
                {
                    prom.set_value(std::move(data));
                    m_executing.erase(it);
                    OnTaskCompleted(id);
                    break;
                }
            }

            lock.unlock();
            m_taskFinished.notify_all();
        }
        break;

        };
    }

}
