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

#include <iostream>

namespace QSim
{

    ServerPoolWorker::ServerPoolWorker() : m_pool() { }

    ServerPoolWorker::ServerPoolWorker(std::size_t threadCnt) : m_pool(threadCnt) { }

    void ServerPoolWorker::OnMessageReceived(std::size_t id, NetworkDataPackage data)
    {
        NetworkDataPackage response;
        response.SetMessageId(ServerPool_InvalidRequest);

        switch (data.GetMessageId())
        {
        case ServerPool_ReserveRequest:
        {
            std::uint32_t reserveCnt = *reinterpret_cast<std::uint32_t*>(data.GetData());
            reserveCnt = TryReserve(id, reserveCnt);

            response.Allocate(8);
            response.SetMessageId(ServerPool_Reserved);
            *reinterpret_cast<std::uint32_t*>(response.GetData()) = reserveCnt;
        } 
        break;
        
        case ServerPool_PostRequest:
        {
            if (TryRun(id, std::move(data)))
                response.SetMessageId(ServerPool_Posted);
            else
                response.SetMessageId(ServerPool_NotPosted);
        } 
        break;
        
        default:
        {
            response.SetMessageId(ServerPool_InvalidRequest);
        } 
        break;
        }

        WriteTo(id, std::move(response));
    }

    std::uint32_t ServerPoolWorker::TryReserve(std::size_t id, std::uint32_t cnt)
    {
        std::size_t available = m_pool.GetAvailableThreads();
        std::size_t reserved = m_tickets.size();
        if (available <= reserved)
            return 0;

        available -= reserved;

        if (cnt > available)
            cnt = available;

        for (std::uint32_t i = 0; i < cnt; i++)
            m_tickets.insert(id);
        
        return cnt;
    }

    bool ServerPoolWorker::TryRun(std::size_t id, NetworkDataPackage data)
    {
        // check for ticket
        auto it = m_tickets.find(id);
        if (it == m_tickets.end())
        {
            // no ticket
            TryReserve(id, 1);
            it = m_tickets.find(id);
            if (it == m_tickets.end())
                return false;
        }

        // remove ticket
        m_tickets.erase(it);

        m_pool.Submit([this, id, data=std::move(data)]()
        { 
            NetworkDataPackage result = DoWork(std::move(data));
            result.SetMessageId(ServerPool_TaskCompleted);
            WriteTo(id, std::move(result));
        });

        return true;
    }

}
