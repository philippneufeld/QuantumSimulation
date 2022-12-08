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

    SocketDataPackage ServerPoolWorker::OnMessageReceived(std::size_t id, SocketDataPackage data)
    {
        if (data.GetSize() < 4)
            return SocketDataPackage();
        
        std::uint32_t desc = *reinterpret_cast<std::uint32_t*>(data.GetData());
        std::uint8_t* payload = (data.GetData() + 4);

        if (desc == ServerPoolQuery_Reserve)
        {
            std::uint32_t reserveCnt = *reinterpret_cast<std::uint32_t*>(payload);
            reserveCnt = TryReserve(id, reserveCnt);

            SocketDataPackage response(8);
            *reinterpret_cast<std::uint32_t*>(response.GetData()) = ServerPoolQuery_Reserve;
            *reinterpret_cast<std::uint32_t*>(response.GetData() + 4) = reserveCnt;
            return response;
        }

        return SocketDataPackage();
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

    bool ServerPoolWorker::TryRun(std::size_t id, SocketDataPackage data)
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

        m_pool.Submit([this, data=std::move(data)](){ return DoWork(data); });

        return true;
    }

}
