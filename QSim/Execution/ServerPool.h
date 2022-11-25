// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>
#include <string>
#include <tuple>

#include "../Util/TCPIP.h"
#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

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
