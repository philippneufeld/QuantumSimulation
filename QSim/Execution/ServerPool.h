// Philipp Neufeld, 2021-2022

#ifndef QSim_Execution_ServerPool_H_
#define QSim_Execution_ServerPool_H_


#include <cstdint>

#include "Progress.h"
#include "ThreadPool.h"

namespace QSim
{

    class SocketDataPackageBin
    {
    public:
        SocketDataPackageBin() : m_size(0), m_pData(nullptr) {}
        ~SocketDataPackageBin() { }

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
    public:


    };

}


#endif
