// Philipp Neufeld, 2021-2022

#include <iostream>
#include <chrono>
#include <list>
#include <Eigen/Dense>

#include <sys/select.h>

#include <QSim/Util/UUID.h>
#include <QSim/Util/Argparse.h>
#include <QSim/Execution/ServerPool.h>

using namespace QSim;
using namespace Eigen;

using namespace std::chrono_literals;


class Worker : public ServerPoolWorker
{
public:
    Worker() : ServerPoolWorker(6) {}

    virtual DataPackagePayload DoWork(DataPackagePayload data) override
    {
        Print(std::string("Starting (") + (char*)data.GetData() + ")");
        std::this_thread::sleep_for(5s);

        return DataPackagePayload();
    }

    void Print(const std::string& str)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        std::cout << str << std::endl;
    }
private:
    std::mutex m_mutex;
};

class Master : public ServerPool
{
public:
    virtual void OnTaskCompleted(UUIDv4 id)
    {
        std::cout << "Task completed: " << id.ToString() << std::endl;
    }
};

int WorkerMain()
{
    Worker worker;
    if(!worker.Run(8000))
        return 1;
    return 0;
}

int MasterMain()
{
    Master pool;
    std::thread thread([&](){ pool.Run(); });
    
    pool.ConnectWorkerHostname("localhost", 8000);
    pool.ConnectWorkerHostname("calca", 8000);
    pool.ConnectWorkerHostname("calcb", 8000);
    pool.ConnectWorkerHostname("calcc", 8000);
    // pool.ConnectWorkerHostname("panama", 8000);
    // pool.ConnectWorker("192.168.2.2", 8000);

    if (pool.GetWorkerCount() == 0)
    {
        std::cout << "No worker found!" << std::endl;
        pool.Stop();
        thread.join();
        return 1;
    }

    for (int i=0; i<50; i++)
    {
        std::string str = "Task " + std::to_string(i);
        DataPackagePayload data(str.size());
        std::copy_n(str.data(), str.size(), (char*)data.GetData());
        pool.Submit(data);
    }

    std::cout << "Waiting..." << std::endl;
    pool.WaitUntilFinished();
    std::cout << "Finished" << std::endl;

    pool.Stop();
    thread.join();

    return 0;
}


int main(int argc, const char** argv)
{
    ArgumentParser parser;
    parser.AddOption("worker", "Start as worker");
    auto args = parser.Parse(argc, argv);
    
    if (args.IsOptionPresent("worker"))
        return WorkerMain();
    else 
        return MasterMain();

    return 0;
}
