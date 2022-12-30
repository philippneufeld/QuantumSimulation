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

class MyServer : public PackageServer
{
public:
    virtual void OnClientConnected(std::size_t id, const std::string& ip) override
    {
        std::cout << "Client connected: " << ip << std::endl;
    }
    virtual void OnClientDisconnected(std::size_t id) override 
    {
        std::cout << "Client disconnected" << std::endl;     
    }

    virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) override
    { 
        char* msg = reinterpret_cast<char*>(data.GetData());
        std::cout << "Message received (" << id << "): " << msg << std::endl;

        constexpr std::uint32_t n = 100000;
        std::vector<double> vec;
        vec.reserve(n);
        for (std::size_t i = 0; i < n; i++)
        {
            vec.push_back(i);
        }

        NetworkDataPackage reply(n*sizeof(double));
        std::copy_n(vec.data(), n, reinterpret_cast<double*>(reply.GetData()));
        this->WriteTo(id, reply);
    }
    
};

class MyClient : public PackageClient
{
public:
    virtual void OnMessageReceived(std::size_t id, NetworkDataPackage data) override
    {
        std::cout << "Received data: " << data.GetSize() / sizeof(double) << std::endl;
    }
};


class Worker : public ServerPoolWorker
{
public:
    virtual DataPackagePayload DoWork(DataPackagePayload data) override
    {
        Print("Starting task....");
        std::this_thread::sleep_for(2s);
        Print((char*)data.GetData());
        std::this_thread::sleep_for(2s);
        Print("Finished task");

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

int ServerMain()
{
    MyServer server;
    if (!server.Run(8000))
        return 1;
    return 0;
}

int ClientMain()
{
    std::string helloServer = "Hello server!";

    MyClient client;

    std::thread thread([&](){ client.Run(); });

    //std::this_thread::sleep_for(50s);

    std::size_t id = client.ConnectHostname("panama", 8000);
    if (id == 0)
    {
        std::cout << "Failed to connect to server" << std::endl;
        client.Stop();
        thread.join();
        return 1;
    }
    
    NetworkDataPackage msg(helloServer.size());
    std::copy_n(helloServer.begin(), helloServer.size(), reinterpret_cast<char*>(msg.GetData()));
    client.WriteTo(id, std::move(msg));
    
    std::this_thread::sleep_for(5s);

    std::cout << "Stopping..." << std::endl;
    client.Stop();
    thread.join();

    return 0;
}

int WorkerMain()
{
    Worker worker;
    if(!worker.Run(8000))
        return 1;
    return 0;
}

int MasterMain()
{
    ServerPool pool;
    std::thread thread([&](){ pool.Run(); });
    
    pool.ConnectWorkerHostname("panama", 8000);
    pool.ConnectWorker("192.168.2.2", 8000);

    if (pool.GetWorkerCount() == 0)
    {
        std::cout << "No worker found!" << std::endl;
        pool.Stop();
        thread.join();
        return 1;
    }

    std::string str = "Hello world!";
    DataPackagePayload data(str.size());
    std::copy_n(str.data(), str.size(), (char*)data.GetData());

    for (int i=0; i<10; i++)
    {
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
