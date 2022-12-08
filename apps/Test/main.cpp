// Philipp Neufeld, 2021-2022

#include <iostream>
#include <chrono>
#include <list>
#include <Eigen/Dense>

#include <sys/select.h>

#include <QSim/Util/Argparse.h>
#include <QSim/Execution/ServerPool.h>

using namespace QSim;
using namespace Eigen;

using namespace std::chrono_literals;

class MyServer : public TCPIPServer
{
public:
    virtual bool OnClientConnected(std::size_t id, const std::string& ip) override
    {
        std::cout << "Client connected: " << ip << std::endl;
        return true; 
    }
    virtual void OnClientDisconnected(std::size_t id) override 
    {
        std::cout << "Client disconnected" << std::endl;     
    }

    virtual SocketDataPackage OnMessageReceived(std::size_t id, SocketDataPackage data) override
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

        SocketDataPackage reply(n*sizeof(double));
        std::copy_n(vec.data(), n, reinterpret_cast<double*>(reply.GetData()));

        return reply;
    }
    
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

    TCPIPClient client;
    if (!client.ConnectHostname("ludwigsburg", 8000))
    // if (!client.ConnectHostname("calcc", 8000))
    {
        std::cout << "Failed to connect to server" << std::endl;
        return 1;
    }
    
    auto reply = client.Query(helloServer.data(), helloServer.size());
    if (reply)
        std::cout << "Received data: " << reply.GetSize() / sizeof(double) << std::endl;

    std::this_thread::sleep_for(5s);

    return 0;
}

int main(int argc, const char** argv)
{
    fd_set set;
    FD_ZERO(&set);
    constexpr auto askjdfh = FD_SETSIZE;

    ArgumentParser parser;
    parser.AddOption("worker", "Start as worker");
    auto args = parser.Parse(argc, argv);
    
    if (args.IsOptionPresent("worker"))
        return ClientMain();
    else 
        return ServerMain();

    return 0;
}
