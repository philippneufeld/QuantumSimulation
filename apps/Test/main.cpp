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
    virtual bool OnClientConnected(const std::string& ip) override
    {
        std::cout << "Client connected: " << ip << std::endl;
        return true; 
    }
    virtual void OnClientDisconnected() override 
    {
        std::cout << "Client disconnected" << std::endl;     
    }

    virtual SocketDataPackage GetWelcomeMessage() override
    { 
        return SocketDataPackage(); 
    }

    virtual SocketDataPackage OnMessageReceived(SocketDataPackage data) override
    { 
        char* msg = reinterpret_cast<char*>(data.GetData());
        std::cout << "Message received: " << msg << std::endl;
        return SocketDataPackage(); 
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

    TCPIPClientSocket client;
    if (!client.ConnectHostname("ludwigsburg", 8000))
    {
        std::cout << "Failed to connect to server" << std::endl;
        return 1;
    }
    client.Send(helloServer.data(), helloServer.size());


    std::this_thread::sleep_for(5s);

    // char buff[4096] = "";
    // client.Recv(buff, 4096);
    // std::cout << buff << std::endl;

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
