// Philipp Neufeld, 2021-2022

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include <QSim/Util/Argparse.h>
#include <QSim/Execution/ServerPool.h>

using namespace QSim;
using namespace Eigen;

using namespace std::chrono_literals;

class TestDataPackage : public SocketDataPackage
{
public:
    TestDataPackage() : m_iPayload(0), m_dPayload(0.0) {}
    TestDataPackage(int i, double d) : m_iPayload(i), m_dPayload(d) {}


    virtual SocketDataPackageBin Serialize() override
    {
        SocketDataPackageBin bin;
        bin.Allocate(sizeof(int) + sizeof(double));
        *reinterpret_cast<int*>(bin.GetData()) = m_iPayload;
        *reinterpret_cast<double*>(bin.GetData() + sizeof(int)) = m_dPayload;
        return bin;
    }

    virtual void Deserialize(const SocketDataPackageBin& data) override
    {
        m_iPayload = *reinterpret_cast<int*>(data.GetData());
        m_dPayload = *reinterpret_cast<double*>(data.GetData() + sizeof(int));
    }

private:
    int m_iPayload;
    double m_dPayload;
};

int ServerMain()
{
    std::string helloClient = "Hello client!";

    TCPIPSocketServer server;

    if (!server.Bind(8000))
        std::cout << "Failed to bind server " << strerror(errno) << std::endl;

    if (!server.Listen(5))
        std::cout << "Failed to listen on port" << std::endl;

    while (true)
    {
        auto [client, clientIp] = server.Accept();

        std::cout << "Connection accepted: " << clientIp << " " << client.IsValid() << std::endl;
        client.Send(helloClient.data(), helloClient.size());
    }
    std::cout << "server" << std::endl;
    return 0;
}

int ClientMain()
{
    std::string helloServer = "Hello server!";

    TCPIPSocketClient client;
    if (!client.ConnectHostname("calcc", 8000))
    {
        std::cout << "Failed to connect to server" << std::endl;
        return 1;
    }
    client.Send(helloServer.data(), helloServer.size());

    char buff[4096] = "";
    client.Recv(buff, 4096);
    std::cout << buff << std::endl;

    return 0;
}

int main(int argc, const char** argv)
{

    ArgumentParser parser;
    parser.AddOption("worker", "Start as worker");
    auto args = parser.Parse(argc, argv);
    
    if (args.IsOptionPresent("worker"))
        return ClientMain();
    else 
        return ServerMain();

    /*ThreadPool pool(5);

    TestDataPackage test1(1, 4.5);
    auto testBin = test1.Serialize();
    
    TestDataPackage test2;
    test2.Deserialize(testBin);

    pool.Submit([](){ std::this_thread::sleep_for(2s); });
    for (int i=0; i<4; i++) pool.Submit([](){ std::this_thread::sleep_for(5s); });

    pool.WaitUntilReadyForTask();
    std::cout << "Ready for task" << std::endl;

    pool.WaitUntilFinished();
    std::cout << "Finished!" << std::endl;*/

    return 0;
}
