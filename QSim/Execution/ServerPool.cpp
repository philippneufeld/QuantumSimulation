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
    //
    // SocketDataPackageBin
    //

    SocketDataPackageBin::SocketDataPackageBin() 
        : m_size(0), m_pData(nullptr) {}
    
    SocketDataPackageBin::~SocketDataPackageBin() 
    { 
        Allocate(0);
    }

    SocketDataPackageBin::SocketDataPackageBin(const SocketDataPackageBin& rhs)
        : SocketDataPackageBin()
    {
        if (rhs.m_size > 0)
        {
            Allocate(rhs.m_size);
            std::copy(rhs.m_pData, rhs.m_pData+rhs.m_size, m_pData);
        }
    }
    
    SocketDataPackageBin::SocketDataPackageBin(SocketDataPackageBin&& rhs)
        : m_pData(rhs.m_pData), m_size(rhs.m_size)
    {
        rhs.m_pData = nullptr;
        rhs.m_size = 0;
    }

    SocketDataPackageBin& SocketDataPackageBin::operator=(const SocketDataPackageBin& rhs)
    {
        SocketDataPackageBin tmp(rhs);
        std::swap(*this, tmp);
        return *this;
    }

    SocketDataPackageBin& SocketDataPackageBin::operator=(SocketDataPackageBin&& rhs)
    {
        std::swap(m_pData, rhs.m_pData);
        std::swap(m_size, rhs.m_size);
        return *this;
    }
    
    bool SocketDataPackageBin::Allocate(std::uint64_t size)
    {
        if (size == m_size)
            return true;

        if (m_pData)
            delete[] m_pData;
        m_pData = nullptr;
        m_size = 0;      
        
        if (size > 0)
        {
            m_pData = new std::uint8_t[size];
            m_size = size;
            if (!m_pData)
            {
                m_pData = nullptr;
                m_size = 0;
            }
        }

        return (size == m_size);
    }


    //
    // TCPIPSocket
    //
    
    TCPIPSocket::TCPIPSocket(int handle)
        : m_handle(handle) {}

    TCPIPSocket::TCPIPSocket()
        : TCPIPSocket(socket(AF_INET, SOCK_STREAM, 0)) {}

    TCPIPSocket::~TCPIPSocket()
    {
        Close();
    }
    
    TCPIPSocket::TCPIPSocket(TCPIPSocket&& rhs)
        : TCPIPSocket()
    {
        std::swap(m_handle, rhs.m_handle);
    }

    TCPIPSocket& TCPIPSocket::operator=(TCPIPSocket&& rhs)
    {
        std::swap(m_handle, rhs.m_handle);
        return *this;
    }

    void TCPIPSocket::Close()
    {
        if (IsValid())
            close(m_handle);
        m_handle = -1;
    }


    //
    // TCPIPSocketServer
    //

    TCPIPSocketServer::TCPIPSocketServer(int handle)
        : TCPIPSocket(handle) {}

    TCPIPSocketServer::TCPIPSocketServer()
        : TCPIPSocket() {}

    bool TCPIPSocketServer::Bind(short port)
    {
        if (!IsValid())
            return false;

        // lets server reuse previously bound address
        const int enable = 1;
        if (setsockopt(GetHandle(), SOL_SOCKET, SO_REUSEADDR, &enable, sizeof(enable)) < 0)
            return false;

        sockaddr_in addr = { 0 };
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY; // bind to all interfaces
        addr.sin_port = htons(port);

        if (bind(GetHandle(), reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
            return false;

        return true;
    }

    bool TCPIPSocketServer::Listen(int queue_length)
    {
        if (!IsValid())
            return false;

        if (listen(GetHandle(), queue_length) < 0)
            return false;

        return true;
    }

    std::tuple<TCPIPConnection, std::string> TCPIPSocketServer::Accept()
    {
        sockaddr_in addr = { 0 };
        addr.sin_family = AF_INET;
        socklen_t len = sizeof(addr);

        int clientHandle = accept(GetHandle(), reinterpret_cast<sockaddr*>(&addr), &len);
        std::string ip = inet_ntoa(addr.sin_addr);

        return std::make_tuple(TCPIPConnection(clientHandle), ip);
    }


    //
    // TCPIPConnection
    //

    TCPIPConnection::TCPIPConnection(int handle)
        : TCPIPSocket(handle) {}

    TCPIPConnection::TCPIPConnection() 
        : TCPIPSocket() {}

    std::size_t TCPIPConnection::Send(const void* data, std::size_t n)
    {
        if (!IsValid())
            return 0;

        return send(GetHandle(), data, n, 0);
    }

    //
    // TCPIPSocketClient
    //

    TCPIPSocketClient::TCPIPSocketClient()
        : TCPIPConnection() {}

    bool TCPIPSocketClient::Connect(const std::string& ip, short port)
    {
        sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);

        std::istringstream iss(ip);
        for (int i=0; i<4; i++)
        {
            std::string tmp;
            std::getline(iss, tmp, '.');
            reinterpret_cast<std::uint8_t*>(&addr.sin_addr.s_addr)[i] = std::atoi(tmp.c_str());
        }

        if (connect(GetHandle(), reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
            return false;

        return true;
    }

    bool TCPIPSocketClient::ConnectHostname(const std::string& hostname, short port)
    {
        constexpr int bufferLen = 256;
        char buffer[bufferLen] = "";
        gethostname(buffer, bufferLen - 1);

        hostent* pRes = gethostbyname(hostname.c_str());
        char** addresses = pRes->h_addr_list;
        if (pRes->h_length < 4)
            return false;

        std::string ip = 
            std::to_string((int)(unsigned char)addresses[0][0]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][1]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][2]) + "." 
            + std::to_string((int)(unsigned char)addresses[0][3]);
        
        return Connect(ip, port);
    }


}
