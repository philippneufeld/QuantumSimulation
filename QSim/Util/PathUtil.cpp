// Philipp Neufeld, 2021-2022

#include "PathUtil.h"

#include <filesystem>
#include <algorithm>
#include <chrono>
#include <ctime>

#if defined(QSim_PLATFORM_LINUX)
#include <unistd.h>
#endif

namespace QSim
{
    
    std::string GetHostname()
    {
#if defined(QSim_PLATFORM_LINUX)
        char hostname[1024] = "";
        gethostname(hostname, sizeof(hostname));
        return hostname;
#else
        return "UNKNOWN"
#endif
    }

    std::string GetHomeDirPath()
    {
        // get home directory by reading the HOME environment variable
        std::string homeDir = ".";
        char* szHomeDir = std::getenv("HOME");
        if (szHomeDir) homeDir.assign(szHomeDir);
        if (!homeDir.empty() && (homeDir.back() == '/' || homeDir.back() == '\\'))
            homeDir.pop_back();
        return homeDir;
    }

    std::string GetMicCellFolder()
    {
        std::string homeDir = GetHomeDirPath();
        std::string host = GetHostname();
        
        // to lowercase
        std::transform(host.begin(), host.end(), host.begin(),
            [](unsigned char c){ return std::tolower(c); });

        if (host == "ludwigsburg" || host == "sost")
            return homeDir + "/MicCells";
        else if (host.size() == 5 && host.substr(0, 4) == "calc")
            return "/mnt/groups/MicCells";
        else
            return homeDir + "/remote_groups/MicCells";
    }

    std::string GetDefaultAppDir(const std::string& appName, bool create)
    {
        std::string path = GetMicCellFolder();
        path.append("/TraceGasSensing/Simulation/QNOSE_QuantumSimulation/");
        if (!appName.empty())
            path.append(appName);
        else
            path.append("Unknown");

        if (create)
            std::filesystem::create_directories(path);

        return path;
    }

    std::string GetDefaultAppDataDir(const std::string& appName, bool create)
    {
        std::string path = GetDefaultAppDir(appName, create) + "/Data";
        if (create)
            std::filesystem::create_directories(path);
        return path;
    }

    std::string GetTimestampString()
    {
        time_t now = std::time(nullptr);
        tm* ptime = std::localtime(&now);
        std::size_t datestamp = (1900 + ptime->tm_year) * 10000 + (ptime->tm_mon + 1) * 100 + ptime->tm_mday;
        std::size_t timestamp = (ptime->tm_hour * 100 + ptime->tm_min) * 100 + ptime->tm_sec;
        return std::to_string(datestamp) + (ptime->tm_hour < 10 ? "-0" : "-") + std::to_string(timestamp);
    }

    std::string GenerateFilename(const std::string& baseName)
    {
        return baseName + (baseName.empty() ? "" : "_") + GetTimestampString() + "_" + GetHostname();
    }

    bool MoveFile(const std::string& from, const std::string& to)
    {
        try 
        {
            std::filesystem::rename(from, to);
        } 
        catch (std::filesystem::filesystem_error&)
        {
            try 
            {
                std::filesystem::copy(from, to);
                std::filesystem::remove(from);
            }
            catch (std::filesystem::filesystem_error&)
            {
                return false;
            }
        }

        return true;
    }

}
