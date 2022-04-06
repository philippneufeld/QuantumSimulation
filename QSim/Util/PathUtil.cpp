// Philipp Neufeld, 2021-2022

#include "PathUtil.h"

#include <filesystem>
#include <chrono>
#include <ctime>

#if defined(QSim_PLATFORM_LINUX)
#include <unistd.h>
#endif

namespace QSim
{
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

    std::string GetHomeDirSubfolderPath(const std::string& name)
    {
        // find remote_home directory and fallback to normal home if it could not be found
        std::string homeDir = GetHomeDirPath();
        std::string remoteHome = homeDir + "/" + name;
        if (!std::filesystem::is_directory(remoteHome))
            remoteHome = homeDir;
        return remoteHome;
    }

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

    std::string GetTimestampString()
    {
        time_t now = std::time(nullptr);
        tm* ptime = std::localtime(&now);
        std::size_t datestamp = (1900 + ptime->tm_year) * 10000 + (ptime->tm_mon + 1) * 100 + ptime->tm_mday;
        std::size_t timestamp = (ptime->tm_hour * 100 + ptime->tm_min) * 100 + ptime->tm_sec;
        return std::to_string(datestamp) + '-' + std::to_string(timestamp);
    }

    std::string GenerateFilename(const std::string& baseName)
    {
        return baseName + (baseName.empty() ? "" : "_") + GetTimestampString() + "_" + GetHostname();
    }

}
