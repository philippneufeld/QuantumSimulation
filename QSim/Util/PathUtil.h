// Philipp Neufeld, 2021-2022

#ifndef QSim_PathUtil_H_
#define QSim_PathUtil_H_

#include "../Platform.h"

#include <string>

namespace QSim
{
    std::string GetHomeDirPath();
    std::string GetMicCellFolder();
    std::string GetDefaultAppDir(const std::string& appName, bool create=true);
    std::string GetDefaultAppDataDir(const std::string& appName, bool create=true);

    std::string GetHostname();
    std::string GetTimestampString();
    std::string GenerateFilename(const std::string& baseName);
    
    bool MoveFile(const std::string& from, const std::string& to);
}

#endif
