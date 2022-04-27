// Philipp Neufeld, 2021-2022

#ifndef QSim_PathUtil_H_
#define QSim_PathUtil_H_

#include "../Platform.h"

#include <string>

namespace QSim
{
    std::string GetHomeDirPath();
    std::string GetHomeDirSubfolderPath(const std::string& name);

    std::string GetHostname();
    std::string GetTimestampString();

    std::string GenerateFilename(const std::string& baseName);

    bool MoveFile(const std::string& from, const std::string& to);
}

#endif
