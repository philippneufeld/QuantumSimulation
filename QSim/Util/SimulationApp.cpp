// Philipp Neufeld, 2021-2022

#include "../Platform.h"

#include <iostream>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <filesystem>

#if defined(QSim_PLATFORM_LINUX)
#include <unistd.h>
#endif

#include "SimulationApp.h"
#include "Argparse.h"

namespace QSim
{

    std::string IncrementName(const std::string& name)
    {
        std::size_t idx = 2;
        auto found = name.find_last_not_of("0123456789") + 1;
        std::string cntstr = name.substr(found);
        std::string base = name.substr(0, found);
        if (!cntstr.empty())
            idx = std::stoull(cntstr) + 1;
        return base + std::to_string(idx);
    }


    SimulationApp::SimulationApp()
    {

    }

    SimulationApp::~SimulationApp()
    {

    }

    int SimulationApp::Run(int argc, const char** argv)
    {
        // get home directory
        std::string homeDir = ".";
        char* szHomeDir = std::getenv("HOME");
        if (szHomeDir) homeDir.assign(szHomeDir);
        if (!homeDir.empty() && (homeDir.back() == '/' || homeDir.back() == '\\'))
            homeDir.pop_back();


        // find remote_home directory and fallback to normal home if it could not be found
        std::string remoteHome = homeDir + "/remote_home";
        if (!std::filesystem::is_directory(remoteHome))
            remoteHome = homeDir;

        // generate defaults
        std::string defaultFileName = GenerateDefaultFilename(argc, argv);
        std::string defaultDir = remoteHome + "/SimulationData/" + ExtractProgramName(argc, argv);

        // parse command line arguents
        ArgumentParser parser;
        parser.AddOptionDefault("f,file", "Absolute or relative path (from data directory) to data file.", defaultFileName);
        parser.AddOptionDefault("d,dir", "Data directory path.", defaultDir);
        parser.AddOption("c,continue", "Name of the simulation to continue with.");
        parser.AddOption("p,plot", "Enable plotting.");
        parser.AddOption("h,help", "Print this help string.");
        auto cmdArgs = parser.Parse(argc, argv);

        if (cmdArgs.IsError())
        {
            std::cout << cmdArgs.GetError() << std::endl;
            return 1;
        }
        else if (cmdArgs.IsOptionPresent("help"))
        {
            std::cout << parser.GetHelpString() << std::endl;
            return 0;
        }

        // replace leading tilde by home directory path
        std::string dirPathStr = cmdArgs.GetOptionStringValue("dir");
        std::string filePathStr = cmdArgs.GetOptionStringValue("file");
        if ((dirPathStr.find("~") == 0 && dirPathStr.size() == 1) || dirPathStr.find("~/") == 0)
            dirPathStr.replace(dirPathStr.begin(), dirPathStr.begin() + 1, homeDir);
        if ((filePathStr.find("~") == 0 && dirPathStr.size() == 1) || dirPathStr.find("~/") == 0)
            filePathStr.replace(filePathStr.begin(), filePathStr.begin() + 1, homeDir);

        // create path for data file
        std::filesystem::path dirPath = dirPathStr;
        std::filesystem::path filePath = filePathStr;
        if (filePath.is_relative())
            filePath = (dirPath / filePath);
        filePath = filePath.lexically_normal();
        dirPath = filePath.parent_path();

        // make sure that the parent directory exists
        if (!std::filesystem::is_directory(dirPath))
        {
            if (!std::filesystem::create_directories(dirPath))
            {
                std::cout << "Failed to create directory " << dirPath << std::endl;
                return 1;
            }
        }

        // create/open datafile
        DataFileOpenFlag openFlag = DataFile_DEFAULT;
        if (cmdArgs.IsOptionPresent("continue"))
            openFlag = DataFile_MUST_EXIST;

        DataFile dataFile;
        if (!dataFile.Open(filePath.string(), openFlag))
        {
            std::cout << "Failed to create file " << filePath << std::endl;
            return 1;
        }

        std::cout << "Opened data file " << filePath << std::endl;

        DataFileGroup root = dataFile.OpenRootGroup();
        std::vector<std::string> archivedSims = root.EnumerateSubgroups();

        DataFileGroup simData;
        if (cmdArgs.IsOptionPresent("continue"))
        {
            std::string cont = cmdArgs.GetOptionStringValue("continue");
            if (!root.DoesSubgroupExist(cont))
            {
                std::cout << "Simulation data for \"" + cont + "\" not found" << std::endl;
                return 1;
            }

            simData = root.GetSubgroup(cont);
        }
        else
        {
            // find unique name
            std::size_t idx = 0;
            for (const auto& name: archivedSims)
            {
                if (name.find("Simulation_") != 0)
                    continue;

                auto found = name.find_last_not_of("0123456789");
                std::string cntstr = name.substr(found + 1);
                if (cntstr.empty())
                    continue;

                std::size_t cnt = std::stoul(cntstr);
                idx = cnt > idx ? cnt : idx;
            }

            std::string name = "Simulation_" + std::to_string(idx + 1);
            simData = root.CreateSubgroup(name);

            if (!simData.IsValid())
            {
                std::cout << "Could not create new data group \"" + name + "\"" << std::endl;
                return 1;
            }

            std::cout << "Created data group \"" + name + "\"" << std::endl;

            // Initialize new simulation
            this->Init(simData);
        }

        if (!this->IsFinished(simData))
            this->Continue(simData);

        if (true || cmdArgs.IsOptionPresent("plot"))
            this->Plot(simData);

        return 0;
    }

    bool SimulationApp::IsFinished(DataFileGroup& simdata)
    {
        return simdata.DoesAttributeExist("Finished");
    }

    void SimulationApp::SetFinished(DataFileGroup& simdata)
    {
        simdata.CreateAttribute("Finished", {});
    }

    std::string SimulationApp::ExtractProgramName(int argc, const char** argv)
    {
        // extract name of the program
        std::string progName;
        if (argc > 0)
        {
            progName.assign(argv[0]);
            
            std::size_t startPos = progName.find_last_of("/\\") + 1;
            progName.erase(progName.begin(), progName.begin() + startPos);
            
            std::size_t endPos = progName.find_last_of(".");
            if (endPos > 0 && endPos < progName.size()) 
                progName.erase(progName.begin() + endPos, progName.end());
        }
        if (progName.empty())
            progName = "Unknwon";

        return progName;
    }

    std::string SimulationApp::GenerateDefaultFilename(int argc, const char** argv)
    {
        std::string defaultFileName;

        defaultFileName += ExtractProgramName(argc, argv);

        // get hostname
#if defined(QSim_PLATFORM_LINUX)
        char hostname[1024] = "";
        gethostname(hostname, sizeof(hostname));
        defaultFileName.push_back('_');
        defaultFileName.append(hostname);
#endif

        // append date
        time_t now = std::time(nullptr);
        tm* ptime = std::localtime(&now);
        std::size_t datestamp = (1900 + ptime->tm_year) * 10000 + (ptime->tm_mon + 1) * 100 + ptime->tm_mday;
        defaultFileName += '_' + std::to_string(datestamp) + ".h5";

        return defaultFileName;
    }

}
