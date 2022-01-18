// Philipp Neufeld, 2021-2022

#include "../Platform.h"

#include <iostream>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <filesystem>

#if defined(QSim_PLATFORM_LINUX)
#include <unistd.h>
#endif

#include "CalcApp2.h"
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
        std::string defaultFileName = GenerateDefaultFilename(argc, argv);
        std::string defaultDir = QSim_DATA_DIR;
        defaultDir += "/" + ExtractProgramName(argc, argv);

        // parse command line arguents
        QSim::ArgumentParser parser;
        parser.AddOptionDefault("f,file", "Absolute or relative path (from data directory) to data file.", defaultFileName);
        parser.AddOptionDefault("d,dir", "Data directoty.", defaultDir);
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

        // create path for data file and make sure that the parent directory exists
        std::filesystem::path dirPath = cmdArgs.GetOptionStringValue("dir");
        std::filesystem::path filePath = cmdArgs.GetOptionStringValue("file");
        
        if (filePath.is_relative())
            filePath = (dirPath / filePath);
        filePath.lexically_normal();
        dirPath = filePath.parent_path();

        if (!std::filesystem::is_directory(dirPath))
        {
            if (!std::filesystem::create_directories(dirPath))
            {
                std::cout << "Failed to create directory " << dirPath << std::endl;
                return 1;
            }
        }

        // create/open datafile
        DataFile3OpenFlag openFlag = DataFile3_DEFAULT;
        if (cmdArgs.IsOptionPresent("continue"))
            openFlag = DataFile3_MUST_EXIST;

        DataFile3 dataFile;
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

            // Initialize new simulation
            this->Init(simData);
        }

        if (!this->IsFinished(simData))
            this->Continue(simData);

        if (cmdArgs.IsOptionPresent("plot"))
            this->Plot(simData);

        return 0;
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
