// Philipp Neufeld, 2021

#include "Argparse.h"

namespace QSim
{
    //
    // ArgumentParserResult
    //

    void ArgumentParserResult::SetOption(const std::string& optName, const std::string& optVal)
    {
        m_options[optName] = optVal;
    }

    bool ArgumentParserResult::IsOptionPresent(std::string optName) const
    {
        if (m_shortNameMap.find(optName) != m_shortNameMap.end())
            optName = m_shortNameMap.at(optName);
        return m_options.find(optName) != m_options.end();
    }

    std::string ArgumentParserResult::GetOptionStringValue(std::string optName) const
    {
        if (m_shortNameMap.find(optName) != m_shortNameMap.end())
            optName = m_shortNameMap.at(optName);
        return m_options.at(optName);
    }

    //
    // ArgumentParser
    //

    bool ArgumentParser::AddOption(const std::string& optName, const std::string& desc)
    {
        return AddOptionHelper(optName, desc, false, "");
    }

    bool ArgumentParser::AddOptionDefault(const std::string& optName, const std::string& desc, const std::string& defaultVal)
    {
        return AddOptionHelper(optName, desc, true, defaultVal);
    }

    bool ArgumentParser::AddOptionHelper(const std::string& optName, const std::string& desc, 
        bool useDefault, const std::string& defaultVal)
    {
        // split optName in long and short name
        std::string name, shortName;
        auto commaIt = std::find(optName.begin(), optName.end(), ',');
        if (commaIt == optName.end())
        {
            name = optName;
            if(name.length() == 1)
                shortName = name;
        }
        else
        {
            // too many commas -> wrong format
            if (std::find(commaIt + 1, optName.end(), ',') != optName.end())
                return false;

            shortName.assign(optName.begin(), commaIt);
            name.assign(commaIt + 1, optName.end());

            if (shortName.length() > 1)
                std::swap(name, shortName);

            if (shortName.length() > 1 || name.empty())
                return false;
        }

        // check if option with the same name is already present
        if (m_options.find(name) != m_options.end())
            return false;
        if (m_options.find(shortName) != m_options.end())
            return false;

        // check if option with the same short name is already present
        if (m_shortNameMap.find(name) != m_shortNameMap.end())
            return false;
        if (m_shortNameMap.find(shortName) != m_shortNameMap.end())
            return false;

        if (shortName.length() == 1)
            m_shortNameMap.insert({shortName, name});

        Internal::ArgumentParserOption option;
        option.description = desc;
        option.useDefaultValue = useDefault;
        option.defaultValue = defaultVal;
        m_options.insert({name, option});

        return true;
    }

    ArgumentParserResult ArgumentParser::Parse(int argc, const char** argv)
    {
        ArgumentParserResult result(m_shortNameMap);

        // set default values
        for (auto& opt: m_options)
        {
            if (opt.second.useDefaultValue)
                result.SetOption(opt.first, opt.second.defaultValue);
        }

        std::string currOption = "";
        for (int i = 1; i < argc; i++)
        {
            std::string tmp = argv[i];

            // Trim
            auto isWhitespace = [](char c) { return c == ' ' || c == '\t'; };
            auto it1 = std::find_if_not(tmp.begin(), tmp.end(), isWhitespace);
            auto it2 = std::find_if_not(tmp.rbegin(), tmp.rend(), isWhitespace).base();
            tmp.assign(it1, it2);

            if (tmp.empty())
                continue;

            if (tmp.find("-") == 0)
            {
                // get option name
                if (tmp.find("--") == 0)
                    currOption = tmp.substr(2);
                else
                {
                    auto it = m_shortNameMap.find(tmp.substr(1));
                    if (it == m_shortNameMap.end())
                        return ArgumentParserResult{"Invalid option name: " + tmp};
                    currOption = it->second;
                }
                
                if (m_options.find(currOption) == m_options.end())
                    return ArgumentParserResult{"Invalid option name: " + tmp};

                result.SetOption(currOption, "");
            }
            else
            {
                if (currOption.empty())
                    return ArgumentParserResult{"Option value must be preceeded by an option name"};

                result.SetOption(currOption, tmp);
                currOption.clear();
            }
        }

        return result;
    }

    std::string ArgumentParser::GetHelpString() const
    {
        std::string help;
        std::size_t descriptionCol = 25;
        std::size_t colMax = 80; 

        for (auto& opt: m_options)
        {
            std::string line = "   ";
            
            std::string shortName;
            for (auto& sn: m_shortNameMap)
            {
                if (sn.second == opt.first)
                {
                    shortName = sn.first;
                    break;
                }
            }
            
            if (!shortName.empty())
                line += "-" + shortName + ", ";
            line += "--" + opt.first;

            std::string spaces = " ";
            if (descriptionCol - 1 > line.length())
                spaces.resize(descriptionCol - line.length(), ' ');
            line += spaces;

            std::string descr = opt.second.description;
            if (opt.second.useDefaultValue)
                descr += " (Default value: \"" + opt.second.defaultValue + "\")";
            line.reserve(line.length() + colMax * (1 + descr.length() / (colMax - descriptionCol)) + 16);

            spaces.resize(descriptionCol, ' ');
            std::size_t handledLen = 0;
            for (std::size_t i = 0; i < descr.length(); i++)
            {
                if (line.length() - handledLen >= colMax)
                {
                    line.push_back('\n');
                    handledLen = line.length();
                    line += spaces;
                }
                line.push_back(descr[i]);
            }

            help += line + "\n";
        }
        
        if (!help.empty())
        {
            if (help.back() == '\n')
                help.pop_back();
        }

        return help;
    }

}
