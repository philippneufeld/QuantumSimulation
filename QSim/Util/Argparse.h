// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_Argparse_H_
#define QSim_Util_Argparse_H_

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>

#include "../Platform.h"

namespace QSim
{

    class ArgumentParserResult
    {
    public:
        ArgumentParserResult(const std::map<std::string, std::string>& shortNameMap) 
            : m_isError(false), m_shortNameMap(shortNameMap) { }
        ArgumentParserResult(const std::string& error) 
            : m_isError(true), m_errorString(error) { }

        // setter
        void SetOption(const std::string& optName, const std::string& optVal);
        
        // getter
        bool IsOptionPresent(std::string optName) const;
        std::string GetOptionStringValue(std::string optName) const;
        template<typename T>
        T GetOptionValue(std::string optName) const;

        // error
        bool IsError() const { return m_isError; }
        const std::string& GetError() const { return m_errorString; }

    private:
        bool m_isError;
        std::string m_errorString;
        std::map<std::string, std::string> m_options;
        std::map<std::string, std::string> m_shortNameMap;
    };

    namespace Internal
    {
        struct ArgumentParserOption
        {
            std::string description;
            bool useDefaultValue;
            std::string defaultValue;
        };
    }

    class ArgumentParser
    {
    public:
        ArgumentParser() = default;

        bool AddOption(const std::string& optName, const std::string& desc);
        bool AddOptionDefault(const std::string& optName, const std::string& desc, const std::string& defaultVal);
    private:
        bool AddOptionHelper(const std::string& optName, const std::string& desc, 
            bool useDefault, const std::string& defaultVal);
        
    public:
        ArgumentParserResult Parse(int argc, const char** argv);

        std::string GetHelpString() const;

    private:
        std::map<std::string, std::string> m_shortNameMap;
        std::map<std::string, Internal::ArgumentParserOption> m_options;
    };


    namespace Internal
    {
        template<typename T>
        struct ArgumentValueParser
        {
            static T Parse(const std::string& str)
            {
                T res{};
                if (!str.empty())
                {
                    std::istringstream ss(str);
                    ss >> res;
                }
                return res;
            }
        };

        template<typename T>
        struct ArgumentValueParser<std::vector<T>>
        {
            static std::vector<T> Parse(const std::string& str)
            {
                std::vector<T> vect{};
                T val{};
                if (!str.empty())
                {
                    std::istringstream ss(str);
                    if (ss.peek() == '[')
                        ss.ignore();

                    for (T val{}; ss >> val;) {
                        vect.push_back(val);    
                        if (ss.peek() == ',')
                            ss.ignore();
                    }
                }
                return vect;
            }
        };
    }


    // template function implementationd
    template<typename T>
    T ArgumentParserResult::GetOptionValue(std::string optName) const
    {
        auto str = GetOptionStringValue(optName);
        return Internal::ArgumentValueParser<T>::Parse(str);
    }

}

#endif
