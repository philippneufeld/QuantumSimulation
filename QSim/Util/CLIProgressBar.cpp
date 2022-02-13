// Philipp Neufeld, 2021-2022

#include "../Platform.h"
#include "CLIProgressBar.h"

#ifdef QSim_PLATFORM_LINUX
#include <sys/ioctl.h> //ioctl() and TIOCGWINSZ
#include <unistd.h> // for STDOUT_FILENO
#endif

#include <iostream>

using namespace std::literals::chrono_literals;

namespace QSim
{

#ifdef QSim_PLATFORM_LINUX
    std::size_t GetCLIWidth()
    {
        struct winsize w;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
        return w.ws_col;
    }
#else
    std::size_t GetCLIWidth()
    {
        return 80;
    }
#endif

    CLIProgBar::CLIProgBar(std::size_t total, const std::string& title)
        : CLIProgBar(std::make_shared<Progress>(total), title) { }
    
    CLIProgBar::CLIProgBar(std::shared_ptr<Progress> pHandle, const std::string& title)
        : m_thread([&]() { WorkerThread(); }), m_stopThread(false), m_waitOnDtor(true),
        m_pHandle(pHandle ? pHandle : std::make_shared<Progress>(0)),
        m_title(title), m_width(GetCLIWidth()), m_progressChar('#') { }

    CLIProgBar::~CLIProgBar()
    {
        if (m_waitOnDtor)
            m_pHandle->WaitUntilFinished();
        m_stopThread = true;
        m_pHandle->SendSignal();
        WaitUntilFinished();
    }

    void CLIProgBar::WaitUntilFinished()
    {
        if (m_thread.joinable())
            m_thread.join();
    }

    void CLIProgBar::Update()
    {
        m_pHandle->SendSignal();
    }

    void CLIProgBar::IncrementCount(std::size_t inc)
    {
        m_pHandle->IncrementCount(inc);
    }

    std::string CLIProgBar::TimeToString(double secs, bool millis)
    {
        secs = secs > 0 ? secs : 0;

        std::string mm = std::to_string(static_cast<int>(secs)/60);
        std::string ss = std::to_string(static_cast<int>(secs)%60);

        if (mm.length() == 1)
            mm.insert(mm.begin(), '0');
        if (ss.length() == 1)
            ss.insert(ss.begin(), '0');

        std::string str;
        str.reserve(16);
        str.append(mm);
        str.push_back(':');
        str.append(ss);

        if (millis)
        {
            std::string ms = std::to_string(static_cast<int>((secs - static_cast<int>(secs)) * 1000));
            ms.insert(ms.begin(), std::min<std::size_t>(3 - ms.size(), 3), '0');
            str.push_back(':');
            str.append(ms);
        }

        return str;
    }

    std::string CLIProgBar::GetProgressText(std::size_t cnt, std::size_t tot)
    {
        std::string str;
        str.reserve(8);
        str.push_back(' ');

        std::size_t prog = (100 * cnt) / tot;
        str.append(std::to_string(prog));
        str.insert(str.begin() + 1, std::min<std::size_t>(4 - str.size(), 4), ' ');
        str.push_back('%');
        
        return str;
    }
    std::string CLIProgBar::GetProgressTextFrac(std::size_t cnt, std::size_t tot)
    {
        std::string str;
        str.reserve(16);
        str.push_back(' ');

        std::string totStr = std::to_string(tot);
        str.append(std::to_string(cnt));
        str.insert(str.begin() + 1, std::min<std::size_t>(totStr.size() + 1 - str.size(), totStr.size() + 1), ' ');
        str.push_back('/');
        str.append(totStr);
        
        return str;
    }

    std::string CLIProgBar::GetTimeText(double elapsedTime, double totalTimeEst, bool finished)
    {  
        std::string str;
        str.reserve(32);
        str.append(" [");

        if (!finished)
        {
            str.append(TimeToString(elapsedTime, false));
            str.push_back('<');
            str.append(TimeToString(totalTimeEst, false));
        }
        else
        {
            str.append(TimeToString(elapsedTime, true));
        }
        str.push_back(']');

        return str;
    }

    std::string CLIProgBar::GenerateBar(std::size_t cnt, std::size_t tot, 
            double elapsedTime, double totalTimeEst) const
    {
        std::string str;
        str.reserve(m_width);
        
        bool finished = !(cnt < tot);
        std::string front = m_title;
        std::string back = GetProgressText(cnt, tot) + GetTimeText(elapsedTime, totalTimeEst, finished);

        // cutoff front and back string if necessary
        std::size_t minWidth = 5;
        if (front.size() + back.size() > m_width - minWidth)
        {
            std::size_t exess = std::min(back.size(), front.size() + back.size() - m_width + minWidth);
            back.erase(back.end() - exess, back.end());

            exess = std::min(front.size(), front.size() + back.size() - m_width + minWidth);
            front.erase(front.end() - exess, front.end());
        }

        str.append(front);
        str.append("[");
        std::size_t barWidth = m_width - std::min(m_width - 2, front.size() + back.size()) - 2;
        std::size_t iProg = barWidth * cnt / tot;
        str.append(iProg, m_progressChar);
        str.append(barWidth - iProg, ' ');
        str.append("]");
        str.append(back);

        return str;
    }

    void CLIProgBar::WorkerThread()
    {
        std::string currBar;
        while (!m_stopThread)
        {
            auto [cnt, tot, elT, estTot] = m_pHandle->WaitForSignal(0.1);
            std::string bar = GenerateBar(cnt, tot, elT, estTot);

            if (bar != currBar)
            {
                std::cout << "\u001b[1000D" << bar << std::flush;
                currBar = bar;
            }

            if (cnt >= tot)
                m_stopThread = true;
        }
        std::cout << std::endl;
    }

}
