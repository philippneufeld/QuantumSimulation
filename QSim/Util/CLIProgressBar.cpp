// Philipp Neufeld, 2021-2022

#include "CLIProgressBar.h"

#ifdef QSim_PLATFORM_LINUX
#include <sys/ioctl.h> //ioctl() and TIOCGWINSZ
#include <unistd.h> // for STDOUT_FILENO
#endif

#include <algorithm>

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

    CLIProgBar::CLIProgBar()
        : m_stopThread(false), m_started(0), m_progress(0), 
        m_width(GetCLIWidth()), m_progressChar('#') { }

    CLIProgBar::~CLIProgBar()
    {
        Stop();
    }

    void CLIProgBar::Start() 
    {
        Stop();
        m_stopThread = false;
        m_thread = std::thread([this]() { this->ThreadFunc(); });
        m_startTs = std::chrono::high_resolution_clock::now();
        m_prevTs = m_startTs;
    }

    void CLIProgBar::WaitUntilFinished()
    {
        if (m_thread.joinable())
            m_thread.join();
        m_thread = std::thread();
    }

    void CLIProgBar::Stop()
    {
        if (m_progress < 1.0)
            m_stopThread = true;
        Update();
        WaitUntilFinished();
    }

    void CLIProgBar::Update()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_wakeUp.notify_all();
    }

    void CLIProgBar::SetWidth(std::size_t width)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_width = width;
        m_wakeUp.notify_all();
    }

    void CLIProgBar::SetProgressCharacter(char progChar)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_progressChar = progChar;
        m_wakeUp.notify_all();
    }

    void CLIProgBar::SetTitle(const std::string& title)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_title = !title.empty() ? title + ": " : title;
        m_wakeUp.notify_all();
    }

    double CLIProgBar::GetProgress() const
    {
        return m_progress;
    }

    void CLIProgBar::SetProgress(double prog) 
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_progress = prog;
        m_prevTs = std::chrono::high_resolution_clock::now();
        m_wakeUp.notify_all();
    }

    std::string CLIProgBar::GetDescription() const
    {
        return m_title;
    }

    std::string CLIProgBar::GetProgressText() const
    {
        std::string strPerc = std::to_string(static_cast<std::size_t>(m_progress * 100));
        
        std::string str;
        str.reserve(5);

        str.append(4 - strPerc.size(), ' ');
        str.append(strPerc);
        str.push_back('%');
        return str;
    }

    std::string CLIProgBar::GetTimeText() const
    {
        auto dt = std::chrono::high_resolution_clock::now() - m_startTs;
        double secs = dt.count() / 1e9;
        std::string curr = SecondsToString(static_cast<std::size_t>(secs));

        // 14 character string
        std::string str;
        str.reserve(14);


        str.append(" [");
        if (m_progress < 1.0)
        {
            double projSecs = m_progress > 0 ? (m_prevTs - m_startTs).count() / m_progress / 1e9 : 0.0;
            std::string tot = SecondsToString(static_cast<std::size_t>(projSecs));

            str.append(curr);
            str.push_back('<');
            str.append(tot);
            
        }
        else
        {  
            curr.reserve(11);
            curr += ':';
            
            int ms = static_cast<int>((secs - static_cast<int>(secs)) * 1000);
            ms = (ms > 999 ? 999 : (ms < 0 ? 0 : ms));
            std::string strms = std::to_string(ms);
            strms.insert(0, 3 - strms.size(), '0');
            curr += strms;

            str.append((11 - curr.size()) / 2, ' ');
            str += curr;
            str.append(13-str.size(), ' ');
        }
        str.push_back(']');

        return str;
    }

    std::string CLIProgBar::SecondsToString(std::size_t secs)
    {
        // 5 character string "mm:ss"
        std::string str;
        str.reserve(10);

        std::string strMins = std::to_string(secs / 60);
        std::string strSecs = std::to_string(secs % 60);
        if (strMins.size() < 2)
            str.append(2 - strMins.size(), '0');
        str.append(strMins);
        str.push_back(':');
        if (strSecs.size() < 2)
            str.append(2 - strSecs.size(), '0');
        str.append(strSecs);

        return str;
    }

    void CLIProgBar::ThreadFunc()
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        
        bool firstIteratoin = true;
        while (!m_stopThread)
        {
            if (!firstIteratoin)
            {
                m_wakeUp.wait_for(lock, 0.25s);
                if (m_stopThread)
                    break;
            }
            firstIteratoin = false;

            PrintBar();

            if (m_progress == 1.0)
                m_stopThread = true;
        }
        std::cout << std::endl;
    }

    void CLIProgBar::PrintBar()
    {
        std::string str;
        str.reserve(m_width);
        
        auto front = GetDescription();
        auto back = GetProgressText() + GetTimeText();

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
        std::size_t iProg = static_cast<std::size_t>(barWidth * m_progress);
        str.append(iProg, m_progressChar);
        str.append(barWidth - iProg, ' ');
        str.append("]");
        str.append(back);

        std::cout << "\u001b[1000D" << str << std::flush;
    }



    CLIProgBarInt::CLIProgBarInt()
        : CLIProgBarInt(100) { }
        
    CLIProgBarInt::CLIProgBarInt(std::size_t total)
        : m_cnt(0), m_total(total) { }

    void CLIProgBarInt::SetTotal(std::size_t total)
    {
        m_total = total;
        this->SetProgress(static_cast<double>(m_cnt) / m_total);
    }
    
    void CLIProgBarInt::SetCount(std::size_t cnt)
    {
        m_cnt = std::min<std::size_t>(cnt, m_total);
        this->SetProgress(static_cast<double>(m_cnt) / m_total);
    }

    void CLIProgBarInt::IncrementCount()
    {
        std::size_t cnt = ++m_cnt; // atomically increment
        std::size_t total = m_total;
        if (cnt < total)
            this->SetProgress(static_cast<double>(cnt) / m_total);
        else
            this->SetProgress(1.0);
    }

    std::string CLIProgBarInt::GetProgressText() const
    {
        auto cnt = std::to_string(m_cnt);
        auto tot = std::to_string(m_total);
        std::size_t len = 2*tot.size() + 2;

        std::string str;
        str.reserve(len);

        str.append(1 + tot.size() - cnt.size(), ' ');
        str.append(cnt);
        str.push_back('/');
        str.append(tot);

        return str;
    }

}
