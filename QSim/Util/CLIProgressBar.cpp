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

    CLIProgBar::CLIProgBar(std::size_t total, const std::string& title)
        : m_stopThread(false), m_started(false), m_total(total), m_title(title),
        m_width(GetCLIWidth()), m_progressChar('#') { }

    CLIProgBar::~CLIProgBar()
    {
        Stop();
    }

    void CLIProgBar::Start() 
    {
        Stop();
        m_stopThread = false;
        m_cnt = 0;
        m_startTs = std::chrono::high_resolution_clock::now();
        m_prevTs = m_startTs;
        m_thread = std::thread([this]() { this->WorkerThread(); });
    }

    void CLIProgBar::WaitUntilFinished()
    {
        if (m_thread.joinable())
            m_thread.join();
        m_thread = std::thread();
    }

    void CLIProgBar::Stop()
    {
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_stopThread = true;
        }
        
        m_wakeUp.notify_all();
        WaitUntilFinished();
    }

    void CLIProgBar::Update()
    {
        m_wakeUp.notify_all();
    }

    void CLIProgBar::SetCount(std::size_t cnt)
    {
        auto ts = std::chrono::high_resolution_clock::now();

        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cnt = cnt < m_total ? cnt : m_total;
            m_prevTs = ts;
        }

        m_wakeUp.notify_all();
    }

    void CLIProgBar::IncrementCount()
    {
        auto ts = std::chrono::high_resolution_clock::now();

        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cnt = m_cnt + 1 < m_total ? m_cnt + 1 : m_total;
            m_prevTs = ts;
        }

        m_wakeUp.notify_all();
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

    std::string CLIProgBar::GetProgressText(std::size_t cnt, bool percentage) const
    {
        std::string str;
        str.reserve(16);
        str.push_back(' ');

        if (percentage)
        {
            std::size_t prog = (100 * cnt) / m_total;
            str.append(std::to_string(prog));
            str.insert(str.begin() + 1, std::min<std::size_t>(4 - str.size(), 4), ' ');
            str.push_back('%');
        }
        else
        {
            std::string tot = std::to_string(m_total);
            str.append(std::to_string(cnt));
            str.insert(str.begin() + 1, std::min<std::size_t>(tot.size() + 1 - str.size(), tot.size() + 1), ' ');
            str.push_back('/');
            str.append(tot);
        }

        return str;
    }

    std::string CLIProgBar::GetTimeText(std::size_t cnt, Timestamp_t cntTs) const
    {
        auto ts = std::chrono::high_resolution_clock::now();
        double elapsedSecs = (ts - m_startTs).count() / 1e9;

        // 14 character string
        std::string str;
        str.reserve(32);
        str.append(" [");

        if (cnt < m_total)
        {
            // estimate total time
            double estSecs = 0.0;
            if (cnt > 0)
                estSecs = ((cntTs - m_startTs).count() / 1e9) * m_total / cnt;

            str.append(TimeToString(elapsedSecs, false));
            str.push_back('<');
            str.append(TimeToString(estSecs, false));
        }
        else
        {
            str.append(TimeToString(elapsedSecs, true));
        }
        str.push_back(']');

        return str;
    }

    std::string CLIProgBar::GenerateBar() const
    {
        std::string str;
        str.reserve(m_width);
        
        // retrieve state
        m_mutex.lock();
        std::size_t cnt = m_cnt < m_total ? m_cnt : m_total;
        Timestamp_t cntTs = m_prevTs;
        m_mutex.unlock();

        std::string front = m_title;
        std::string back = GetProgressText(cnt, true) + GetTimeText(cnt, cntTs);

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
        std::size_t iProg = barWidth * cnt / m_total;
        str.append(iProg, m_progressChar);
        str.append(barWidth - iProg, ' ');
        str.append("]");
        str.append(back);

        return str;
    }


    void CLIProgBar::WorkerThread()
    {
        std::string printedStr;

        while (true)
        {
            std::string currBar;

            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if(m_cnt >= m_total)
                    m_stopThread = true; // Print one more time then quit
            }

            std::string bar = GenerateBar();
            if (bar != currBar)
            {
                std::cout << "\u001b[1000D" << bar << std::flush;
                currBar = bar;
            }

            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if (m_stopThread) break;
                m_wakeUp.wait_for(lock, 0.25s);
                if (m_stopThread) break;
            }
        }

        std::cout << std::endl;
    }

}
