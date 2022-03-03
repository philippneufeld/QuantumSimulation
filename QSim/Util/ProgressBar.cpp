// Philipp Neufeld, 2021-2022

#include "../Platform.h"
#include "ProgressBar.h"

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

    ProgressBar::ProgressBar(std::size_t total, const std::string& title)
        : m_stopThread(false), 
        m_startTs(std::chrono::high_resolution_clock::now()), 
        m_total(total), m_cnt(0), 
        m_title(title), m_width(GetCLIWidth()), m_progressChar('#') 
        { 
            m_currTs = m_startTs;
            m_thread = std::thread([&]() { WorkerThread(); });
        }

    ProgressBar::~ProgressBar()
    {
        if (m_thread.joinable())
            m_thread.join();
    }

    void ProgressBar::WaitUntilFinished()
    {
        std::unique_lock<std::mutex> lock(m_mutex);  
        m_wakeUp.wait(lock, [&](){ return m_cnt == m_total; });
        lock.unlock();

        if (m_thread.joinable())
            m_thread.join();
    }

    void ProgressBar::Update()
    {
        m_wakeUp.notify_all();
    }

    void ProgressBar::IncrementCount(std::size_t inc)
    {
        auto ts = std::chrono::high_resolution_clock::now();
        
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cnt = inc < m_total - m_cnt ? m_cnt + inc : m_total;
        m_currTs = ts;
        lock.unlock();

        m_wakeUp.notify_all();
    }

    std::string ProgressBar::TimeToString(double secs, bool millis)
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

    std::string ProgressBar::GetProgressText(std::size_t cnt, std::size_t tot)
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
    std::string ProgressBar::GetProgressTextFrac(std::size_t cnt, std::size_t tot)
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

    std::string ProgressBar::GetTimeText(double elapsedTime, double remainingTime, bool finished)
    {  
        std::string str;
        str.reserve(32);
        str.append(" [");

        if (!finished)
        {
            str.append(TimeToString(elapsedTime, false));
            str.push_back('<');
            str.append(TimeToString(remainingTime, false));
        }
        else
        {
            str.append(TimeToString(elapsedTime, true));
        }
        str.push_back(']');

        return str;
    }

    std::string ProgressBar::GenerateBar(std::size_t width, std::size_t cnt, std::size_t tot, 
            double elapsedTime, double remainigTime, bool printEst) const
    {
        std::string str;
        str.reserve(width);
        
        std::string front = m_title;
        std::string back = GetProgressText(cnt, tot) + GetTimeText(elapsedTime, remainigTime, !printEst);

        // cutoff front and back string if necessary
        std::size_t minWidth = 5;
        width = std::max(width, minWidth);

        if (front.size() + back.size() > width - minWidth)
        {
            std::size_t exess = std::min(back.size(), front.size() + back.size() - width + minWidth);
            back.erase(back.end() - exess, back.end());

            exess = std::min(front.size(), front.size() + back.size() - width + minWidth);
            front.erase(front.end() - exess, front.end());
        }

        str.append(front);
        str.append("[");
        std::size_t barWidth = width - std::min(width - 2, front.size() + back.size()) - 2;
        std::size_t iProg = barWidth * cnt / tot;
        str.append(iProg, m_progressChar);
        str.append(barWidth - iProg, ' ');
        str.append("]");
        str.append(back);

        return str;
    }

    void ProgressBar::WorkerThread()
    {
        std::string currBar;
        
        std::unique_lock<std::mutex> lock(m_mutex);

        while (!m_stopThread)
        {
            // store counts to local variables
            std::size_t cnt = m_cnt;
            std::size_t tot = m_total;
            std::size_t width = m_width;

            // check exit condition (finish the loop iteration anyways)
            if (cnt >= tot)
                m_stopThread = true;

            // calculate elapsed and estimated total time
            auto ts = std::chrono::high_resolution_clock::now();
            double currSpan = (m_currTs - m_startTs).count() / 1e9;
            double elT = (cnt == tot ? currSpan : (ts - m_startTs).count() / 1e9);
            double estTot = (cnt > 0 ? currSpan * tot / cnt : 2 * currSpan * tot);
            
            // release the lock (needs not to be held during bar string creation)
            lock.unlock();

            // generate new bar string and print it if it differs 
            // from the previously printed string
            std::string bar = GenerateBar(width, cnt, tot, elT, estTot - elT, cnt < tot);
            if (bar != currBar)
            {
                std::cout << "\u001b[1000D" << bar << std::flush;
                currBar = bar;
            }

            // reaquire lock
            lock.lock(); 
            
            if (m_stopThread)
                break;
            m_wakeUp.wait_for(lock, 100ms);
        }

        std::cout << std::endl;
    }

}
