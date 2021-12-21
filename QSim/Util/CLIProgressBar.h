// Philipp Neufeld, 2021

#ifndef QSim_CLIProgressBar_H_
#define QSim_CLIProgressBar_H_

#include "../Platform.h"

#include <string>
#include <iostream>
#include <atomic>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <chrono>

namespace QSim
{

    std::size_t GetCLIWidth();

    class CLIProgBar
    {
    public:
        CLIProgBar();
        virtual ~CLIProgBar();

        void SetWidth(std::size_t width);
        void SetProgressCharacter(char progChar);
        void SetTitle(const std::string& title);

        double GetProgress() const;
        void SetProgress(double prog);
        
        void Start();
        void WaitUntilFinished();
        void Stop();
        void Update();

        virtual std::string GetDescription() const;
        virtual std::string GetProgressText() const;
        virtual std::string GetTimeText() const;

    private:
        static std::string SecondsToString(std::size_t secs);
        void ThreadFunc();

    private:
        mutable std::mutex m_mutex;

        std::thread m_thread;
        std::condition_variable m_wakeUp;
        bool m_stopThread;
        std::chrono::high_resolution_clock::time_point m_startTs;
        std::chrono::high_resolution_clock::time_point m_prevTs;

        std::string m_title = "";
        
        char m_progressChar;
        bool m_started;
        double m_progress;
        std::size_t m_width;
    };


    class CLIProgBarInt : public CLIProgBar
    {
    public:
        CLIProgBarInt();
        CLIProgBarInt(std::size_t total);

        virtual std::string GetProgressText() const override;

        void SetTotal(std::size_t total);
        void SetCount(std::size_t cnt);

        void IncrementCount();

    private:
        std::atomic<std::size_t> m_total;
        std::atomic<std::size_t> m_cnt;
    };

}

#endif
