// Philipp Neufeld, 2021-2022

#include <iostream>
#include <thread>
#include <queue>
#include <condition_variable>
#include <future>

#include <Eigen/Dense>

#include <QSim/Util/DataFile.h>
#include <QSim/Rydberg/RydbergDiatomic.h>

class IOThread
{
    using Data_t = std::tuple<int, double, 
        Eigen::VectorXd, Eigen::MatrixXd, 
        Eigen::Matrix<double, Eigen::Dynamic, 4>>;
public:
    IOThread(const std::string& path);
    ~IOThread();

    bool Start(const std::vector<QSim::RydbergDiatomicState_t>& basis, 
        double energyCenter, double dE);
    void Stop();

    void StoreData(int i, double eField,
        const Eigen::VectorXd& energies, 
        const Eigen::MatrixXd& states, 
        const Eigen::Matrix<double, Eigen::Dynamic, 4>& character);
    std::future<Data_t> LoadData(int i);

    void WaitUntilFinished();

private:
    void ThreadProc();
    std::string GenerateGroupName(int i) const;

private:
    std::string m_path;

    std::thread m_thread;
 
    bool m_stopThread;
    enum class CmdType { Load, Store };
    std::queue<CmdType> m_commandQueue;
    std::queue<Data_t> m_storageQueue;
    std::queue<std::pair<int, std::promise<Data_t>>> m_loadingQueue;
    std::mutex m_mutex;
    std::condition_variable m_cond;
};
