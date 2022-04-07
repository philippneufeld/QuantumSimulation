// Philipp Neufeld, 2021-2022

#include <iostream>
#include <thread>
#include <queue>
#include <condition_variable>

#include <Eigen/Dense>

#include <QSim/Util/DataFile.h>
#include <QSim/Rydberg/RydbergDiatomic.h>

class StorageThread
{
    using Data_t = std::tuple<int, double, Eigen::VectorXd, Eigen::MatrixXd>;
public:
    StorageThread(const std::string& path, 
        const QSim::RydbergDiatomicState_t& state, 
        const std::vector<QSim::RydbergDiatomicState_t>& basis, 
        double dE, std::size_t cnt);
    ~StorageThread();

    void AddData(int i, double eField, const Eigen::VectorXd& energies, const Eigen::MatrixXd& states);
    void WaitUntilFinished();

private:
    void ThreadProc();

private:
    QSim::DataFile m_file;

    std::thread m_thread;
    std::mutex m_mutex;
    std::condition_variable m_cond;
    std::queue<Data_t> m_dataQueue;

    const std::size_t m_totalCnt;
};
