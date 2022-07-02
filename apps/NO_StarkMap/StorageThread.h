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
    using Data_t = std::tuple<int, double, 
        Eigen::VectorXd, Eigen::MatrixXd, 
        Eigen::Matrix<double, Eigen::Dynamic, 4>>;
public:
    StorageThread(const std::string& path, 
        const QSim::RydbergDiatomicState_t& state, 
        const std::vector<QSim::RydbergDiatomicState_t>& basis, 
        double dE, std::size_t cnt);
    ~StorageThread();

    void StoreData(int i, double eField, 
        const Eigen::VectorXd& energies, 
        const Eigen::MatrixXd& states, 
        const Eigen::Matrix<double, Eigen::Dynamic, 4>& character);
    Data_t LoadData(int i);
    void AdjustStates(int i, const Eigen::MatrixXd& states);
    void WaitUntilFinished();

private:
    void ThreadProc();
    std::string GenerateGroupName(int i) const;

private:
    std::string m_path;

    std::thread m_thread;
    std::mutex m_mutexQueue, m_mutexFile;
    std::condition_variable m_cond;
    std::queue<Data_t> m_dataQueue;

    const std::size_t m_totalCnt;
};
