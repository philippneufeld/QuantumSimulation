// Philipp Neufeld, 2021-2022

#include <array>
#include "StorageThread.h"

using namespace QSim;
using namespace Eigen;

std::array<double, 5> StateToArray(const RydbergDiatomicState_t& state)
{
    auto createArray = [](auto&&... x) { return std::array<double, 5>{static_cast<double>(x)...}; };
    return std::apply(createArray, state);
}

StorageThread::StorageThread(const std::string& path, const RydbergDiatomicState_t& state, const std::vector<RydbergDiatomicState_t>& basis, double dE, std::size_t cnt)
    : m_totalCnt(cnt)
{
    m_file.Open(path, DataFile_MUST_NOT_EXIST);
    auto root = m_file.OpenRootGroup();
    
    // write stark map parameters
    root.CreateAttribute("State", { 5 });
    root.StoreAttribute("State", StateToArray(state).data());
    root.CreateAttribute("Energy_Range", { 1 });
    root.StoreAttribute("Energy_Range", &dE);

    // store basis
    auto basisStorage = root.CreateDataset("basis", {basis.size(), 5});
    MatrixXd basisMat(basis.size(), 5);
    for (int i = 0; i < basis.size(); i++)
        basisMat.row(i) = Map<Matrix<double, 5, 1>>(StateToArray(state).data());
    basisStorage.StoreMatrix(basisMat);

    m_thread = std::thread([&](){ this->ThreadProc(); });
}

StorageThread::~StorageThread() 
{ 
    WaitUntilFinished();
    m_file.Close(); 
}

void StorageThread::AddData(int i, double eField, const VectorXd& energies, const MatrixXd& states)
{
    std::unique_lock<std::mutex> lock(m_mutex);
    m_dataQueue.emplace(i, eField, energies, states);
    m_cond.notify_all();
}

void StorageThread::WaitUntilFinished() 
{
    m_cond.notify_all();
    if (m_thread.joinable())
        m_thread.join();
}

void StorageThread::ThreadProc()
{
    int cnt = 0;
    int groupNameLen = std::to_string(m_totalCnt).size();
    auto root = m_file.OpenRootGroup();

    std::unique_lock<std::mutex> lock(m_mutex);
    while(true)
    {
        m_cond.wait(lock, [&](){ return !m_dataQueue.empty() || cnt >= m_totalCnt; });
        
        // move to local vector (this way the queue lock can be released)
        std::vector<Data_t> data;
        data.reserve(m_dataQueue.size());
        while (!m_dataQueue.empty())
        {
            data.push_back(std::move(m_dataQueue.front()));
            m_dataQueue.pop();
        }
        lock.unlock();

        // store data to file
        for (const Data_t& dat : data)
        {
            const auto& [i, eField, energies, states] = dat;

            // create data group
            std::string groupName = std::to_string(i);
            groupName.insert(0, groupNameLen - groupName.size(), '0');
            auto group = root.CreateSubgroup(groupName);

            // define datasets
            group.CreateAttribute("Electric_Field", { 1 });
            auto energyStorage = group.CreateDataset("Energies", { static_cast<std::size_t>(energies.size()) });
            auto stateStorage = group.CreateDataset("States", { static_cast<std::size_t>(states.rows()), 
                static_cast<std::size_t>(states.cols()) });

            // store data
            group.StoreAttribute("Electric_Field", &eField);
            energyStorage.Store(energies.data());
            stateStorage.StoreMatrix(states);

            cnt++;
        }

        if (cnt >= m_totalCnt)
            break;

        // reaquire lock
        lock.lock();
    }

    // finished writing all data
    m_file.Close();
}

