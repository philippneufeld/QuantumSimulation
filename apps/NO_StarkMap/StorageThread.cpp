// Philipp Neufeld, 2021-2022

#include <array>
#include "StorageThread.h"

using namespace QSim;
using namespace Eigen;

Matrix<double, 5, 1> StateToArray(const RydbergDiatomicState_t& state)
{
    auto createArray = [](auto&&... x) { return std::array<double, 5>{static_cast<double>(x)...}; };
    return Map<Matrix<double, 5, 1>>(std::apply(createArray, state).data());
}

StorageThread::StorageThread(const std::string& path, 
    const RydbergDiatomicState_t& state, const std::vector<RydbergDiatomicState_t>& basis, 
    double dE, std::size_t cnt)
    : m_totalCnt(cnt)
{
    m_file.Open(path, DataFile_MUST_NOT_EXIST);
    auto root = m_file.OpenRootGroup();
    
    // write stark map parameters
    root.SetAttribute("State", StateToArray(state));
    root.SetAttribute("Energy_Range", dE);

    // store basis
    MatrixXd basisMat(basis.size(), 5);
    for (int i = 0; i < basis.size(); i++)
        basisMat.row(i) = StateToArray(basis[i]);
    root.CreateDataset("Basis", basisMat);

    m_thread = std::thread([&](){ this->ThreadProc(); });
}

StorageThread::~StorageThread() 
{ 
    WaitUntilFinished();
    m_file.Close(); 
}

void StorageThread::AddData(int i, double eField, const VectorXd& energies, 
    const MatrixXd& states, const Matrix<double, Dynamic, 4>& character)
{
    std::unique_lock<std::mutex> lock(m_mutex);
    m_dataQueue.emplace(i, eField, energies, states, character);
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
            const auto& [i, eField, energies, states, character] = dat;

            // create data group
            std::string groupName = std::to_string(i);
            groupName.insert(0, groupNameLen - groupName.size(), '0');
            auto group = root.CreateSubgroup(groupName);

            // store data
            group.SetAttribute("Electric_Field", eField);
            group.CreateDataset("Energies", energies);
            group.CreateDataset("States", states);
            group.CreateDataset("Character", character);

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

