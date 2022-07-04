// Philipp Neufeld, 2021-2022

#include <array>
#include "IOThread.h"

using namespace QSim;
using namespace Eigen;

Matrix<double, 5, 1> StateToArray(const RydbergDiatomicState_t& state)
{
    auto createArray = [](auto&&... x) { return std::array<double, 5>{static_cast<double>(x)...}; };
    return Map<Matrix<double, 5, 1>>(std::apply(createArray, state).data());
}

IOThread::IOThread(const std::string& path)
    : m_path(path), m_stopThread(false) { }

IOThread::~IOThread() 
{ 
    Stop();
}

bool IOThread::Start(const QSim::RydbergDiatomicState_t& state, 
        const std::vector<QSim::RydbergDiatomicState_t>& basis, double dE)
{
    DataFile file;
    if (!file.Open(m_path, DataFile_DEFAULT))
        return false;
    
    // write stark map parameters
    auto root = file.OpenRootGroup();
    root.SetAttribute("State", StateToArray(state));
    root.SetAttribute("Energy_Range", dE);

    // store basis
    MatrixXd basisMat(basis.size(), 5);
    for (int i = 0; i < basis.size(); i++)
        basisMat.row(i) = StateToArray(basis[i]);
    root.CreateDataset("Basis", basisMat);

    m_stopThread = false;
    m_thread = std::thread([&](){ this->ThreadProc(); });

    return true;
}

void IOThread::Stop()
{
    m_stopThread = true;
    m_cond.notify_all();
    if (m_thread.joinable())
        m_thread.join();
}

void IOThread::StoreData(int i, double eField, const VectorXd& energies, 
    const MatrixXd& states, const Matrix<double, Dynamic, 4>& character)
{
    std::unique_lock<std::mutex> lock(m_mutex);
    m_commandQueue.push(CmdType::Store);
    m_storageQueue.emplace(i, eField, energies, states, character);
    lock.unlock();

    m_cond.notify_all();
}

std::future<typename IOThread::Data_t> IOThread::LoadData(int i)
{
    std::promise<Data_t> promise;
    std::future<Data_t> future = promise.get_future();

    std::unique_lock<std::mutex> lock(m_mutex);
    m_commandQueue.push(CmdType::Load);
    m_loadingQueue.emplace(i, std::move(promise));
    lock.unlock();

    m_cond.notify_all();

    return std::move(future);
}

void IOThread::WaitUntilFinished()
{
    std::unique_lock<std::mutex> lock(m_mutex);
    m_cond.wait(lock, [&](){ return m_commandQueue.empty() || m_stopThread; });
}

void IOThread::ThreadProc()
{
    std::unique_lock<std::mutex> lock(m_mutex);
    while(!m_stopThread)
    {
        m_cond.wait(lock, [&](){ return !m_commandQueue.empty() || m_stopThread; });
        
        if (m_stopThread) break;

        // check what operation to do next
        auto op = m_commandQueue.front();
        m_commandQueue.pop();

        if (op == CmdType::Store)
        {
            // retrieve data from queue
            Data_t data = std::move(m_storageQueue.front());
            m_storageQueue.pop();
            lock.unlock();

            DataFile file;
            file.Open(m_path, DataFile_DEFAULT);
            auto root = file.OpenRootGroup();

            const auto& [i, eField, energies, states, character] = data;
            std::string groupName = GenerateGroupName(i);
            if (!root.DoesSubgroupExist(groupName))
            {
                auto group = root.CreateSubgroup(groupName);
                group.SetAttribute("Electric_Field", eField);
                group.CreateDataset("Energies", energies);
                group.CreateDataset("States", states);
                group.CreateDataset("Character", character);
            }
            else
            {
                auto group = root.GetSubgroup(GenerateGroupName(i)); 
                group.SetAttribute("Electric_Field", eField);
                group.GetDataset("Energies").Set(energies);
                group.GetDataset("States").Set(states);
                group.GetDataset("Character").Set(character);
            }

        }
        else if(op == CmdType::Load)
        {
            auto [i, promise] = std::move(m_loadingQueue.front());
            m_loadingQueue.pop();
            lock.unlock();

            DataFile file;
            file.Open(m_path, DataFile_DEFAULT);
            auto root = file.OpenRootGroup();
            auto group = root.GetSubgroup(GenerateGroupName(i)); 
            auto data = std::make_tuple(
                i, group.GetAttribute<double>("Electric_Field"),
                group.GetDataset("Energies").Get<VectorXd>(),
                group.GetDataset("States").Get<MatrixXd>(),
                group.GetDataset("Character").Get<MatrixXd>());
            
            promise.set_value(std::move(data));
        }

        m_cond.notify_all();

        // reaquire lock
        lock.lock();
    }
}

std::string IOThread::GenerateGroupName(int i) const
{
    int groupNameLen = 9;
    std::string groupName = std::to_string(i);
    groupName.insert(0, groupNameLen - groupName.size(), '0');
    return groupName;
}

