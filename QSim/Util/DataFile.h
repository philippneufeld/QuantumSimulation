// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_DataFile_H_
#define QSim_Util_DataFile_H_

#include <string>
#include <vector>
#include <hdf5.h>

#include "../Math/Matrix.h"

namespace QSim
{

    class DataFileObject
    {
    protected:
        DataFileObject(hid_t hid);

    public:
        DataFileObject();
        virtual ~DataFileObject();

        DataFileObject(const DataFileObject&);
        DataFileObject(DataFileObject&& rhs);
        DataFileObject& operator=(const DataFileObject&);
        DataFileObject& operator=(DataFileObject&& rhs);

        bool IsValid() const;
        operator bool() const { return IsValid(); }
        bool operator!() const { return !IsValid(); }

        // attributes
        bool DoesAttributeExist(const std::string& name) const;
        bool CreateAttribute(const std::string& name, const std::vector<std::size_t>& dims);
        std::vector<std::size_t> GetAttributeDims(const std::string& name) const;
        bool LoadAttribute(const std::string& name, double* data) const;
        bool StoreAttribute(const std::string& name, const double* data);

    protected:
        hid_t GetNative() const { return m_hid; }

    private:
        hid_t m_hid;
    };


    class DataFileDataset : public DataFileObject
    {
        friend class DataFileGroup;
    protected:
        DataFileDataset(hid_t dataset);

    public:
        DataFileDataset() = default;

        std::vector<std::size_t> GetDims() const;
        bool Load(double* data) const;
        bool Store(const double* data);

        TDynamicMatrix<double> Load2DMatrix() const;
        TDynamicRowVector<double> Load1DRowVector() const;
        TDynamicColVector<double> Load1DColVector() const;
    };

    class DataFileGroup : public DataFileObject
    {
        friend class DataFile;
    protected:
        DataFileGroup(hid_t group);

    public:
        DataFileGroup() = default;

        // subgroups
        bool DoesSubgroupExist(const std::string& name) const;
        DataFileGroup CreateSubgroup(const std::string& name);
        DataFileGroup GetSubgroup(const std::string& name) const;
        std::vector<std::string> EnumerateSubgroups() const;
        
        // datasets
        bool DoesDatasetExist(const std::string& name) const;
        DataFileDataset CreateDataset(const std::string& name, const std::vector<std::size_t>& dims);
        DataFileDataset GetDataset(const std::string& name) const;
        std::vector<std::string> EnumerateDatasets() const;
        
    private:
        static herr_t EnumGroupsHelper(hid_t group, const char *name, const H5L_info_t *info, void *data);
        static herr_t EnumDatasetHelper(hid_t group, const char *name, const H5L_info_t *info, void *data);
    };


    enum DataFileOpenFlag
    {
        DataFile_DEFAULT = 0,
        DataFile_MUST_EXIST = 1, 
        DataFile_MUST_NOT_EXIST = 2,
        DataFile_TRUNCATE = 4,
    };

    class DataFile
    {
    public:
        DataFile();
        ~DataFile();

        DataFile(const DataFile&) = delete;
        DataFile(DataFile&& rhs);
        DataFile& operator=(const DataFile&) = delete;
        DataFile& operator=(DataFile&& rhs);

        bool Open(const std::string& name, DataFileOpenFlag flag);
        void Close();

        bool IsOpen() const;
        operator bool() const { return IsOpen(); }
        bool operator!() const { return !IsOpen(); }

        DataFileGroup OpenRootGroup();

    private:
        hid_t m_file;
    };

}

#endif
