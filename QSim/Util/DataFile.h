// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_DataFile_H_
#define QSim_Util_DataFile_H_

#include <string>
#include <vector>
#include <hdf5.h>

#include <Eigen/Dense>
#include "../Math/MatrixTraits.h"

namespace QSim
{

    //
    // DataFileObject
    //

    class DataFileObject
    {
        // SFINAE switch for GetAttribute function
        // match vectors and and row-major matrices for SO=RowMajor
        // match col-major matrices for SO=ColMajor
        template<typename Ty, int SO>
        struct EnableIfGetAttr;
        template<typename Ty, int N, int M, int NM, int MM, int SO>
        struct EnableIfGetAttr<Eigen::Matrix<Ty, N, M, SO, NM, MM>, Eigen::RowMajor>
            : std::enable_if<SO == Eigen::RowMajor || N == 1 || M == 1, Eigen::Matrix<Ty, N, M, SO, NM, MM>> {};
        template<typename Ty, int N, int M, int NM, int MM, int SO>
        struct EnableIfGetAttr<Eigen::Matrix<Ty, N, M, SO, NM, MM>, Eigen::ColMajor>
            : std::enable_if<!(SO == Eigen::RowMajor || N == 1 || M == 1), Eigen::Matrix<Ty, N, M, SO, NM, MM>> {};
        template<typename Ty, int SO>
        using EnableIfGetAttr_t = typename EnableIfGetAttr<Ty, SO>::type;

        // conversion of matrix type to row major storage order
        template<typename Ty> 
        struct MakeRowMajor
        {
            using type = Ty;
        };
        template<typename Ty, int N, int M, int NM, int MM> 
        struct MakeRowMajor<Eigen::Matrix<Ty, N, M, Eigen::ColMajor, NM, MM>>
        {
            using type = Eigen::Matrix<Ty, N, M, Eigen::RowMajor, NM, MM>;
        };
        template<typename Ty> 
        using MakeRowMajor_t = typename MakeRowMajor<Ty>::type;

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
        std::vector<std::size_t> GetAttributeDims(const std::string& name) const;

        bool GetAttributeRaw(const std::string& name, double* data) const;
        bool GetAttributeRawChecked(const std::string& name, const std::vector<std::size_t>& dims, double* data) const;
        template<typename Ty>
        std::enable_if_t<std::is_convertible_v<double, Ty>, Ty> GetAttribute(const std::string& name);
        template<typename Ty>
        EnableIfGetAttr_t<Ty, Eigen::RowMajor> GetAttribute(const std::string& name);
        template<typename Ty>
        EnableIfGetAttr_t<Ty, Eigen::ColMajor> GetAttribute(const std::string& name);
        
        bool SetAttribute(const std::string& name, double value);
        bool SetAttribute(const std::string& name, const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>& value);
        bool SetAttributeRaw(const std::string& name, const std::vector<std::size_t>& dims, const double* data);

        bool CreateAttribute(const std::string& name, const std::vector<std::size_t>& dims);
        bool LoadAttribute(const std::string& name, double* data) const;
        bool StoreAttribute(const std::string& name, const double* data);

    protected:
        hid_t GetNative() const { return m_hid; }

        static std::vector<std::size_t> GetAttributeDimsHelper(hid_t attr);

    private:
        hid_t m_hid;
    };

    template<typename Ty>
    std::enable_if_t<std::is_convertible_v<double, Ty>, Ty> 
        DataFileObject::GetAttribute(const std::string& name)
    {
        double val = 0.0;
        if (!GetAttributeRawChecked(name, {1}, &val))
            return std::numeric_limits<double>::quiet_NaN();
        return val;
    }

    template<typename Ty>
    DataFileObject::EnableIfGetAttr_t<Ty, Eigen::RowMajor> 
        DataFileObject::GetAttribute(const std::string& name)
    {
        auto dims = GetAttributeDims(name);
        if (dims.size() != 2 || dims[0]*dims[1] == 0)
            return Ty{};

        Ty mat(dims[0], dims[1]);
        if (!GetAttributeRaw(name, &(mat(0, 0))))
            return Ty{};

        return mat;
    }

    template<typename Ty>
    DataFileObject::EnableIfGetAttr_t<Ty, Eigen::ColMajor> 
        DataFileObject::GetAttribute(const std::string& name)
    {
        return this->template GetAttribute<MakeRowMajor_t<Ty>>(name);
    }


    //
    // DataFileDataset
    //

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

        bool StoreMatrix(Eigen::MatrixXd mat);
        Eigen::MatrixXd LoadMatrix() const;
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
