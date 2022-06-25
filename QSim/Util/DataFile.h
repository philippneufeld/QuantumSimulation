// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_DataFile_H_
#define QSim_Util_DataFile_H_

#include <string>
#include <vector>
#include <complex>
#include <hdf5.h>

#include <Eigen/Dense>
#include "../Math/MatrixTraits.h"
#include "ScopeGuard.h"


namespace QSim
{

    //
    // DataFileType
    //

    template<typename Ty>
    class DataFileType;
    template<typename Ty>
    class DataFileType<std::complex<Ty>>
    {
    public:
        DataFileType() 
        {
            auto dft = DataFileType<Ty>{};
            m_msTy = dft.GetStorageType();
            m_mnTy = dft.GetNativeType();
            H5Iinc_ref(m_msTy);
            H5Iinc_ref(m_mnTy);

            m_sTy = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<Ty>));
            H5Tinsert(m_sTy, "real", 0, m_msTy);
            H5Tinsert(m_sTy, "imaginary", sizeof(Ty), m_msTy);

            m_nTy = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<Ty>));
            H5Tinsert(m_nTy, "real", 0, m_mnTy);
            H5Tinsert(m_nTy, "imaginary", sizeof(Ty), m_mnTy);
        }

        ~DataFileType()
        {
            H5Tclose(m_sTy);
            H5Tclose(m_nTy);
            H5Idec_ref(m_msTy);
            H5Idec_ref(m_mnTy);
        }

        hid_t GetStorageType() const { return m_sTy; }
        hid_t GetNativeType() const { return m_nTy; }

    private:
        hid_t m_msTy, m_mnTy;
        hid_t m_sTy, m_nTy;
    };

    template<>
    class DataFileType<double>
    {
    public:
        hid_t GetStorageType() const { return H5T_IEEE_F64LE; }
        hid_t GetNativeType() const { return H5T_NATIVE_DOUBLE; }
    };


    //
    // DataFileObject
    //

    class DataFileObject
    {
    protected:
        // SFINAE switch for GetAttribute function
        // match vectors and and row-major matrices for SO=RowMajor
        // match col-major matrices for SO=ColMajor
        template<typename Ty, int SO, typename Ret=Ty>
        struct EnableIfMat;
        template<typename Ty, int N, int M, int NM, int MM, int SO, typename Ret>
        struct EnableIfMat<Eigen::Matrix<Ty, N, M, SO, NM, MM>, Eigen::RowMajor, Ret>
            : std::enable_if<SO == Eigen::RowMajor || N == 1 || M == 1, Ret> {};
        template<typename Ty, int N, int M, int NM, int MM, int SO, typename Ret>
        struct EnableIfMat<Eigen::Matrix<Ty, N, M, SO, NM, MM>, Eigen::ColMajor, Ret>
            : std::enable_if<!(SO == Eigen::RowMajor || N == 1 || M == 1), Ret> {};
        template<typename Ty, int SO, typename Ret=Ty>
        using EnableIfMat_t = typename EnableIfMat<Ty, SO, Ret>::type;

        template<typename Ty, typename Ret>
        struct EnableIfScalar : std::enable_if<true, Ret> {};
        template <typename Ty, int N, int M, int NM, int MM, typename Ret>
        struct EnableIfScalar<Eigen::Matrix<Ty, N, M, NM, MM>, Ret>;
        template<typename Ty, typename Ret=Ty>
        using EnableIfScalar_t = typename EnableIfScalar<Ty, Ret>::type;

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

        template<typename Ty>
        EnableIfScalar_t<Ty> GetAttribute(const std::string& name) const;
        template<typename Ty>
        EnableIfMat_t<Ty, Eigen::RowMajor> GetAttribute(const std::string& name) const;
        template<typename Ty>
        EnableIfMat_t<Ty, Eigen::ColMajor> GetAttribute(const std::string& name) const;
        template<typename Ty>
        bool GetAttributeRaw(const std::string& name, Ty* data) const;
        template<typename Ty>
        bool GetAttributeRawChecked(const std::string& name, const std::vector<std::size_t>& dims, Ty* data) const;
        
        template<typename Ty>
        EnableIfScalar_t<Ty, bool> SetAttribute(const std::string& name, Ty value);
        template<typename Ty>
        EnableIfMat_t<Ty, Eigen::RowMajor, bool> SetAttribute(const std::string& name, const Ty& value);
        template<typename Ty>
        EnableIfMat_t<Ty, Eigen::ColMajor, bool> SetAttribute(const std::string& name, const Ty& value);
        template<typename Ty>
        bool SetAttributeRaw(const std::string& name, const std::vector<std::size_t>& dims, const Ty* data);

    protected:
        hid_t GetNative() const { return m_hid; }
        static std::vector<std::size_t> GetAttributeDimsHelper(hid_t attr);

    private:
        hid_t m_hid;
    };

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

        template<typename Ty> std::enable_if_t<std::is_same_v<Ty, double>, Ty> Get() const;
        template<typename Ty> EnableIfMat_t<Ty, Eigen::RowMajor> Get() const;
        template<typename Ty> EnableIfMat_t<Ty, Eigen::ColMajor> Get() const;
        template<typename Ty> bool GetRaw(Ty* data) const;
        template<typename Ty> bool GetRawChecked(const std::vector<std::size_t>& dims, Ty* data) const;
        
        template<typename Ty> EnableIfScalar_t<Ty, bool> Set(Ty value);
        template<typename Ty> EnableIfMat_t<Ty, Eigen::RowMajor, bool> Set(const Ty& value);
        template<typename Ty> EnableIfMat_t<Ty, Eigen::ColMajor, bool> Set(const Ty& value);
        template<typename Ty> bool SetRaw(const std::vector<std::size_t>& dims, const Ty* data);
    };

    //
    // DataFileGroup
    //

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
        template<typename Ty> 
        EnableIfScalar_t<Ty, DataFileDataset> CreateDataset(const std::string& name, Ty value);
        template<typename Ty> 
        EnableIfMat_t<Ty, Eigen::RowMajor, DataFileDataset> CreateDataset(const std::string& name, const Ty& value);
        template<typename Ty> 
        EnableIfMat_t<Ty, Eigen::ColMajor, DataFileDataset> CreateDataset(const std::string& name, const Ty& value);
        template<typename Ty> 
        DataFileDataset CreateDatasetRaw(const std::string& name, const std::vector<std::size_t>& dims, const Ty* data);

        bool DoesDatasetExist(const std::string& name) const;
        DataFileDataset GetDataset(const std::string& name) const;
        std::vector<std::string> EnumerateDatasets() const;

    private:
        static herr_t EnumGroupsHelper(hid_t group, const char *name, const H5L_info_t *info, void *data);
        static herr_t EnumDatasetHelper(hid_t group, const char *name, const H5L_info_t *info, void *data);
    };


    //
    // DataFile
    //

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


    //
    // Template function definitions (DataFileObject)
    //

    template<typename Ty>
    DataFileObject::EnableIfScalar_t<Ty>
        DataFileObject::GetAttribute(const std::string& name) const
    {
        Ty val = 0.0;
        if (!GetAttributeRawChecked(name, {1}, &val))
            return std::numeric_limits<Ty>::quiet_NaN();
        return val;
    }

    template<typename Ty>
    DataFileObject::EnableIfMat_t<Ty, Eigen::RowMajor> 
        DataFileObject::GetAttribute(const std::string& name) const
    {
        auto dims = GetAttributeDims(name);
        if (dims.size() != 2 || dims[0]*dims[1] == 0)
            return Ty{};

        Ty mat(dims[0], dims[1]);
        return (GetAttributeRaw(name, &(mat(0, 0))) ? mat : Ty{});
    }

    template<typename Ty>
    DataFileObject::EnableIfMat_t<Ty, Eigen::ColMajor> 
        DataFileObject::GetAttribute(const std::string& name) const
    {
        return this->template GetAttribute<MakeRowMajor_t<Ty>>(name);
    }

    template<typename Ty>
    bool DataFileObject::GetAttributeRaw(const std::string& name, Ty* data) const
    {
        hid_t attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
        if (attr < 0) return false;
        auto attrGuard = CreateScopeGuard([=](){ H5Aclose(attr); });
        
        auto type = DataFileType<std::decay_t<Ty>>{};
        if (H5Aread(attr, type.GetNativeType(), data) < 0)
            return false;
        
        return true;
    }

    template<typename Ty>
    bool DataFileObject::GetAttributeRawChecked(const std::string& name, 
        const std::vector<std::size_t>& dims, Ty* data) const
    {
        hid_t attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
        if (attr < 0) return false;
        auto attrGuard = CreateScopeGuard([=](){ H5Aclose(attr); });
        
        if (GetAttributeDimsHelper(attr) != dims)
                return false;

        auto type = DataFileType<std::decay_t<Ty>>{};
        if (H5Aread(attr, type.GetNativeType(), data) < 0)
            return false;
        
        return true;
    }

    template<typename Ty>
    DataFileObject::EnableIfScalar_t<Ty, bool> 
        DataFileObject::SetAttribute(const std::string& name, Ty value)
    {
        return SetAttributeRaw(name, {1, }, &value);
    }

    template<typename Ty>
    DataFileObject::EnableIfMat_t<Ty, Eigen::RowMajor, bool> 
        DataFileObject::SetAttribute(const std::string& name, const Ty& value)
    {
        // data has to be in row-major format to be stored
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return SetAttributeRaw(name, {rows, cols}, &(value(0,0)));
    }

    template<typename Ty>
    DataFileObject::EnableIfMat_t<Ty, Eigen::ColMajor, bool> 
        DataFileObject::SetAttribute(const std::string& name, const Ty& value)
    {
        return SetAttribute(name, static_cast<const MakeRowMajor_t<Ty>&>(value));
    }

    template<typename Ty>
    bool DataFileObject::SetAttributeRaw(const std::string& name, const std::vector<std::size_t>& dims, const Ty* data)
    {
        hid_t attr = H5I_INVALID_HID;
        auto attrGuard = CreateScopeGuard([&](){ H5Aclose(attr); });

        auto type = DataFileType<std::decay_t<Ty>>{};

        // make sure an attribute with the correct dimensionality is available
        if (!DoesAttributeExist(name))
        {
            // create dataspace
            std::vector<hsize_t> hdims(dims.begin(), dims.end());
            hid_t dspace = H5Screate_simple(hdims.size(), hdims.data(), hdims.data());
            if (dspace < 0) return false;
            auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });

            // create the actual attribute
            attr = H5Acreate2(m_hid, name.c_str(), type.GetStorageType(), dspace, H5P_DEFAULT, H5P_DEFAULT);
            if (attr < 0) return false;
        }
        else
        {
            // open attribute and validate dimensions
            attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
            if (GetAttributeDimsHelper(attr) != dims)
                return false;
        }

        // write data
        if(H5Awrite(attr, type.GetNativeType(), data) < 0)
            return false;

        return true;
    }

    //
    // Template function definitions (DataFileDataset)
    //

    template<typename Ty> 
    std::enable_if_t<std::is_same_v<Ty, double>, Ty> DataFileDataset::Get() const
    {
        double val = 0.0;
        if (!GetRawChecked({1}, &val))
            return std::numeric_limits<double>::quiet_NaN();
        return val;
    }

    template<typename Ty> 
    DataFileDataset::EnableIfMat_t<Ty, Eigen::RowMajor> DataFileDataset::Get() const
    {
        auto dims = GetDims();
        if (dims.size() != 2 || dims[0]*dims[1] == 0)
            return Ty{};

        Ty mat(dims[0], dims[1]);
        return (GetRaw(&(mat(0, 0))) ? mat : Ty{});
    }

    template<typename Ty> 
    DataFileDataset::EnableIfMat_t<Ty, Eigen::ColMajor> DataFileDataset::Get() const
    {
        return this->template Get<MakeRowMajor_t<Ty>>();
    }

    template<typename Ty> 
    bool DataFileDataset::GetRaw(Ty* data) const
    {
        auto type = DataFileType<std::decay_t<Ty>>{};
        return !(H5Dread(GetNative(), type.GetNativeType(), 
            H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0);
    }
    
    template<typename Ty> 
    bool DataFileDataset::GetRawChecked(const std::vector<std::size_t>& dims, Ty* data) const
    {
        // validate dimensions
        if (dims != GetDims())
            return false;

        return GetRaw(data);
    }

    template<typename Ty>
    DataFileDataset::EnableIfScalar_t<Ty, bool> DataFileDataset::Set(Ty value)
    {
        return SetRaw({1}, &value);
    }

    template<typename Ty>
    DataFileDataset::EnableIfMat_t<Ty, Eigen::RowMajor, bool> DataFileDataset::Set(const Ty& value)
    {
        // data has to be in row-major format to be stored
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return SetRaw({rows, cols}, &(value(0,0)));
    }

    template<typename Ty>
    DataFileDataset::EnableIfMat_t<Ty, Eigen::ColMajor, bool> DataFileDataset::Set(const Ty& value)
    {
        return Set(static_cast<const MakeRowMajor_t<Ty>&>(value));
    }

    template<typename Ty>
    bool DataFileDataset::SetRaw(const std::vector<std::size_t>& dims, const Ty* data)
    {
        auto type = DataFileType<std::decay_t<Ty>>{};

        // validate dimensions
        if (dims != GetDims())
            return false;

        if(H5Dwrite(GetNative(), type.GetNativeType(), H5S_ALL, 
                    H5S_ALL, H5P_DEFAULT, data) < 0)
            return false;

        if(H5Dflush(GetNative()) > 0)
            return false;
 
        return true;
    }

    
    //
    // Template function definitions (DataFileGroup)
    //

    template<typename Ty>
    DataFileGroup::EnableIfScalar_t<Ty, DataFileDataset> 
        DataFileGroup::CreateDataset(const std::string& name, Ty value)
    {
        return CreateDatasetRaw(name, {1}, &value);
    }

    template<typename Ty> 
    DataFileGroup::EnableIfMat_t<Ty, Eigen::RowMajor, DataFileDataset> 
        DataFileGroup::CreateDataset(const std::string& name, const Ty& value)
    {
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return CreateDatasetRaw(name, {rows, cols}, &(value(0,0)));
    }

    template<typename Ty> 
    DataFileGroup::EnableIfMat_t<Ty, Eigen::ColMajor, DataFileDataset> 
        DataFileGroup::CreateDataset(const std::string& name, const Ty& value)
    {
        return CreateDataset(name, static_cast<const MakeRowMajor_t<Ty>&>(value));
    }

    template<typename Ty>
    DataFileDataset DataFileGroup::CreateDatasetRaw(const std::string& name, 
        const std::vector<std::size_t>& dims, const Ty* data)
    {
        if (DoesDatasetExist(name))
            return H5I_INVALID_HID;

        auto type = DataFileType<std::decay_t<Ty>>{};

        std::vector<hsize_t> hdims(dims.begin(), dims.end());
        hid_t dspace = H5Screate_simple(hdims.size(), hdims.data(), hdims.data());
        if (dspace < 0) return H5I_INVALID_HID;
        auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });
        
        DataFileDataset dataset = H5Dcreate2(GetNative(), name.c_str(), 
            type.GetStorageType(), dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        if (!dataset) 
            return H5I_INVALID_HID;

        if (!dataset.SetRaw(dims, data))
            return H5I_INVALID_HID;

        return dataset;
    }



}

#endif
