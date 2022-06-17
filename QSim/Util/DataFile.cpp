// Philipp Neufeld, 2021-2022

#include <utility>

#include "DataFile.h"
#include "ScopeGuard.h"

namespace QSim
{

    //
    // DataFileObject
    //

    DataFileObject::DataFileObject()
        : m_hid(H5I_INVALID_HID) { }

    DataFileObject::DataFileObject(hid_t hid) 
        : m_hid(hid) { }
    
    DataFileObject::DataFileObject(const DataFileObject& rhs)
        : m_hid(rhs.m_hid)
    {
        if (m_hid >= 0)
            H5Iinc_ref(m_hid);
    }

    DataFileObject::DataFileObject(DataFileObject&& rhs)
        : m_hid(rhs.m_hid)
    {
        rhs.m_hid = H5I_INVALID_HID;
    }

    DataFileObject& DataFileObject::operator=(const DataFileObject& rhs)
    {
        return this->operator=(std::move(DataFileObject(rhs)));
    }

    DataFileObject& DataFileObject::operator=(DataFileObject&& rhs)
    {
        std::swap(m_hid, rhs.m_hid);
        return *this;
    }

    DataFileObject::~DataFileObject() 
    {
        if (m_hid >= 0)
            H5Idec_ref(m_hid);
        m_hid = H5I_INVALID_HID;
    }

    bool DataFileObject::IsValid() const
    {
        return !(m_hid < 0);
    }

    bool DataFileObject::DoesAttributeExist(const std::string& name) const
    {
        return (H5Aexists(m_hid, name.c_str()) > 0);
    }

    std::vector<std::size_t> DataFileObject::GetAttributeDims(const std::string& name) const
    {
        hid_t attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
        auto attrGuard = CreateScopeGuard([=](){ H5Aclose(attr); });
        return GetAttributeDimsHelper(attr);
    }

    bool DataFileObject::GetAttributeRaw(const std::string& name, double* data) const
    {
        hid_t attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
        if (attr < 0) return false;
        auto attrGuard = CreateScopeGuard([=](){ H5Aclose(attr); });
        
        if (H5Aread(attr, H5T_NATIVE_DOUBLE, data) < 0)
            return false;
        
        return true;
    }

    bool DataFileObject::GetAttributeRawChecked(const std::string& name, 
        const std::vector<std::size_t>& dims, double* data) const
    {
        hid_t attr = H5Aopen(m_hid, name.c_str(), H5P_DEFAULT);
        if (attr < 0) return false;
        auto attrGuard = CreateScopeGuard([=](){ H5Aclose(attr); });
        
        if (GetAttributeDimsHelper(attr) != dims)
                return false;

        if (H5Aread(attr, H5T_NATIVE_DOUBLE, data) < 0)
            return false;
        
        return true;
    }

    bool DataFileObject::SetAttribute(const std::string& name, double value)
    {
        return SetAttributeRaw(name, {1, }, &value);
    }

    bool DataFileObject::SetAttribute(const std::string& name, const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>& value)
    {
        // data has to be in row-major format to be stored
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return SetAttributeRaw(name, {rows, cols}, &(value(0,0)));
    }
    
    bool DataFileObject::SetAttributeRaw(const std::string& name, const std::vector<std::size_t>& dims, const double* data)
    {
        hid_t attr = H5I_INVALID_HID;
        auto attrGuard = CreateScopeGuard([&](){ H5Aclose(attr); });

        // make sure an attribute with the correct dimensionality is available
        if (!DoesAttributeExist(name))
        {
            // create dataspace
            std::vector<hsize_t> hdims(dims.begin(), dims.end());
            hid_t dspace = H5Screate_simple(hdims.size(), hdims.data(), hdims.data());
            if (dspace < 0) return false;
            auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });

            // create the actual attribute
            attr = H5Acreate2(m_hid, name.c_str(), H5T_IEEE_F64LE, dspace, H5P_DEFAULT, H5P_DEFAULT);
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
        if(H5Awrite(attr, H5T_NATIVE_DOUBLE, data) < 0)
            return false;

        return true;
    }

    std::vector<std::size_t> DataFileObject::GetAttributeDimsHelper(hid_t attr)
    {
        if (attr < 0) return std::vector<std::size_t>{};
        
        hid_t dspace = H5Aget_space(attr);
        if (dspace < 0) return std::vector<std::size_t>{};
        auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });

        int ndims = H5Sget_simple_extent_ndims(dspace);
        if (ndims <= 0) return std::vector<std::size_t>{};

        std::vector<hsize_t> dims(ndims);
        if (H5Sget_simple_extent_dims(dspace, &dims[0], nullptr) < 0)
            return std::vector<std::size_t>{};

        return std::vector<std::size_t>(dims.begin(), dims.end());
    }


    //
    // DataFileDataset
    //

    DataFileDataset::DataFileDataset(hid_t dataset)
        : DataFileObject(dataset) { }

    std::vector<std::size_t> DataFileDataset::GetDims() const
    {
        hid_t dspace = H5Dget_space(GetNative());
        if (dspace < 0) return std::vector<std::size_t>{};
        auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });

        int ndims = H5Sget_simple_extent_ndims(dspace);
        if (ndims <= 0) return std::vector<std::size_t>{};

        std::vector<hsize_t> dims(ndims);
        if (H5Sget_simple_extent_dims(dspace, &dims[0], nullptr) < 0)
            return std::vector<std::size_t>{};

        return std::vector<std::size_t>(dims.begin(), dims.end());
    }

    bool DataFileDataset::GetRaw(double* data) const
    {
        return !(H5Dread(GetNative(), H5T_NATIVE_DOUBLE, 
            H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0);
    }
    
    bool DataFileDataset::GetRawChecked(const std::vector<std::size_t>& dims, double* data) const
    {
        // validate dimensions
        if (dims != GetDims())
            return false;

        return GetRaw(data);
    }

    bool DataFileDataset::Set(double value)
    {
        return SetRaw({1}, &value);
    }

    bool DataFileDataset::Set(const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>& value)
    {
        // data has to be in row-major format to be stored
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return SetRaw({rows, cols}, &(value(0,0)));
    }

    bool DataFileDataset::SetRaw(const std::vector<std::size_t>& dims, const double* data)
    {
        // validate dimensions
        if (dims != GetDims())
            return false;

        if(H5Dwrite(GetNative(), H5T_NATIVE_DOUBLE, H5S_ALL, 
                    H5S_ALL, H5P_DEFAULT, data) < 0)
            return false;

        if(H5Dflush(GetNative()) > 0)
            return false;
 
        return true;
    }


    DataFileGroup::DataFileGroup(hid_t group) 
        : DataFileObject(group) { }
    
    bool DataFileGroup::DoesSubgroupExist(const std::string& name) const
    {
        if (H5Lexists(GetNative(), name.c_str(), H5P_DEFAULT) <= 0)
            return false;

        bool valid = false;    
        H5E_BEGIN_TRY { valid = GetSubgroup(name).IsValid(); } H5E_END_TRY;
        return valid;
    }

    DataFileGroup DataFileGroup::CreateSubgroup(const std::string& name)
    {
        if (DoesSubgroupExist(name))
            return H5I_INVALID_HID;
            
        return H5Gcreate2(GetNative(), name.c_str(), 
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    DataFileGroup DataFileGroup::GetSubgroup(const std::string& name) const   
    {
        return H5Gopen2(GetNative(), name.c_str(), H5P_DEFAULT);
    }

    std::vector<std::string> DataFileGroup::EnumerateSubgroups() const
    {
        std::vector<std::string> subgroups;
        H5Literate(GetNative(), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, &EnumGroupsHelper, &subgroups);
        return subgroups;
    }

    DataFileDataset DataFileGroup::CreateDataset(const std::string& name, double value)
    {
        return CreateDatasetRaw(name, {1}, &value);
    }

    DataFileDataset DataFileGroup::CreateDataset(const std::string& name, 
        const Eigen::Matrix<double, -1, -1, Eigen::RowMajor>& value)
    {
        std::size_t rows = static_cast<std::size_t>(value.rows());
        std::size_t cols = static_cast<std::size_t>(value.cols());
        return CreateDatasetRaw(name, {rows, cols}, &(value(0,0)));
    }

    DataFileDataset DataFileGroup::CreateDatasetRaw(const std::string& name, 
        const std::vector<std::size_t>& dims, const double* data)
    {
        if (DoesDatasetExist(name))
            return H5I_INVALID_HID;

        std::vector<hsize_t> hdims(dims.begin(), dims.end());
        hid_t dspace = H5Screate_simple(hdims.size(), hdims.data(), hdims.data());
        if (dspace < 0) return H5I_INVALID_HID;
        auto dspaceGuard = CreateScopeGuard([=](){ H5Sclose(dspace); });
        
        DataFileDataset dataset = H5Dcreate2(GetNative(), name.c_str(), 
            H5T_IEEE_F64LE, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        if (!dataset) 
            return H5I_INVALID_HID;

        if (!dataset.SetRaw(dims, data))
            return H5I_INVALID_HID;

        return dataset;
    }

    bool DataFileGroup::DoesDatasetExist(const std::string& name) const
    {
        if (H5Lexists(GetNative(), name.c_str(), H5P_DEFAULT) <= 0)
            return false;

        bool valid = false;    
        H5E_BEGIN_TRY { valid = GetDataset(name).IsValid(); } H5E_END_TRY;
        return valid;
    }

    DataFileDataset DataFileGroup::GetDataset(const std::string& name) const
    {
        return H5Dopen2(GetNative(), name.c_str(), H5P_DEFAULT);
    }

    std::vector<std::string> DataFileGroup::EnumerateDatasets() const
    {
        std::vector<std::string> subgroups;
        H5Literate(GetNative(), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, &EnumDatasetHelper, &subgroups);
        return subgroups;
    }

    herr_t DataFileGroup::EnumGroupsHelper(hid_t group, const char *name, const H5L_info_t *info, void *data)
    {
        DataFileGroup child = H5Gopen2(group, name, H5P_DEFAULT);
        if (child.IsValid())
            (*static_cast<std::vector<std::string>*>(data)).push_back(name);
        return 0;
    }

    herr_t DataFileGroup::EnumDatasetHelper(hid_t group, const char *name, const H5L_info_t *info, void *data)
    {
        DataFileDataset child = H5Dopen2(group, name, H5P_DEFAULT);
        if (child.IsValid())
            (*static_cast<std::vector<std::string>*>(data)).push_back(name);
        return 0;
    }


    //
    // DataFile
    //

    DataFile::DataFile() 
        : m_file(H5I_INVALID_HID) {}

    DataFile::DataFile(DataFile&& rhs)
        : m_file(rhs.m_file)
    {
        m_file = H5I_INVALID_HID;
    }

    DataFile& DataFile::operator=(DataFile&& rhs)
    {
        std::swap(m_file, rhs.m_file);
        return *this;
    }

    DataFile::~DataFile()
    {
        Close();
    }

    bool DataFile::Open(const std::string& name, DataFileOpenFlag flag)
    {
        if (IsOpen())
            return false;

        if (flag & (DataFile_MUST_NOT_EXIST | DataFile_TRUNCATE))
        {
            // create empty file
            m_file = H5Fcreate(name.c_str(), 
                (flag & DataFile_MUST_NOT_EXIST) ? H5F_ACC_EXCL : H5F_ACC_TRUNC, 
                H5P_DEFAULT, H5P_DEFAULT);
        }
        else if(flag & DataFile_MUST_EXIST)
        {
            // open existing file
            m_file = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (m_file < 0 && !(flag & DataFile_MUST_EXIST))
                m_file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        else
        {
            // default open mode
            // open exisiting file -> if it does not exist create new file
            H5E_BEGIN_TRY {
                m_file = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            } H5E_END_TRY;

            if (m_file < 0 && !(flag & DataFile_MUST_EXIST))
                m_file = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }

        return true;
    }

    void DataFile::Close()
    {
        if (m_file >= 0)
            H5Fclose(m_file);
        m_file = -1;
    }

    bool DataFile::IsOpen() const
    {
        return !(m_file < 0);
    }

    DataFileGroup DataFile::OpenRootGroup()
    {
        if (m_file < 0)
            return H5I_INVALID_HID;

        return H5Gopen2(m_file, "/", H5P_DEFAULT);
    }

}
