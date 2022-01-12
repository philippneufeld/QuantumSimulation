// Philipp Neufeld, 2021-2022

#include "DataFile2.h"

#include <algorithm>

namespace QSim
{
    DataFileOSHandle::DataFileOSHandle(std::FILE* pFile)
        : m_pFile(pFile) { }

    DataFileOSHandle::~DataFileOSHandle()
    {
        if(m_pFile)
            std::fclose(m_pFile);
    }

    DataFileOSHandle::DataFileOSHandle(DataFileOSHandle&& rhs)
        : m_pFile(rhs.m_pFile)
    {
        rhs.m_pFile = nullptr;
    }

    DataFileOSHandle& DataFileOSHandle::operator=(DataFileOSHandle&& rhs)
    {
        std::swap(m_pFile, rhs.m_pFile);
        return *this;
    }

    std::size_t DataFileOSHandle::Write(const void* ptr, std::size_t n) const
    {
        return std::fwrite(ptr, 1, n, m_pFile);
    }

    std::size_t DataFileOSHandle::Read(void* ptr, std::size_t n) const
    {
        return std::fread(ptr, 1, n, m_pFile);
    }

    bool DataFileOSHandle::SetCursorPos(DataFileOffset idx) const
    {
        long off = static_cast<long>(idx);
        return (std::fseek(m_pFile, off, SEEK_SET) == 0);
    }

    DataFileOffset DataFileOSHandle::GetCursorPos() const
    {
        return static_cast<DataFileOffset>(std::ftell(m_pFile));
    }

    bool DataFileOSHandle::IsEOF() const
    {
        return (std::feof(m_pFile) != 0);
    }


    // uuid that was generated using the Linux command uuidgen
    const Internal::DataFileUuid DataFileDriver::Uuid_v = {
            0x3CFA8C2D, 0x63D7, 0x4871, 0x89, 0x70, 
            0x90, 0xD5, 0x12, 0x53, 0x30, 0x30};

    DataFileDriver::DataFileDriver()
        : m_pFile(nullptr), m_header{}
    {
    }

    bool DataFileDriver::OpenFile(const std::string& path)
    {
        // release file
        m_pFile.reset();

        // open file for update
        FILE* pFile = std::fopen(path.c_str(), "r+b");
        if (pFile)
        {
            m_pFile = std::make_shared<DataFileOSHandle>(pFile);
            if (!ReadFileHeader())
            {
                m_pFile.reset();
                return false;
            }
        }
        else
        {
            // create new file
            pFile = std::fopen(path.c_str(), "w+b");
            if (!pFile)
                return false;
            m_pFile = std::make_shared<DataFileOSHandle>(pFile);

            // write file header
            m_header.verMaj = MajorVersion_v;
            m_header.verMin = MinorVersion_v;
            m_header.blockSize = DefaultBlockSize_v;
            m_header.handleTblAddr = sizeof(m_header);
            if (!WriteFileHeader())
            {
                m_pFile.reset();
                return false;
            }
        }
        
        return !!m_pFile;
    }

    std::weak_ptr<DataFileOSHandle> DataFileDriver::GetOSHandle() const
    {
        return m_pFile;
    }

    bool DataFileDriver::WriteFileHeader()
    {
        // set right uuid
        auto dstPtr = reinterpret_cast<std::uint8_t*>(&m_header.uuid);
        auto srcPtr = reinterpret_cast<const std::uint8_t*>(&Uuid_v);
        std::copy_n(srcPtr, sizeof(m_header.uuid), dstPtr);
        
        if (!m_pFile)
            return false;

        if (!m_pFile->SetCursorPos(0))
            return false;

        if (m_pFile->Write(&m_header, sizeof(m_header)) != sizeof(m_header))
            return false;

        return true;
    }

    bool DataFileDriver::ReadFileHeader()
    {
        // reset header data
        m_header = {};
        if (!m_pFile)
            return false;

        if (!m_pFile->SetCursorPos(0))
            return false;

        if (m_pFile->Read(&m_header, sizeof(m_header)) != sizeof(m_header))
            return false;

        // check if the uuid is right
        auto ptr1 = reinterpret_cast<std::uint8_t*>(&m_header.uuid);
        auto ptr2 = reinterpret_cast<const std::uint8_t*>(&Uuid_v);
        if (!std::equal(ptr1, ptr1 + sizeof(m_header.uuid), ptr2))
            return false;

        // file too old -> unsupported
        if ((m_header.verMaj < MajVerSupport_v) ||
            (m_header.verMaj == MajVerSupport_v && m_header.verMin < MinVerSupport_v))
            return false;

        return true;
    }

}
