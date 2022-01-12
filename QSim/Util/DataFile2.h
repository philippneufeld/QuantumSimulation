// Philipp Neufeld, 2021-2022

#ifndef QSim_Util_DataFile2_H_
#define QSim_Util_DataFile2_H_

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>


namespace QSim
{

    using DataFileOffset = std::uint32_t;

    // 
    // RAII Wrapper class for OS file handle
    //
    class DataFileOSHandle
    {
    public:
        DataFileOSHandle(std::FILE* pFile);
        ~DataFileOSHandle();

        DataFileOSHandle(const DataFileOSHandle&) = delete;
        DataFileOSHandle(DataFileOSHandle&& rhs);

        DataFileOSHandle& operator=(const DataFileOSHandle&) = delete;
        DataFileOSHandle& operator=(DataFileOSHandle&& rhs);

        operator std::FILE*() { return m_pFile; }
        operator const std::FILE*() const { return m_pFile; }

        // raw read/write
        std::size_t Write(const void* ptr, std::size_t n) const;
        std::size_t Read(void* ptr, std::size_t n) const;

        // file cursor positioning
        bool SetCursorPos(DataFileOffset idx) const;
        DataFileOffset GetCursorPos() const;

        bool IsEOF() const;

    private:
        std::FILE* m_pFile;
    };

    
    namespace Internal
    {
        // unique id struct
        struct DataFileUuid
        {
            std::uint32_t m_data1;
            std::uint16_t m_data2;
            std::uint16_t m_data3;
            std::uint8_t m_data4[8];
        };

        struct DataFileHeader
        {
            DataFileUuid uuid;
            std::uint32_t verMaj;
            std::uint32_t verMin;
            std::uint32_t blockSize;
            DataFileOffset handleTblAddr;
        };
    }


    class DataFileDriver
    {
        static const Internal::DataFileUuid Uuid_v;

        static constexpr std::uint32_t MajorVersion_v = 1;
        static constexpr std::uint32_t MinorVersion_v = 0;
        static constexpr std::uint32_t MajVerSupport_v = 1;
        static constexpr std::uint32_t MinVerSupport_v = 0;

        static constexpr std::uint32_t DefaultBlockSize_v = 1024;

    public:
        DataFileDriver();

        bool OpenFile(const std::string& path);

        std::weak_ptr<DataFileOSHandle> GetOSHandle() const;

        bool WriteFileHeader();
        bool ReadFileHeader();

        bool WriteHandleTable(DataFileOffset addr);
        bool ReadHandleTable(DataFileOffset addr);

        std::size_t WriteData(const void* data, std::size_t size);
        bool ReadData(DataFileOffset addr, void* data, std::size_t size);

    private:
        std::shared_ptr<DataFileOSHandle> m_pFile;

        Internal::DataFileHeader m_header;
    };

}

#endif
