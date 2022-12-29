// Philipp Neufeld, 2021-2022

#include "UUID.h"

#include <random>
#include <limits>
#include <algorithm>
#include <cstring>

#if defined(QSim_SSSE3)
#include <tmmintrin.h>
#endif

#if defined(QSim_SSE4_1)
#include <smmintrin.h>
#endif

#if defined(QSim_AVX2)
#include <immintrin.h>
#endif

#if defined(__BIG_ENDIAN__)
#error "Please implement UUIDv4 for big endian systems"
#endif

namespace QSim
{

    template<typename size_type>
    struct UUIDv4Hash;
    template<>
    struct UUIDv4Hash<std::uint32_t>
    {
        static std::uint32_t Hash(const void* uuid)
        {
            const std::uint32_t* ptr = reinterpret_cast<const std::uint32_t*>(uuid);
            return ((*ptr ^ *(ptr+1)) ^ (*(ptr+2) ^ *(ptr+3)));
        }
    };
    template<>
    struct UUIDv4Hash<std::uint64_t>
    {
        static std::uint64_t Hash(const void* uuid)
        {
            const std::uint64_t* ptr = reinterpret_cast<const std::uint64_t*>(uuid);
            return (*ptr ^ *(ptr+1));
        }
    };


    UUIDv4::UUIDv4(UUIDv4Register uuid)
    {
        StoreToBufferLEReg(uuid, m_data);
    }

    UUIDv4::UUIDv4() : UUIDv4(UUIDv4::Generate()) {}

    bool UUIDv4::operator==(const UUIDv4& rhs) const
    {
#if defined(QSim_SSE4_1)
        __m128i cmp = _mm_xor_si128(Load(), rhs.Load());
        return !!_mm_testz_si128(cmp, cmp);
#elif defined(QSim_SSE2)
        __m128i cmp = _mm_cmpeq_epi8(Load(), rhs.Load());
        int s = _mm_movemask_epi8(cmp);
        return (s == 0xffff);
#else
        const std::uint64_t* x = reinterpret_cast<const std::uint64_t*>(m_data);
        const std::uint64_t* y = reinterpret_cast<const std::uint64_t*>(rhs.m_data);
        return *x == *y && *(x + 1) == *(y + 1);
#endif
    }

    bool UUIDv4::operator<(const UUIDv4& rhs) const
    {
        const std::uint64_t* x = reinterpret_cast<const std::uint64_t*>(m_data);
        const std::uint64_t* y = reinterpret_cast<const std::uint64_t*>(rhs.m_data);
        return *x < *y || (*x == *y && *(x + 1) < *(y + 1));
    }

    bool UUIDv4::operator!=(const UUIDv4& rhs) const
    {
        return !(*this == rhs);
    }

    bool UUIDv4::operator<=(const UUIDv4& rhs) const
    {
        return !(rhs < *this);
    }

    bool UUIDv4::operator>(const UUIDv4& rhs) const
    {
        return (rhs < *this);
    }

    bool UUIDv4::operator>=(const UUIDv4& rhs) const
    {
        return !(*this < rhs);
    }
 
    std::size_t UUIDv4::Hash() const
    {
        return UUIDv4Hash<std::size_t>::Hash(m_data);
    }

    void UUIDv4::ToString(char* buffer) const
    {
#if defined(QSim_SSE2)
        __m128i uuid = Load();
        __m128i uuids = _mm_srli_epi64(uuid, 4);
        
#if defined(QSim_AVX2)
        const __m256i asciiDigitOff = _mm256_set1_epi8(48);
        const __m256i asciiLetterOffset = _mm256_set1_epi8(55);
#else
        const __m128i asciiDigitOff = _mm_set1_epi8(48);
        const __m128i asciiLetterOffset = _mm_set1_epi8(55);
#endif

        const __m128i mask = _mm_set1_epi8(0x0F);
        __m128i lo = _mm_unpacklo_epi8(uuids, uuid);
        __m128i hi = _mm_unpackhi_epi8(uuids, uuid);
        lo = _mm_and_si128(mask, lo);
        hi = _mm_and_si128(mask, hi);

#if defined(QSim_AVX2)
        const __m256i ten = _mm256_set1_epi8(10);
        __m256i hilo = _mm256_set_m128i(hi, lo);
        __m256i digitMask = _mm256_cmpgt_epi8(ten, hilo);
        __m256i asciiOff = _mm256_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMask);
        hilo = _mm256_add_epi8(hilo, asciiOff);
        lo = _mm256_castsi256_si128(hilo);
        hi = _mm256_extracti128_si256(hilo, 1);
#else
        const __m128i ten = _mm_set1_epi8(10);
        __m128i digitMaskLo = _mm_cmpgt_epi8(ten, lo);
        __m128i digitMaskHi = _mm_cmpgt_epi8(ten, hi);
#if defined(QSim_SSE4_1)
        __m128i asciiOffLo = _mm_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMaskLo);
        __m128i asciiOffHi = _mm_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMaskHi);
#else
        __m128i asciiDigitOffLo = _mm_and_si128(digitMaskLo, asciiDigitOff);
        __m128i asciiLetterOffLo = _mm_andnot_si128(digitMaskLo, asciiLetterOffset);
        __m128i asciiDigitOffHi = _mm_and_si128(digitMaskHi, asciiDigitOff);
        __m128i asciiLetterOffHi = _mm_andnot_si128(digitMaskHi, asciiLetterOffset);
        __m128i asciiOffLo = _mm_or_si128(asciiDigitOffLo, asciiLetterOffLo);
        __m128i asciiOffHi = _mm_or_si128(asciiDigitOffHi, asciiLetterOffHi);
#endif
        lo = _mm_add_epi8(lo, asciiOffLo);
        hi = _mm_add_epi8(hi, asciiOffHi);
#endif

        const __m128i zero = _mm_setzero_si128();
        const __m128i full = _mm_cmpeq_epi32(zero, zero);
        const __m128i mask1 = _mm_set_epi64x(0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFF00000000ull);
        const __m128i mask2 = _mm_xor_si128(full, mask1);
        const __m128i los1 = _mm_srli_si128(lo, 8);
        const __m128i los2 = _mm_srli_si128(lo, 12);
        
        _mm_maskmoveu_si128(hi, mask1, buffer + 20);
        _mm_maskmoveu_si128(hi, mask2, buffer + 19);
        _mm_storeu_si32(buffer + 14, los2);
        _mm_storeu_si32(buffer + 9, los1);
        _mm_storeu_si64(buffer, lo);
#else
        UUIDv4Register uuid = Load();
        auto index = [](auto i) { return i%2==0 ? i+1 : i-1; };
        auto tochar = [](auto i) { return i>9 ? 55 + i : 48 + i; };
        for (int i=0; i<8; i++)
            buffer[i] = tochar((uuid.lo >> (4*index(i))) & 0x0f);
        for (int i=0; i<4; i++)
            buffer[9+i] = tochar((uuid.lo >> (32+4*index(i))) & 0x0f);
        for (int i=0; i<4; i++)
            buffer[14+i] = tochar((uuid.lo >> (48+4*index(i))) & 0x0f);
        for (int i=0; i<4; i++)
            buffer[19+i] = tochar((uuid.hi >> (4*index(i))) & 0x0f);
        for (int i=0; i<12; i++)
            buffer[24+i] = tochar((uuid.hi >> (16+4*index(i))) & 0x0f);
#endif

        buffer[8] = '-';
        buffer[13] = '-';
        buffer[18] = '-';
        buffer[23] = '-';
    }

    void UUIDv4::ToString(std::string& buffer) const
    {
        buffer.resize(36, '0');
        ToString(buffer.data());
    }

    std::string UUIDv4::ToString() const
    {
        std::string str;
        ToString(str);
        return str;
    }

    UUIDv4 UUIDv4::FromString(const char* uuidStr)
    {
#if defined(QSim_SSE2)
#if defined(QSim_AVX2)
        const __m256i asciiDigitOff = _mm256_set1_epi8(48);
        const __m256i asciiLetterOffset = _mm256_set1_epi8(55);
#else
        const __m128i asciiDigitOff = _mm_set1_epi8(48);
        const __m128i asciiLetterOffset = _mm_set1_epi8(55);
#endif

        __m128i lo1 = _mm_loadu_si64(uuidStr);
        __m128i lo2 = _mm_loadu_si32(uuidStr + 9);
        __m128i lo3 = _mm_loadu_si32(uuidStr + 14);
        lo2 = _mm_slli_si128(lo2, 8);
        lo3 = _mm_slli_si128(lo3, 12);
        __m128i lo = _mm_or_si128(lo1, lo2);
        lo = _mm_or_si128(lo, lo3);

        __m128i mask = _mm_cvtsi32_si128(0xFFFFFFFFul);
        __m128i hi1 = _mm_loadu_si32(uuidStr + 19);
        __m128i hi2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(uuidStr + 20));
        hi2 = _mm_andnot_si128(mask, hi2);
        __m128i hi = _mm_or_si128(hi1, hi2);

#if defined(QSim_AVX2)
        const __m256i asciiAm1 = _mm256_set1_epi8(64);
        __m256i hilo = _mm256_set_m128i(hi, lo);
        __m256i digitMask = _mm256_cmpgt_epi8(asciiAm1, hilo);
        __m256i asciiOff = _mm256_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMask);
        hilo = _mm256_sub_epi8(hilo, asciiOff);
        lo = _mm256_castsi256_si128(hilo);
        hi = _mm256_extracti128_si256(hilo, 1);
#else
        const __m128i asciiA = _mm_set1_epi8(65);
        __m128i digitMaskLo = _mm_cmplt_epi8(lo, asciiA);
        __m128i digitMaskHi = _mm_cmplt_epi8(hi, asciiA);
#if defined(QSim_SSE4_1)
        __m128i asciiOffLo = _mm_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMaskLo);
        __m128i asciiOffHi = _mm_blendv_epi8(asciiLetterOffset, asciiDigitOff, digitMaskHi);
#else
        __m128i asciiDigitOffLo = _mm_and_si128(digitMaskLo, asciiDigitOff);
        __m128i asciiLetterOffLo = _mm_andnot_si128(digitMaskLo, asciiLetterOffset);
        __m128i asciiDigitOffHi = _mm_and_si128(digitMaskHi, asciiDigitOff);
        __m128i asciiLetterOffHi = _mm_andnot_si128(digitMaskHi, asciiLetterOffset);
        __m128i asciiOffLo = _mm_or_si128(asciiDigitOffLo, asciiLetterOffLo);
        __m128i asciiOffHi = _mm_or_si128(asciiDigitOffHi, asciiLetterOffHi);
#endif
        lo = _mm_sub_epi8(lo, asciiOffLo);
        hi = _mm_sub_epi8(hi, asciiOffHi);
#endif

        __m128i uuid = PackUUID(hi, lo);

        // flip bits (0x10325476 -> 0x01234567)
        const __m128i flipmask = _mm_set_epi64x(0x0f0f0f0f0f0f0f0full, 0x0f0f0f0f0f0f0f0full);
        __m128i uuid1 = _mm_and_si128(flipmask, uuid);
        __m128i uuid2 = _mm_andnot_si128(flipmask, uuid);
        uuid1 = _mm_slli_epi32(uuid1, 4);
        uuid2 = _mm_srli_epi32(uuid2, 4);
        uuid = _mm_or_si128(uuid1, uuid2);

        return uuid;
#else
        UUIDv4Register uuid;
        uuid.lo = 0;
        uuid.hi = 0;
        auto index = [](auto i) { return i%2==0 ? i+1 : i-1; };
        auto fromchar = [](auto c) { return static_cast<std::uint64_t>(c>64 ? c - 55 : c - 48); };
        for (int i=0; i<8; i++)
            uuid.lo |= (fromchar(uuidStr[i]) & 0x0f) << (4*index(i));
        for (int i=0; i<4; i++)
            uuid.lo |= (fromchar(uuidStr[9+i]) & 0x0f) << (32+4*index(i));
        for (int i=0; i<4; i++)
            uuid.lo |= (fromchar(uuidStr[14+i]) & 0x0f) << (48+4*index(i));
        for (int i=0; i<4; i++)
            uuid.hi |= (fromchar(uuidStr[19+i]) & 0x0f) << (4*index(i));
        for (int i=0; i<12; i++)
            uuid.hi |= (fromchar(uuidStr[24+i]) & 0x0f) << (16+4*index(i));
        return uuid;
#endif
    }

    UUIDv4 UUIDv4::FromString(const std::string& uuidStr)
    {
        return FromString(uuidStr.c_str());
    }

    void UUIDv4::StoreToBufferLE(void* buffer) const
    {
#if defined(QSim_SSE2)
        __m128i uuid = _mm_loadu_si128(reinterpret_cast<const __m128i*>(m_data));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(buffer), uuid);
#else
        std::copy_n(m_data, 16, static_cast<std::uint8_t*>(buffer));
#endif
    }
    
    void UUIDv4::StoreToBufferBE(void* buffer) const
    {
#if defined(QSim_SSE2)
        __m128i uuid = _mm_loadu_si128(reinterpret_cast<const __m128i*>(m_data));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(buffer), SwapEndianess(uuid));
#else
        for (int i=0; i<16; i++) 
            static_cast<std::uint8_t*>(buffer)[i] = m_data[15-i];
#endif
    }

    UUIDv4 UUIDv4::LoadFromBufferLE(const void* buffer)
    {
        return LoadFromBufferLEReg(buffer);
    }
    
    UUIDv4 UUIDv4::LoadFromBufferBE(const void* buffer)
    {
        return LoadFromBufferBEReg(buffer);
    }

    UUIDv4Register UUIDv4::LoadFromBufferLEReg(const void* buffer)
    {
#if defined(QSim_SSE2)
        return _mm_loadu_si128(reinterpret_cast<const __m128i*>(buffer));
#else
        return *reinterpret_cast<const UUIDv4Register*>(buffer);
#endif
    }
    
    UUIDv4Register UUIDv4::LoadFromBufferBEReg(const void* buffer)
    {
        return SwapEndianess(LoadFromBufferLEReg(buffer));
    }

    void UUIDv4::StoreToBufferLEReg(UUIDv4Register uuid, void* buffer)
    {
#if defined(QSim_SSE2)
        _mm_storeu_si128(reinterpret_cast<__m128i*>(buffer), uuid);
#else
        std::copy_n(reinterpret_cast<const std::uint8_t*>(&uuid), 16, 
            reinterpret_cast<std::uint8_t*>(buffer));
#endif
    }

    void UUIDv4::StoreToBufferBEReg(UUIDv4Register uuid, void* buffer)
    {
        StoreToBufferLEReg(SwapEndianess(uuid), buffer);
    } 

    UUIDv4Register UUIDv4::Load() const
    {
        return LoadFromBufferLEReg(m_data);
    }
    
    void UUIDv4::Store(UUIDv4Register uuid)
    {
        StoreToBufferLEReg(uuid, m_data);
    }
    
    UUIDv4Register UUIDv4::Generate()
    {
        // Generates version 4 variant 1 UUIDs
        static std::random_device rdev;
        static std::mt19937_64 engine(rdev());
        static std::uniform_int_distribution<std::uint64_t> dist(
            std::numeric_limits<std::uint64_t>::min(), 
            std::numeric_limits<std::uint64_t>::max());

#if defined(QSim_SSE2)
        const __m128i and_mask = _mm_set_epi64x(0xFFFFFFFFFFFFFF3Full, 0xFF0FFFFFFFFFFFFFull);
        const __m128i or_mask =  _mm_set_epi64x(0x0000000000000080ull, 0x0040000000000000ull); // v4.1
        __m128i v = _mm_set_epi64x(dist(engine), dist(engine));
        __m128i uuid = _mm_or_si128(_mm_and_si128(v, and_mask), or_mask);
#else
        UUIDv4Register uuid;
        uuid.hi = (dist(engine) & 0xFFFFFFFFFFFFFF3Full) | 0x0000000000000080ull;
        uuid.lo = (dist(engine) & 0xFF0FFFFFFFFFFFFFull) | 0x0040000000000000ull;
#endif

        return uuid;
    }

    UUIDv4Register UUIDv4::SwapEndianess(UUIDv4Register v)
    {
#if defined(QSim_SSSE3)
        return _mm_shuffle_epi8(v, _mm_set_epi64x(0x0001020304050607ull, 0x08090A0B0C0D0E0Full));
#elif defined(QSim_SSE2)
        v = _mm_shuffle_epi32(v, _MM_SHUFFLE(0,1,2,3));
        v = _mm_shufflelo_epi16(v, _MM_SHUFFLE(2,3,0,1));
        v = _mm_shufflehi_epi16(v, _MM_SHUFFLE(2,3,0,1));

        const __m128i mask = _mm_set_epi64x(0xFF00FF00FF00FF00ull, 0xFF00FF00FF00FF00ull);
        __m128i tmp1 = _mm_slli_epi16(v, 8);
        __m128i tmp2 = _mm_srli_epi16(v, 8);
        tmp1 = _mm_and_si128(mask, tmp1);
        tmp2 = _mm_andnot_si128(mask, tmp2);

        return _mm_or_si128(tmp1, tmp2);
#else
    UUIDv4Register res;
    for (int i=0; i<16; i++) 
            reinterpret_cast<std::uint8_t*>(&res)[i] = reinterpret_cast<std::uint8_t*>(&v)[15-i];
    return res;
#endif
    }

#if defined(QSim_SSE2)
    __m128i UUIDv4::PackUUID(__m128i hi, __m128i lo)
    {
#if defined(QSim_SSSE3)
        const __m128 zerof = _mm_setzero_ps();
        const __m128i shufmask = _mm_set_epi64x(0x0f0d0b0907050301ull, 0x0e0c0a0806040200);

        __m128i hishuf = _mm_shuffle_epi8(hi, shufmask);
        __m128 hishuff = _mm_castsi128_ps(hishuf);
        __m128i hi1 = _mm_castps_si128(_mm_movehl_ps(hishuff, zerof));
        __m128i hi2 = _mm_castps_si128(_mm_movelh_ps(zerof, hishuff));
        hi1 = _mm_slli_epi32(hi1, 4);
        hi = _mm_or_si128(hi1, hi2);

        __m128i loshuf = _mm_shuffle_epi8(lo, shufmask);
        __m128 loshuff = _mm_castsi128_ps(loshuf);
        __m128i lo1 = _mm_castps_si128(_mm_movehl_ps(zerof, loshuff));
        __m128i lo2 = _mm_castps_si128(_mm_movelh_ps(loshuff, zerof));
        lo1 = _mm_slli_epi32(lo1, 4);
        lo = _mm_or_si128(lo1, lo2);
        
        return _mm_or_si128(lo, hi);
#else
        const __m128i zero = _mm_setzero_si128();
        __m128i mask = _mm_set_epi64x(0xFFFF0000FFFF0000ull, 0xFFFF0000FFFF0000ull);
        
        __m128i lo1 =  _mm_unpacklo_epi8(lo, _mm_setzero_si128());
        __m128i lo11 = _mm_andnot_si128(mask, lo1);
        __m128i lo12 = _mm_srli_epi32(_mm_and_si128(mask, lo1), 12);
        lo11 = _mm_packus_epi16(lo11, lo12);
        lo12 = _mm_shuffle_epi32(lo11, _MM_SHUFFLE(1, 0, 3, 2));
        lo1 = _mm_or_si128(lo11, lo12);

        __m128i hi1 =  _mm_unpackhi_epi8(lo, _mm_setzero_si128());
        __m128i hi11 = _mm_andnot_si128(mask, hi1);
        __m128i hi12 = _mm_srli_epi32(_mm_and_si128(mask, hi1), 12);
        hi11 = _mm_packus_epi16(hi11, hi12);
        hi12 = _mm_shuffle_epi32(hi11, _MM_SHUFFLE(1, 0, 3, 2));
        hi1 = _mm_or_si128(hi11, hi12);

        __m128i lo2 =  _mm_unpacklo_epi8(hi, _mm_setzero_si128());
        __m128i lo21 = _mm_andnot_si128(mask, lo2);
        __m128i lo22 = _mm_srli_epi32(_mm_and_si128(mask, lo2), 12);
        lo21 = _mm_packus_epi16(lo21, lo22);
        lo22 = _mm_shuffle_epi32(lo21, _MM_SHUFFLE(1, 0, 3, 2));
        lo2 = _mm_or_si128(lo21, lo22);

        __m128i hi2 =  _mm_unpackhi_epi8(hi, _mm_setzero_si128());
        __m128i hi21 = _mm_andnot_si128(mask, hi2);
        __m128i hi22 = _mm_srli_epi32(_mm_and_si128(mask, hi2), 12);
        hi21 = _mm_packus_epi16(hi21, hi22);
        hi22 = _mm_shuffle_epi32(hi21, _MM_SHUFFLE(1, 0, 3, 2));
        hi2 = _mm_or_si128(hi21, hi22);

        __m128i res1 = _mm_unpacklo_epi64(lo1, hi1);
        __m128i res2 = _mm_unpacklo_epi64(lo2, hi2);
        return _mm_packus_epi16(res1, res2);
#endif
    }
#endif



}
