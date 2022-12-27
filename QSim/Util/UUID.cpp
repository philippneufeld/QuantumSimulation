// Philipp Neufeld, 2021-2022

#include "UUID.h"

#include <random>
#include <limits>
#include <algorithm>

#include <iostream>

// SSE and SSE2 are always available
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2

#include <tmmintrin.h> // SSE2


namespace QSim
{

    UUIDv4::UUIDv4(__m128i uuid)
    {
        _mm_storeu_si128(reinterpret_cast<__m128i*>(m_data), uuid);
    }

    UUIDv4::UUIDv4() : UUIDv4(UUIDv4::Generate()) {}

    bool UUIDv4::operator==(const UUIDv4& rhs) const
    {
        const __m128i uuid1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(m_data));
        const __m128i uuid2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(rhs.m_data));
        __m128i cmp = _mm_cmpeq_epi8(uuid1, uuid2);
        int s = _mm_movemask_epi8(cmp);
        return (s == 0xffff);
    }

    void UUIDv4::ToString(char* buffer) const
    {
        __m128i uuid = _mm_loadu_si128(reinterpret_cast<const __m128i*>(m_data));
        __m128i uuids = _mm_srli_epi64(uuid, 4);
        
        const __m128i asciiDigitOff = _mm_set1_epi8(48);
        const __m128i asciiLetterOffset = _mm_set1_epi8(55);

        const __m128i mask = _mm_set1_epi8(0x0F);
        __m128i lo = _mm_unpacklo_epi8(uuids, uuid);
        __m128i hi = _mm_unpackhi_epi8(uuids, uuid);
        lo = _mm_and_si128(mask, lo);
        hi = _mm_and_si128(mask, hi);

        const __m128i ten = _mm_set1_epi8(10);
        __m128i digitMaskLo = _mm_cmplt_epi8(lo, ten);
        __m128i digitMaskHi = _mm_cmplt_epi8(hi, ten);
        __m128i asciiDigitOffLo = _mm_and_si128(digitMaskLo, asciiDigitOff);
        __m128i asciiLetterOffLo = _mm_andnot_si128(digitMaskLo, asciiLetterOffset);
        __m128i asciiDigitOffHi = _mm_and_si128(digitMaskHi, asciiDigitOff);
        __m128i asciiLetterOffHi = _mm_andnot_si128(digitMaskHi, asciiLetterOffset);
        __m128i asciiOffLo = _mm_or_si128(asciiDigitOffLo, asciiLetterOffLo);
        __m128i asciiOffHi = _mm_or_si128(asciiDigitOffHi, asciiLetterOffHi);

        lo = _mm_add_epi8(lo, asciiOffLo);
        hi = _mm_add_epi8(hi, asciiOffHi);

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
        const __m128i asciiDigitOff = _mm_set1_epi8(48);
        const __m128i asciiLetterOffset = _mm_set1_epi8(55);

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

        const __m128i asciiA = _mm_set1_epi8(65);
        __m128i digitMaskLo = _mm_cmplt_epi8(lo, asciiA);
        __m128i digitMaskHi = _mm_cmplt_epi8(hi, asciiA);
        __m128i asciiDigitOffLo = _mm_and_si128(digitMaskLo, asciiDigitOff);
        __m128i asciiLetterOffLo = _mm_andnot_si128(digitMaskLo, asciiLetterOffset);
        __m128i asciiDigitOffHi = _mm_and_si128(digitMaskHi, asciiDigitOff);
        __m128i asciiLetterOffHi = _mm_andnot_si128(digitMaskHi, asciiLetterOffset);
        __m128i asciiOffLo = _mm_or_si128(asciiDigitOffLo, asciiLetterOffLo);
        __m128i asciiOffHi = _mm_or_si128(asciiDigitOffHi, asciiLetterOffHi);
        
        lo = _mm_sub_epi8(lo, asciiOffLo);
        hi = _mm_sub_epi8(hi, asciiOffHi);

        const __m128 zerof = _mm_setzero_ps();
        const __m128i shufmask = _mm_set_epi64x(0x0f0d0b0907050301ull, 0x0e0c0a0806040200);
        __m128i hishuf = _mm_shuffle_epi8(hi, shufmask);
        __m128 hishuff = _mm_castsi128_ps(hishuf);
        hi1 = _mm_castps_si128(_mm_movehl_ps(hishuff, zerof));
        hi2 = _mm_castps_si128(_mm_movelh_ps(zerof, hishuff));
        hi1 = _mm_slli_epi32(hi1, 4);
        hi = _mm_or_si128(hi1, hi2);

        __m128i loshuf = _mm_shuffle_epi8(lo, shufmask);
        __m128 loshuff = _mm_castsi128_ps(loshuf);
        lo1 = _mm_castps_si128(_mm_movehl_ps(zerof, loshuff));
        lo2 = _mm_castps_si128(_mm_movelh_ps(loshuff, zerof));
        lo1 = _mm_slli_epi32(lo1, 4);
        lo = _mm_or_si128(lo1, lo2);

        __m128i uuid = _mm_or_si128(lo, hi);

        const __m128i flipmask = _mm_set_epi64x(0x0f0f0f0f0f0f0f0full, 0x0f0f0f0f0f0f0f0full);
        __m128i uuid1 = _mm_and_si128(flipmask, uuid);
        __m128i uuid2 = _mm_andnot_si128(flipmask, uuid);
        uuid1 = _mm_slli_epi32(uuid1, 4);
        uuid2 = _mm_srli_epi32(uuid2, 4);
        uuid = _mm_or_si128(uuid1, uuid2);

        return uuid;
    }

    UUIDv4 UUIDv4::FromString(const std::string& uuidStr)
    {
        return FromString(uuidStr.c_str());
    }

    __m128i UUIDv4::Generate()
    {
        // Generates version 4 variant 1 UUIDs
        static std::random_device rdev;
        static std::mt19937_64 engine(rdev());
        static std::uniform_int_distribution<std::uint64_t> dist(
            std::numeric_limits<std::uint64_t>::min(), 
            std::numeric_limits<std::uint64_t>::max());

        const __m128i and_mask = _mm_set_epi64x(0xFFFFFFFFFFFFFF3Full, 0xFF0FFFFFFFFFFFFFull);
        const __m128i or_mask =  _mm_set_epi64x(0x0000000000000080ull, 0x0040000000000000ull); // v4.1
        __m128i v = _mm_set_epi64x(dist(engine), dist(engine));
        __m128i uuid = _mm_or_si128(_mm_and_si128(v, and_mask), or_mask);

        // TODO: Big endian systems handling

        return uuid;
    }

}
