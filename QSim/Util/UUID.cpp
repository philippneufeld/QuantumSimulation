// Philipp Neufeld, 2021-2022

#include "UUID.h"

#include <random>
#include <limits>
#include <algorithm>

#include <iostream>

// SSE and SSE2 are always available
#include <xmmintrin.h> // SSE
#include <emmintrin.h> // SSE2

namespace QSim
{

    UUIDv4::UUIDv4() 
    {
        UUIDv4::Generate(m_data);
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
        const __m128i mask1 = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFF00000000);
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

    void UUIDv4::Generate(std::uint8_t* buffer)
    {
        // Generates version 4 variant 1 UUIDs
        static std::random_device rdev;
        static std::mt19937_64 engine(/*rdev()*/0);
        static std::uniform_int_distribution<std::uint64_t> dist(
            std::numeric_limits<std::uint64_t>::min(), 
            std::numeric_limits<std::uint64_t>::max());

        const __m128i and_mask = _mm_set_epi64x(0xFFFFFFFFFFFFFF3Full, 0xFF0FFFFFFFFFFFFFull);
        const __m128i or_mask =  _mm_set_epi64x(0x0000000000000080ull, 0x0040000000000000ull); // v4.1
        __m128i v = _mm_set_epi64x(dist(engine), dist(engine));
        __m128i uuid = _mm_or_si128(_mm_and_si128(v, and_mask), or_mask);

        // TODO: Big endian systems handling

        _mm_storeu_si128(reinterpret_cast<__m128i*>(buffer), uuid);
    }

}
