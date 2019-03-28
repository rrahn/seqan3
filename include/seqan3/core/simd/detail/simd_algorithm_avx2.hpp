// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithm implementation for AVX2.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <immintrin.h>

#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>

namespace seqan3::detail
{

/*!\brief Implementation of seqan3::simd::unpack_hi for avx2.
 * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[int] first   The vector whose values come before the `second`.
 * \param[int] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <size_t scalar_size, typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 32
//!\endcond
constexpr simd_t unpack_hi(simd_t const & lhs, simd_t const & rhs)
{
    [[maybe_unused]] __m256i tmp_lo{};
    [[maybe_unused]] __m256i tmp_hi{};

    if constexpr (scalar_size == 1)
    {
        tmp_lo = _mm256_unpacklo_epi8(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi8(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 2)
    {
        tmp_lo = _mm256_unpacklo_epi16(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi16(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 4)
    {
        tmp_lo = _mm256_unpacklo_epi32(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi32(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 8)
    {
        tmp_lo = _mm256_unpacklo_epi64(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi64(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }

    return reinterpret_cast<simd_t>(_mm256_permute2f128_si256(tmp_lo, tmp_hi, 0x31));
}

/*!\brief Implementation of seqan3::simd::unpack_lo for avx2.
 * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[int] first   The vector whose values come before the `second`.
 * \param[int] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <size_t scalar_size, typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 32
//!\endcond
constexpr simd_t unpack_lo(simd_t const & lhs, simd_t const & rhs)
{
    [[maybe_unused]] __m256i tmp_lo{};
    [[maybe_unused]] __m256i tmp_hi{};

    if constexpr (scalar_size == 1)
    {
        tmp_lo = _mm256_unpacklo_epi8(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi8(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 2)
    {
        tmp_lo = _mm256_unpacklo_epi16(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi16(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 4)
    {
        tmp_lo = _mm256_unpacklo_epi32(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi32(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else if constexpr (scalar_size == 8)
    {
        tmp_lo = _mm256_unpacklo_epi64(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
        tmp_hi = _mm256_unpackhi_epi64(reinterpret_cast<__m256i const &>(lhs), reinterpret_cast<__m256i const &>(rhs));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }

    return reinterpret_cast<simd_t>(_mm256_permute2f128_si256(tmp_lo, tmp_hi, 0x20));
}

template <Simd simd_t, size_t size>
inline void transpose_matrix_avx2(std::array<simd_t, size> & matrix)
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    constexpr std::array<uint8_t, 32> bitRev{ 0,  8,  4, 12,  2, 10,  6, 14,  1,  9,  5, 13,  3, 11,  7, 15,
                                             16, 24, 20, 28, 18, 26, 22, 30, 17, 25, 21, 29, 19, 27, 23, 31};

    // transpose a 32x32 byte matrix
    __m256i tmp1[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i]    = _mm256_unpacklo_epi8(reinterpret_cast<__m256i &>(matrix[2*i]),
                                          reinterpret_cast<__m256i &>(matrix[2*i+1]));
        tmp1[i+16] = _mm256_unpackhi_epi8(reinterpret_cast<__m256i &>(matrix[2*i]),
                                          reinterpret_cast<__m256i &>(matrix[2*i+1]));
    }
    __m256i  tmp2[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i]    = _mm256_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+16] = _mm256_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i]    = _mm256_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+16] = _mm256_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i]    = _mm256_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+16] = _mm256_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        matrix[bitRev[i]]    = reinterpret_cast<simd_t>(_mm256_permute2x128_si256(tmp2[2*i], tmp2[2*i+1], 0x20));
        matrix[bitRev[i+16]] = reinterpret_cast<simd_t>(_mm256_permute2x128_si256(tmp2[2*i], tmp2[2*i+1], 0x31));
    }
}

} // namespace seqan3::detail
