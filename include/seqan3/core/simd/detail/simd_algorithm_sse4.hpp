// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithm implementation for SSE4.
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

/*!\brief Implementation of seqan3::simd::unpack_hi for sse4.
 * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[int] first   The vector whose values come before the `second`.
 * \param[int] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <size_t scalar_size, Simd simd_t>
//!\cond
    requires simd_traits<std::remove_reference_t<simd_t>>::max_length == 16
//!\endcond
constexpr simd_t unpack_hi(simd_t const & first, simd_t const & second)
{
    if constexpr (scalar_size == 1)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi8(reinterpret_cast<__m128i const &>(first),
                                                          reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi16(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi32(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi64(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

/*!\brief Implementation of seqan3::simd::unpack_lo for sse4.
 * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[int] first   The vector whose values come before the `second`.
 * \param[int] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <size_t scalar_size, typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 16
//!\endcond
constexpr simd_t unpack_lo(simd_t const & first, simd_t const & second)
{
    if constexpr (scalar_size == 1)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi8(reinterpret_cast<__m128i const &>(first),
                                                          reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi16(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi32(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi64(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

template <Simd simd_t, size_t size>
inline void transpose_matrix_sse4(std::array<simd_t, size> & matrix)
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    constexpr std::array<char, 16> bitRev{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};

    // transpose a 16x16 byte matrix
    //
    // matrix =
    // A0 A1 A2 ... Ae Af
    // B0 B1 B2 ... Be Bf
    // ...
    // P0 P1 P2 ... Pe Pf
    __m128i tmp1[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi8(reinterpret_cast<__m128i &>(matrix[2*i]),
                                      reinterpret_cast<__m128i &>(matrix[2*i+1]));
        tmp1[i+8] = _mm_unpackhi_epi8(reinterpret_cast<__m128i &>(matrix[2*i]),
                                      reinterpret_cast<__m128i &>(matrix[2*i+1]));
    }
    // tmp1[0]  = A0 B0 A1 B1 ... A7 B7
    // tmp1[1]  = C0 D0 C1 D1 ... C7 D7
    // ...
    // tmp1[7]  = O0 P0 O1 P1 ... O7 P7
    // tmp1[8]  = A8 B8 A9 B9 ... Af Bf
    // ...
    // tmp1[15] = O8 P8 O9 P9 ... Of Pf
    __m128i tmp2[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp2[i]   = _mm_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+8] = _mm_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    // tmp2[0]  = A0 B0 C0 D0 ... A3 B3 C3 D3
    // tmp2[1]  = E0 F0 G0 H0 ... E3 F3 G3 H3
    // ...
    // tmp2[3]  = M0 N0 O0 P0 ... M3 N3 O3 P3
    // tmp2[4]  = A8 B8 C8 D8 ... Ab Bb Cb Db
    // ...
    // tmp2[7]  = M8 N8 O8 P8 ... Mb Nb Ob Pb
    // tmp2[8]  = A4 B4 C4 D4 ... A7 B7 C7 D7
    // ..
    // tmp2[12] = Ac Bc Cc Dc ... Af Bf Cf Df
    // ...
    // tmp2[15] = Mc Nc Oc Pc ... Mf Nf Of Pf
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+8] = _mm_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    // tmp1[0]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    // ...
    // tmp1[4]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    for (int i = 0; i < 8; ++i)
    {
        matrix[bitRev[i]]   = reinterpret_cast<simd_t>(_mm_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]));
        matrix[bitRev[i+8]] = reinterpret_cast<simd_t>(_mm_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]));
    }
}

} // namespace seqan3::detail
