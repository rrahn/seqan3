// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains algorithms to modify seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <immintrin.h>
#include <utility>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_sse4.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_avx2.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_avx512.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd::fill.
//!\ingroup simd
template <Simd simd_t, size_t... I>
constexpr simd_t fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>)
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd::iota.
//!\ingroup simd
template <Simd simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{static_cast<scalar_t>(offset + I)...};
}

/*!\brief Helper function to extract a part of the given simd vector.
 * \ingroup simd
 * \tparam divisor The divisor to select the chunk size.
 * \tparam simd_t  The simd type; must model seqan3::Simd.
 *
 * \param[in] src  The source vector to extract from.
 * \param[in] mask The control mask to select which chunk is extracted.
 * \returns The destination vector containing the extracted part.
 *
 * \details
 *
 * Extracts the specified part of the source simd vector and stores it in the first chunk starting at offset 0 in
 * the destination vector.
 */
template <size_t divisor, Simd simd_t>
constexpr simd_t extract_impl(simd_t const & src, uint8_t const mask)
{
    simd_t dst{};
    constexpr size_t chunk = simd_traits<simd_t>::length / divisor;
    size_t offset = chunk * mask;
    for (size_t i = offset; i < offset + chunk; ++i)
        dst[i - offset] = src[i];

    return dst;
}

/*!\brief Transforms a batch of ranges into simd types (Structure-of-Arrays) and writes them into the output iterator.
 * \tparam simd_t        The simd type; must model seqan3::simd::simd_concept.
 * \tparam output_iter_t The output iterator type; must model std::OutputIterator with `simd_t`.
 * \tparam range_t       The input range type; must model std::ranges::ForwardRange and its reference type must model
 *                       std::ranges::ViewableRange and std::ranges::ForwardRange over seqan3::Alphabet.
 * \param out     The output iterator.
 * \param seq_rng The range over the ranges to transform.
 * \ingroup simd
 *
 * \details
 *
 * This function assumes that the batch contains simd_traits<simd_t>::length many sequences which have the same
 * size. At most simd_traits<simd_t>::length characters of the sequence are transformed or less if the ranges are
 * smaller. If less sequences are provided or the sequences have different size the behaviour is undefined.
 */
template <Simd simd_t,
          std::OutputIterator<simd_t> output_iter_t,
          std::ranges::InputRange range_t>
//!\cond
    requires std::ranges::ForwardRange<reference_t<range_t>> &&
             std::ranges::ViewableRange<reference_t<range_t>> &&
             seqan3::Alphabet<reference_t<reference_t<range_t>>>
//!\endcond
constexpr void transform_batch_to_soa(output_iter_t out, range_t && seq_rng)
{
    // Check if the length is identical in debug mode.
    if constexpr (std::ranges::SizedRange<range_t>)
        assert(std::ranges::size(seq_rng) == static_cast<size_t>(simd_traits<simd_t>::length));

    // stack memory for the cached iterator - sentinel pairs.
    using it_pair_t = std::pair<std::ranges::iterator_t<reference_t<range_t>>,
                                std::ranges::sentinel_t<reference_t<range_t>>>;
    std::array<it_pair_t, simd_traits<simd_t>::length> iter_cache;

    // Cache the iterator, sentinel pairs of the underlying ranges.
    for (auto && [rng, idx] : std::view::zip(seq_rng, std::view::iota(0)))
        iter_cache[idx] = {std::ranges::begin(rng), std::ranges::end(rng)};

    // Iterate over the length of the array or the maximal size of the underlying range.
    // We assume all ranges have the same length.
    simd_t simd{};
    for (size_t j = 0; (j < iter_cache.size()) || (iter_cache[j].first == iter_cache[j].second); ++j, ++out)
    {
        // Fill the simd value with the ranks of the respective alphabets.
        for (size_t i = 0; i < iter_cache.size(); ++i)
        {
            simd[i] = to_rank(*iter_cache[i].first);
            ++(iter_cache[i].first);
        }
        // Store the simd vector.
        *out = simd;
    }
}

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief Fills a seqan3::simd::simd_type vector with a scalar value.
 * \tparam    simd_t The simd type which satisfies seqan3::simd::Simd.
 * \param[in] scalar The scalar value to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/fill.cpp
 */
template <Simd simd_t>
constexpr simd_t fill(typename simd_traits<simd_t>::scalar_type const scalar)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd::simd_type vector with the scalar values offset, offset+1, offset+2, ...
 * \tparam    simd_t The simd type which satisfies seqan3::simd::Simd.
 * \param[in] offset The scalar offset to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/iota.cpp
 */
template <Simd simd_t>
constexpr simd_t iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

/*!\brief Load simd_t size bits of integral data from memory.
 * \tparam    simd_t    The simd type; must model seqan3::simd::Simd.
 * \param[in] mem_addr  The memory address to load from. Does not need to be aligned on any particular boundary.
 * \ingroup simd
 */
template <Simd simd_t>
constexpr simd_t load(void const * mem_addr)
{
    assert(mem_addr != nullptr);

    if constexpr (simd_traits<simd_t>::max_length == 16)
        return reinterpret_cast<simd_t>(_mm_loadu_si128(reinterpret_cast<__m128i const *>(mem_addr)));
    else if constexpr (simd_traits<simd_t>::max_length == 32)
        return reinterpret_cast<simd_t>(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(mem_addr)));
    else if constexpr (simd_traits<simd_t>::max_length == 64)
        return reinterpret_cast<simd_t>(_mm512_loadu_si512(mem_addr));
    else
        static_assert(simd_traits<simd_t>::max_length <= 64, "Unsupported simd type.");
}

/*!\brief Sets the given values to the seqan3::simd::simd_type vector.
 * \tparam    simd_t    The simd type; must model seqan3::simd::Simd.
 * \param[in] scalar... The scalar values to set.
 * \ingroup simd
 */
template <Simd simd_t, typename ... scalar_t>
constexpr simd_t set(scalar_t const & ... scalar)
{
    static_assert((std::ConvertibleTo<scalar_t, typename simd_traits<simd_t>::scalar_type> && ...),
                  "The provided types must be convertible to the scalar type of the simd vector.");
    static_assert(sizeof...(scalar) == simd_traits<simd_t>::length,
                  "The number of set values must be equal to the length of the simd vector.");

    return simd_t{scalar...};
}

/*!\brief Transposes a quadratic matrix with simd vectors.
 * \tparam        simd_t The simd type; must model seqan3::simd::Simd and size of scalar type must be one byte.
 * \tparam        size   The size of the matrix; must be square of `simd_traits<simd_t>::length`.
 * \param[in,out] matrix The matrix to transpose.
 * \ingroup simd
 */
template <Simd simd_t, size_t size>
//!\cond
    requires (detail::is_builtin_simd_v<simd_t>) && (simd_traits<simd_t>::length == simd_traits<simd_t>::max_length)
//!\endcond
constexpr void transpose(std::array<simd_t, size> & matrix)
{
    static_assert(size == simd_traits<simd_t>::length, "Only a quadratic byte matrix can be transposed.");

    if constexpr (simd_traits<simd_t>::length == 16)
        detail::transpose_matrix_sse4(matrix);
    else if constexpr (simd_traits<simd_t>::length == 32)
        detail::transpose_matrix_avx2(matrix);
    else
        static_assert(simd_traits<simd_t>::length <= 32,
                      "Transpose is not supported for this SIMD extension.");
}

/*!\brief Upcasts the given vector into the target vector.
 * \tparam        simd_t The simd type; must model seqan3::simd::Simd and must be a builtin simd type.
 * \param[in,out] matrix The matrix to transpose.
 * \ingroup simd
 */
template <Simd target_simd_t, Simd source_simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<target_simd_t>) && (detail::is_builtin_simd_v<source_simd_t>)
//!\endcond
constexpr target_simd_t upcast(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
                  "The length of the target simd type must be greater or equal than the length of the source simd type.");

    if constexpr (simd_traits<source_simd_t>::length == simd_traits<target_simd_t>::length)
    {
        static_assert(simd_traits<target_simd_t>::max_length == simd_traits<source_simd_t>::max_length,
                      "Target vector has different byte size.");
        return reinterpret_cast<target_simd_t>(src);  // Same packing so we do not cast.
    }
    if constexpr (simd_traits<source_simd_t>::max_length == 16) // SSE4
    {
        static_assert(simd_traits<target_simd_t>::max_length == 16, "Target vector has different byte size.");

        if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi8 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 8) // to epi16
                return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi16(reinterpret_cast<__m128i const &>(src)));
            if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
                return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi32(reinterpret_cast<__m128i const &>(src)));
            if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
                return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi64(reinterpret_cast<__m128i const &>(src)));
        }
        else if constexpr (simd_traits<source_simd_t>::length == 8) // cast from epi16 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
                return reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(reinterpret_cast<__m128i const &>(src)));
            if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
                return reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi64(reinterpret_cast<__m128i const &>(src)));
        }
        else  // cast from epi32 to epi64
        {
            static_assert(simd_traits<source_simd_t>::length == 4, "Expected 32 bit scalar type.");
            return reinterpret_cast<target_simd_t>(_mm_cvtepi32_epi64(reinterpret_cast<__m128i const &>(src)));
        }
    }
    else if constexpr (simd_traits<source_simd_t>::max_length == 32) // AVX2
    {
        static_assert(simd_traits<target_simd_t>::max_length == 32, "Target vector has different byte size.");
        __m128i const & tmp = _mm256_castsi256_si128(reinterpret_cast<__m256i const &>(src));
        if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi8 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 16) // to epi16
                return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi16(tmp));
            if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
                return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi32(tmp));
            if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
                return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi64(tmp));
        }
        else if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi16 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
                return reinterpret_cast<target_simd_t>(_mm256_cvtepi16_epi32(tmp));
            if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
                return reinterpret_cast<target_simd_t>(_mm256_cvtepi16_epi64(tmp));
        }
        else // cast from epi32 to epi64
        {
            static_assert(simd_traits<source_simd_t>::length == 8, "Expected 32 bit scalar type.");
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi32_epi64(tmp));
        }
    }
    else if constexpr (simd_traits<source_simd_t>::max_length == 64) // AVX512
    {
        static_assert(simd_traits<target_simd_t>::max_length == 64, "Target vector has different byte size.");
        __m512i const & tmp = reinterpret_cast<__m512i const &>(src);
        if constexpr (simd_traits<source_simd_t>::length == 64) // cast from epi8 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 32) // to epi16
                return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi16(_mm512_castsi512_si256(tmp)));
            if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
                return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi32(_mm512_castsi512_si128(tmp)));
            if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
                return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi64(_mm512_castsi512_si128(tmp)));
        }
        else if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi16 ...
        {
            if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
                return reinterpret_cast<target_simd_t>(_mm512_cvtepi16_epi32(_mm512_castsi512_si256(tmp)));
            if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
                return reinterpret_cast<target_simd_t>(_mm512_cvtepi16_epi64(_mm512_castsi512_si128(tmp)));
        }
        else // cast from epi32 to epi64
        {
            static_assert(simd_traits<source_simd_t>::length == 16, "Expected 32 bit scalar type.");
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi32_epi64(_mm512_castsi512_si256(tmp)));
        }
    }
    else
    {
        static_assert(simd_traits<source_simd_t>::max_length <= 32, "Simd type is not supported.");
    }
}

template <Simd target_simd_t, Simd source_simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<target_simd_t>) && (detail::is_builtin_simd_v<source_simd_t>)
//!\endcond
constexpr target_simd_t upcast_auto(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
                  "The length of the target type must be less than or equal to the length of the source type.");
    target_simd_t target{};

    for (size_t i = 0; i < simd_traits<target_simd_t>::length; ++i)
        target[i] = src[i];

    return target;
}

template <Simd simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<simd_t>)
//!\endcond
constexpr simd_t extract_halve(simd_t const & src, uint8_t const mask)
{
    return detail::extract_impl<2>(src, mask & 0x01);
}

template <Simd simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<simd_t>)
//!\endcond
constexpr simd_t extract_quarter(simd_t const & src, uint8_t const mask)
{
    return detail::extract_impl<4>(src, mask & 0x03);
}

template <Simd simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<simd_t>)
//!\endcond
constexpr simd_t extract_eighth(simd_t const & src, uint8_t const mask)
{
    return detail::extract_impl<8>(src, mask & 0x07);
}

template <Simd simd_t>
//!\cond
    requires (detail::is_builtin_simd_v<simd_t>)
//!\endcond
constexpr bool test_all_ones(simd_t const & src)
{
    bool res = true;
    for (int8_t i = 0; i < static_cast<int8_t>(simd_traits<simd_t>::length); ++i)
        res &= (src[i] & ~typename simd_traits<simd_t>::scalar_type{0});

    return res;
}

/*!\brief Unpacks and interleaves the high halves of both simd vectors.
 * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[int] first   The vector whose values come before the `second`.
 * \param[int] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
// template <size_t scalar_size, Simd simd_t>
// constexpr simd_t unpack_hi(simd_t const & first, simd_t const & second)
// {
//     static_assert(std::Integral<typename simd_traits<simd_t>::scalar_type>,
//                   "Only integral scalar types are supported.");
//     return detail::unpack_hi<scalar_size>(first, second);
// }
//
// /*!\brief Unpacks and interleaves the high halves of both simd vectors.
//  * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
//  * \param[int] first   The vector whose values come before the `second`.
//  * \param[int] second  The vector whose values come after the `first`.
//  * \ingroup simd
//  */
// template <Simd simd_t>
// constexpr simd_t unpack_hi(simd_t const & first, simd_t const & second)
// {
//     static_assert(std::Integral<typename simd_traits<simd_t>::scalar_type>,
//                   "Only integral scalar types are supported.");
//     return unpack_hi<sizeof(typename simd_traits<simd_t>::scalar_type)>(first, second);
// }
//
// /*!\brief Unpacks and interleaves the low halves of both simd vectors.
//  * \tparam scalar_size The integer width for packing. Either 1, 2, 4, or 8.
//  * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
//  * \param[int] first   The vector whose values come before the `second`.
//  * \param[int] second  The vector whose values come after the `first`.
//  * \ingroup simd
//  */
// template <size_t scalar_size, Simd simd_t>
// constexpr simd_t unpack_lo(simd_t const & first, simd_t const & second)
// {
//     static_assert(std::Integral<typename simd_traits<simd_t>::scalar_type>,
//                   "Only integral scalar types are supported.");
//     return detail::unpack_lo<scalar_size>(first, second);
// }
//
// /*!\brief Unpacks and interleaves the low halves of both simd vectors.
//  * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
//  * \param[int] first   The vector whose values come before the `second`.
//  * \param[int] second  The vector whose values come after the `first`.
//  * \ingroup simd
//  */
// template <Simd simd_t>
// constexpr simd_t unpack_lo(simd_t const & first, simd_t const & second)
// {
//     static_assert(std::Integral<typename simd_traits<simd_t>::scalar_type>,
//                   "Only integral scalar types are supported.");
//     return unpack_lo<sizeof(typename simd_traits<simd_t>::scalar_type)>(first, second);
// }

} // inline namespace simd

} // namespace seqan3
