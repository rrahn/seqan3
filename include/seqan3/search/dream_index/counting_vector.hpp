// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::counting_vector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <seqan3/std/bit>
#include <seqan3/std/concepts>
#include <type_traits>
#include <vector>

#include <sdsl/bit_vectors.hpp>

#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/simd/simd.hpp>

namespace seqan3::detail
{

/*!\brief The base type for the counting vector.
 * \ingroup submodule_dream_index
 * \tparam value_t The value type to use.
 * \tparam simd_enabled A bool value wether the simd is enabled.
 *
 * \details
 *
 * Selects the correct vector type to be used as the base class for the seqan3::counting_vector.
 */
template <typename value_t, bool simd_enabled>
using counting_vector_base_t =
    std::vector<value_t, seqan3::aligned_allocator<value_t,
                                                   (simd_enabled) ? alignof(simd_type_t<value_t>)
                                                                  : alignof(std::max_align_t)>>;


} // namespace seqan3::detail

namespace seqan3
{

/*!\brief A data structure that behaves like a std::vector and can be used to consolidate the results of multiple calls
 *        to seqan3::interleaved_bloom_filter::membership_agent::bulk_contains.
 * \ingroup submodule_dream_index
 * \tparam value_t The type of the count. Must model std::integral.
 *
 * \details
 *
 * When using the seqan3::interleaved_bloom_filter::membership_agent::bulk_contains operation, a common use case is to
 * add up, for example, the results for all k-mers in a query. This yields, for each bin, the number of k-mers of a
 * query that are in the respective bin. Such information can be used to apply further filtering or abundance estimation
 * based on the k-mer counts.
 *
 * The seqan3::counting_vector offers an easy way to add up the individual
 * seqan3::interleaved_bloom_filter::membership_agent::binning_bitvector by offering an `+=` operator.
 *
 * The `value_t` template parameter should be chosen in a way that no overflow occurs if all calls to `bulk_contains`
 * return a hit for a specific bin. For example, `uint8_t` will suffice when processing short Illumina reads, whereas
 * long reads will require at least `uint32_t`.
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_vector.cpp
 */
template <std::integral value_t, bool simd_enabled = false>
class counting_vector : public detail::counting_vector_base_t<value_t, simd_enabled>
{
private:
    //!\brief The base type.
    using base_t = detail::counting_vector_base_t<value_t, simd_enabled>;

    //!\brief A flag indicating wether the AVX512BW instruction set is available.
    static constexpr bool has_avx512bw =
    #ifdef __AVX512BW__
        true;
    #else
        false;
    #endif // __AVX512BW__

    //!\brief A flag indicating wether the AVX512F instruction set is available.
    static constexpr bool has_avx512f =
    #ifdef __AVX512F__
        true;
    #else
        false;
    #endif // __AVX512F__

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_vector() = default; //!< Defaulted.
    counting_vector(counting_vector const &) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector const &) = default; //!< Defaulted.
    counting_vector(counting_vector &&) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector &&) = default; //!< Defaulted.
    ~counting_vector() = default; //!< Defaulted.

    using base_t::base_t;
    //!\}

    /*!\brief Bin-wise adds the bits of a seqan3::interleaved_bloom_filter::membership_agent::binning_bitvector.
     * \tparam rhs_t The type of the right-hand side.
     *         Must be seqan3::interleaved_bloom_filter::membership_agent::binning_bitvector.
     * \param rhs The seqan3::interleaved_bloom_filter::membership_agent::binning_bitvector.
     * \attention The counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    template <typename rhs_t>
    //!\cond
        requires std::is_base_of<sdsl::bit_vector, rhs_t>::value && !simd_enabled
    //!\endcond
    counting_vector & operator+=(rhs_t const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        // Each iteration can handle 64 bits, so we need to iterate `((rhs.size() + 63) >> 6` many times
        for (size_t batch = 0, bin = 0; batch < ((rhs.size() + 63) >> 6); bin = 64 * ++batch)
        {
            size_t tmp = rhs.get_int(batch * 64); // get 64 bits starting at position `batch * 64`
            if (tmp ^ (1ULL<<63)) // This is a special case, because we would shift by 64 (UB) in the while loop.
            {
                while (tmp > 0)
                {
                    // Jump to the next 1 and increment the corresponding vector entry.
                    uint8_t step = std::countr_zero(tmp);
                    bin += step++;
                    tmp >>= step;
                    ++(*this)[bin++];
                }
            }
            else
            {
                ++(*this)[bin + 63];
            }
        }
        return *this;
    }

    //!\brief overload
    template <typename rhs_t>
    //!\cond
        requires std::is_base_of<sdsl::bit_vector, rhs_t>::value && simd_enabled
    //!\endcond
    counting_vector & operator+=(rhs_t const & rhs)
    {
        static_assert((sizeof(value_t) <= 2) ? has_avx512bw : has_avx512f,
                      "Simd counting was requested but the required instructions set was not available. "
                      "Requires AVX512BW for one byte and two byte counters or AVX512F for larger counters.");

        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        // Each iteration can handle 64 bits, so we need to iterate `((rhs.size() + 63) >> 6` many times
        for (size_t batch = 0; batch < ((rhs.size() + 63) >> 6); ++batch)
        {
            size_t const bin_position = batch * 64;
            uint64_t const bit_mask = rhs.get_int(bin_position);
            value_t * data_ptr = this->data() + bin_position;

            mask_increment(data_ptr, bit_mask);
        }
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan3::counting_vector`s.
     * \param rhs The other seqan3::counting_vector.
     * \attention The seqan3::counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }

private:
    /*!\brief Uses the bit mask and increments the corresponding counter starting at the data pointer.
     * \param data_ptr The pointer to the first counter to increment.
     * \param bit_mask A 64 bit wide mask indicating which of the 64 counters should be incremented.
     *
     * \details
     *
     * Uses _mm512_mask_add_epiX instructions to increment the subsequent 64 counters starting at `data_ptr`,
     * depending on the bits set in the 64 bit mask.
     * If the size of the couner type is bigger than 1 byte the corresponding mask_add operation is repeated
     * until all 64 subsequent counters have been processed.
     */
    void mask_increment(value_t * data_ptr, uint64_t const bit_mask) noexcept
        requires simd_enabled
    {
        using simd_t = simd_type_t<value_t>;

        constexpr size_t simd_length{simd_traits<simd_t>::length};
        constexpr simd_t one_vector{simd::fill<simd_t>(1)};

        for (uint32_t shift = 0; shift < sizeof(value_t); ++shift)
        {
            uint32_t const epi_shift = shift * simd_length;
            uint64_t const epi_bit_mask = bit_mask >> epi_shift;
            __m512i * src = reinterpret_cast<__m512i *>(data_ptr + epi_shift);

            if constexpr (simd_length == 64) // epi8 - requires AVX512BW
                *src = _mm512_mask_add_epi8(*src, epi_bit_mask, *src, reinterpret_cast<__m512i const &>(one_vector));
            else if constexpr (simd_length == 32) // epi16 - requires AVX512BW
                *src = _mm512_mask_add_epi16(*src, epi_bit_mask, *src, reinterpret_cast<__m512i const &>(one_vector));
            else if constexpr (simd_length == 16) // epi32 - requires AVX512F
                *src = _mm512_mask_add_epi32(*src, epi_bit_mask, *src, reinterpret_cast<__m512i const &>(one_vector));
            else // epi64 - requires AVX512F
                *src = _mm512_mask_add_epi64(*src, epi_bit_mask, *src, reinterpret_cast<__m512i const &>(one_vector));
        }
    }
};

}  // namespace seqan3
