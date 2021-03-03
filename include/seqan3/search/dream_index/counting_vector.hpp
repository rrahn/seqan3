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
template <std::integral value_t>
class counting_vector : public std::vector<value_t>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t>;
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
        requires std::is_base_of<sdsl::bit_vector, rhs_t>::value
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
};

}  // namespace seqan3
