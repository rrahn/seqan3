// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::interleaved_bloom_filter.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/bit>

#include <sdsl/bit_vectors.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/dream_index/counting_vector.hpp>

namespace seqan3
{

/*!\addtogroup submodule_dream_index
 * \{
 */
//!\brief Determines if the Interleaved Bloom Filter is compressed.
enum data_layout : bool
{
    uncompressed, //!< The Interleaved Bloom Filter is uncompressed.
    compressed    //!< The Interleaved Bloom Filter is compressed.
};

//!\brief A strong type that represents the number of bins for the seqan3::interleaved_bloom_filter.
struct bin_count : public detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the number of bits for each bin in the seqan3::interleaved_bloom_filter.
struct bin_size : public detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the number of hash functions for the seqan3::interleaved_bloom_filter.
struct hash_function_count : public detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the bin index for the seqan3::interleaved_bloom_filter.
struct bin_index : public detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>::strong_type;
};

/*!\brief The IBF binning directory. A data structure that efficiently answers set-membership queries for multiple bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See seqan3::data_layout.
 * \implements seqan3::cerealisable
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) that can answer these queries efficiently.
 *
 * ### Interleaved Bloom Filter (IBF)
 *
 * The Interleaved Bloom Filter is a probabilistic data structure that extends the
 * [Bloom Filter](https://en.wikipedia.org/wiki/Bloom_filter).
 * A Bloom Filter can be thought of as a bitvector of length `n` and `h` hash functions and is used to determine set
 * membership. To insert data, the data is hashed by the `h` hash functions (returning values in `[0, n)`) and the
 * corresponding `h` positions in the bitvector are set to `1`. To query data, i.e. to determine whether the query
 * belongs to the set the Bloom Filter was built for, the query is hashed by the same `h` hash functions and the
 * corresponding positions are checked. If all `h` positions contain a `1`, the query is (probably) in the data set.
 * Since the Bloom Filter has variable length, the hashing is not bijective, i.e. it may return true for a set
 * membership query even though the query was never inserted into the Bloom Filter. Note that the Bloom Filter
 * will always return `true` if the query was inserted, i.e. there may be false positives, but no false negatives.
 *
 * The Interleaved Bloom Filter now applies the concept of a Bloom Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in `b` data sets/bins.
 * Conceptually, a Bloom Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting `b` Bloom Filters are then interleaved such that the `i`'th bit if each Bloom Filter are
 * adjacent to each other:
 * ```
 * Bloom Filter 0       Bloom Filter 1      Bloom Filter 2      Bloom Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where `x.y` denotes the `y`'th bit of the `x`'th Bloom Filter.
 * ```
 * Interleaved Bloom Filter
 * |0.0|1.0|2.0|3.0|0.1|1.1|2.1|3.1|0.2|1.2|2.2|3.2|0.3|1.3|2.3|3.3|
 * ```
 * A query can now be searched in all `b` bins by computing the `h` hash functions, retrieving the `h` sub-bitvectors of
 * length `b` starting at the positions indicated by the hash functions. The bitwise AND of these sub-bitvectors yields
 * the binningvector, a bitvector of length `b` where the `i`'th bit indicates set membership in the `i`'th bin.
 *
 * ### Querying
 * To query the Interleaved Bloom Filter for a value, call seqan3::interleaved_bloom_filter::membership_agent() and use
 * the returned seqan3::interleaved_bloom_filter::membership_agent.
 *
 * To count the occurrences of a range of values in the Interleaved Bloom Filter, call
 * seqan3::interleaved_bloom_filter::counting_agent() and use
 * the returned seqan3::interleaved_bloom_filter::counting_agent_type.
 *
 * ### Compression
 *
 * The Interleaved Bloom Filter can be compressed by passing `data_layout::compressed` as template argument.
 * The compressed `seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>` can only be constructed from a
 * `seqan3::interleaved_bloom_filter`, in which case the underlying bitvector is compressed.
 * The compressed Interleaved Bloom Filter is immutable, i.e. only querying is supported.
 *
 * ### Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * Additionally, concurrent calls to `emplace` are safe iff each thread handles a multiple of wordsize (=64) many bins.
 * For example, calls to `emplace` from multiple threads are safe if `thread_1` accesses bins 0-63, `thread_2` bins
 * 64-127, and so on.
 */
template <data_layout data_layout_mode_ = data_layout::uncompressed>
class interleaved_bloom_filter
{
private:
    //!\cond
    template <data_layout data_layout_mode>
    friend class interleaved_bloom_filter;
    //!\endcond

    //!\brief The underlying datatype to use.
    using data_type = std::conditional_t<data_layout_mode_ == data_layout::uncompressed,
                                         sdsl::bit_vector,
                                         sdsl::sd_vector<>>;

    //!\brief The number of bins specified by the user.
    size_t bins{};
    //!\brief The number of bins stored in the IBF (next multiple of 64 of `bins`).
    size_t technical_bins{};
    //!\brief The size of each bin in bits.
    size_t bin_size_{};
    //!\brief The number of bits to shift the hash value before doing multiplicative hashing.
    size_t hash_shift{};
    //!\brief The number of 64-bit integers needed to store `bins` many bits (e.g. `bins = 50` -> `bin_words = 1`).
    size_t bin_words{};
    //!\brief The number of hash functions.
    size_t hash_funs{};
    //!\brief The bitvector.
    data_type data{};
    //!\brief Precalculated seeds for multiplicative hashing. We use large irrational numbers for a uniform hashing.
    static constexpr std::array<size_t, 5> hash_seeds{13572355802537770549ULL, // 2**64 / (e/2)
                                                      13043817825332782213ULL, // 2**64 / sqrt(2)
                                                      10650232656628343401ULL, // 2**64 / sqrt(3)
                                                      16499269484942379435ULL, // 2**64 / (sqrt(5)/2)
                                                      4893150838803335377ULL}; // 2**64 / (3*pi/5)

    /*!\brief Perturbs a value and fits it into the vector.
     * \param h The value to process.
     * \param seed The seed to use.
     * \returns A hashed value representing a position within the bounds of `data`.
     * \sa https://probablydance.com/2018/06/16/
     * \sa https://lemire.me/blog/2016/06/27
     */
    inline constexpr size_t hash_and_fit(size_t h, size_t const seed) const
    {
        h *= seed;
        assert(hash_shift < 64);
        h ^= h >> hash_shift; // XOR and shift higher bits into lower bits
        h *= 11400714819323198485ULL; // = 2^64 / golden_ration, to expand h to 64 bit range
        // Use fastrange (integer modulo without division) if possible.
#ifdef __SIZEOF_INT128__
        h = static_cast<uint64_t>((static_cast<__uint128_t>(h) * static_cast<__uint128_t>(bin_size_)) >> 64);
#else
        h %= bin_size_;
#endif
        h *= technical_bins;
        return h;
    }

public:
    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr data_layout data_layout_mode = data_layout_mode_;

    class membership_agent; // documented upon definition below
    template <std::integral value_t, bool = false>
    class counting_agent_type; // documented upon definition below

    //!\brief The counting agent type in simd mode.
    template <std::integral value_t>
    using simd_counting_agent_type = counting_agent_type<value_t, true>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    interleaved_bloom_filter() = default; //!< Defaulted.
    interleaved_bloom_filter(interleaved_bloom_filter const &) = default; //!< Defaulted.
    interleaved_bloom_filter & operator=(interleaved_bloom_filter const &) = default; //!< Defaulted.
    interleaved_bloom_filter(interleaved_bloom_filter &&) = default; //!< Defaulted.
    interleaved_bloom_filter & operator=(interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~interleaved_bloom_filter() = default; //!< Defaulted.

    /*!\brief Construct an uncompressed Interleaved Bloom Filter.
     * \param bins_ The number of bins.
     * \param size The bitvector size.
     * \param funs The number of hash functions. Default 2. At least 1, at most 5.
     *
     * \attention This constructor can only be used to construct **uncompressed** Interleaved Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_constructor.cpp
     */
    interleaved_bloom_filter(seqan3::bin_count bins_,
                             seqan3::bin_size size,
                             seqan3::hash_function_count funs = seqan3::hash_function_count{2u})
    //!\cond
        requires (data_layout_mode == data_layout::uncompressed)
    //!\endcond
    {
        bins = bins_.get();
        bin_size_ = size.get();
        hash_funs = funs.get();

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (hash_funs == 0 || hash_funs > 5)
            throw std::logic_error{"The number of hash functions must be > 0 and <= 5."};
        if (bin_size_ == 0)
            throw std::logic_error{"The size of a bin must be > 0."};

        hash_shift = std::countl_zero(bin_size_);
        bin_words = (bins + 63) >> 6; // = ceil(bins/64)
        technical_bins  = bin_words << 6; // = bin_words * 64
        data = sdsl::bit_vector(technical_bins * bin_size_);
    }

    /*!\brief Construct a compressed Interleaved Bloom Filter.
     * \param[in] ibf The uncompressed seqan3::interleaved_bloom_filter.
     *
     * \attention This constructor can only be used to construct **compressed** Interleaved Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_constructor_compressed.cpp
     */
    interleaved_bloom_filter(interleaved_bloom_filter<data_layout::uncompressed> const & ibf)
    //!\cond
        requires (data_layout_mode == data_layout::compressed)
    //!\endcond
    {
        std::tie(bins, technical_bins, bin_size_, hash_shift, bin_words, hash_funs) =
            std::tie(ibf.bins, ibf.technical_bins, ibf.bin_size_, ibf.hash_shift, ibf.bin_words, ibf.hash_funs);

        data = sdsl::sd_vector<>{ibf.data};
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Inserts a value into a specific bin.
     * \param[in] value The raw numeric value to process.
     * \param[in] bin The bin index to insert into.
     *
     * \attention This function is only available for **uncompressed** Interleaved Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_emplace.cpp
     */
    void emplace(size_t const value, bin_index const bin)
    //!\cond
        requires (data_layout_mode == data_layout::uncompressed)
    //!\endcond
    {
        assert(bin.get() < bins);
        for (size_t i = 0; i < hash_funs; ++i)
        {
            size_t idx = hash_and_fit(value, hash_seeds[i]);
            idx += bin.get();
            assert(idx < data.size());
            data[idx] = 1;
        };
    }

    /*!\brief Increases the number of bins stored in the Interleaved Bloom Filter.
     * \param[in] new_bins_ The new number of bins.
     * \throws std::invalid_argument If passed number of bins is smaller than current number of bins.
     *
     * \attention This function is only available for **uncompressed** Interleaved Bloom Filters.
     * \attention The new number of bins must be greater or equal to the current number of bins.
     * \attention This function invalidates all seqan3::interleaved_bloom_filter::membership_agent constructed for this
     * Interleaved Bloom Filter.
     *
     * \details
     *
     * The resulting `seqan3::interleaved_bloom_filter` has an increased size proportional to the increase in the
     * `bin_words` (the number of 64-bit words needed to represent `bins` many bins), e.g.
     * resizing a `seqan3::interleaved_bloom_filter` with 40 bins to 73 bins also increases the `bin_words` from 1 to
     * 2 and hence the new `seqan3::interleaved_bloom_filter` will be twice the size.
     * This increase in size is necessary to avoid invalidating all computed hash functions.
     * If you want to add more bins while keeping the size constant, you need to rebuild the
     * `seqan3::interleaved_bloom_filter`.
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_increase_bin_number_to.cpp
     */
    void increase_bin_number_to(bin_count const new_bins_)
    //!\cond
        requires (data_layout_mode == data_layout::uncompressed)
    //!\endcond
    {
        size_t new_bins = new_bins_.get();

        if (new_bins < bins)
            throw std::invalid_argument{"The number of new bins must be >= the current number of bins."};

        // Equivalent to ceil(new_bins / 64)
        size_t new_bin_words = (new_bins + 63) >> 6;

        bins = new_bins;

        if (new_bin_words == bin_words) // No need for internal resize if bin_words does not change.
            return;

        size_t new_technical_bins = new_bin_words << 6;
        size_t new_bits = bin_size_ * new_technical_bins;

        size_t idx_{new_bits}, idx{data.size()};
        size_t delta = new_technical_bins - technical_bins + 64;

        data.resize(new_bits);

        for (size_t i = idx_, j = idx; j > 0; i -= new_technical_bins, j -= technical_bins)
        {
            size_t stop = i - new_technical_bins;

            for (size_t ii = i - delta, jj = j - 64; stop && ii >= stop; ii -= 64, jj -= 64)
            {
                uint64_t old = data.get_int(jj);
                data.set_int(jj, 0);
                data.set_int(ii, old);
            }
        }

        bin_words = new_bin_words;
        technical_bins = new_technical_bins;
    }
    //!\}

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns a seqan3::interleaved_bloom_filter::membership_agent to be used for lookup.
     * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan3::interleaved_bloom_filter::membership_agent`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_construction.cpp
     * \sa seqan3::interleaved_bloom_filter::membership_agent::bulk_contains
     */
    membership_agent membership_agent() const
    {
        return typename interleaved_bloom_filter<data_layout_mode>::membership_agent{*this};
    }

    /*!\brief Returns a seqan3::interleaved_bloom_filter::counting_agent_type to be used for counting.
     * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan3::interleaved_bloom_filter::counting_agent_type`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent_construction.cpp
     * \sa seqan3::interleaved_bloom_filter::counting_agent_type::bulk_count
     */
    template <typename value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }

    /*!\brief Returns a seqan3::interleaved_bloom_filter::simd_counting_agent_type to be used for counting.
     * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan3::interleaved_bloom_filter::simd_counting_agent_type`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * This simd counting agent allows to use simd instructions to increment the counters.
     * It requires the AVX512 instruction set. If the binary is compiled with another instruction set,
     * then a static assertion will be raised.
     *
     * \sa seqan3::interleaved_bloom_filter::counting_agent
     */
    template <typename value_t = uint16_t>
    simd_counting_agent_type<value_t> simd_counting_agent() const
    {
        static_assert(simd_traits<simd_type_t<value_t>>::max_length == 64,
                      "The simd counting agent requires the AVX512 instruction set! "
                      "Please make sure that the binary is compiled with the correct instruction set.");

        return simd_counting_agent_type<value_t>{*this};
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    /*!\brief Returns the number of hash functions used in the Interleaved Bloom Filter.
     * \returns The number of hash functions.
     */
    size_t hash_function_count() const noexcept
    {
        return hash_funs;
    }

    /*!\brief Returns the number of bins that the Interleaved Bloom Filter manages.
     * \returns The number of bins.
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }

    /*!\brief Returns the size of a single bin that the Interleaved Bloom Filter manages.
     * \returns The size in bits of a single bin.
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return data.size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Test for equality.
     * \param[in] lhs A `seqan3::interleaved_bloom_filter`.
     * \param[in] rhs `seqan3::interleaved_bloom_filter` to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    friend bool operator==(interleaved_bloom_filter const & lhs, interleaved_bloom_filter const & rhs) noexcept
    {
        return std::tie(lhs.bins, lhs.technical_bins, lhs.bin_size_, lhs.hash_shift, lhs.bin_words, lhs.hash_funs,
                        lhs.data) ==
               std::tie(rhs.bins, rhs.technical_bins, rhs.bin_size_, rhs.hash_shift, rhs.bin_words, rhs.hash_funs,
                        rhs.data);
    }

    /*!\brief Test for inequality.
     * \param[in] lhs A `seqan3::interleaved_bloom_filter`.
     * \param[in] rhs `seqan3::interleaved_bloom_filter` to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(interleaved_bloom_filter const & lhs, interleaved_bloom_filter const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(bins);
        archive(technical_bins);
        archive(bin_size_);
        archive(hash_shift);
        archive(bin_words);
        archive(hash_funs);
        archive(data);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the seqan3::interleaved_bloom_filter.
 * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to on `ibf` invalidates the
 * membership_agent.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/membership_agent_construction.cpp
 */
template <data_layout data_layout_mode>
class interleaved_bloom_filter<data_layout_mode>::membership_agent
{
private:
    //!\brief The type of the augmented seqan3::interleaved_bloom_filter.
    using ibf_t = interleaved_bloom_filter<data_layout_mode>;

    // std::vector<uint64_t> _subvector_buffer{};
    std::vector<uint64_t> _bloom_filter_indices{};
    //!\brief A pointer to the augmented seqan3::interleaved_bloom_filter.
    ibf_t const * ibf_ptr{nullptr};
    uint32_t _hash_function_count{};

public:
    class binning_bitvector;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default; //!< Defaulted.
    membership_agent(membership_agent const &) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default; //!< Defaulted.
    ~membership_agent() = default; //!< Defaulted.

    /*!\brief Construct a membership_agent from a seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    membership_agent(ibf_t const & ibf) :
        ibf_ptr{std::addressof(ibf)},
        _hash_function_count{static_cast<uint32_t>(ibf.hash_funs)}
    {
        result_buffer_vector.resize(1);
        result_buffer_vector.front().resize(ibf_ptr->bin_count());
    };
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector<binning_bitvector> result_buffer_vector;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of a given value.
     * \param[in] value The raw value to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_bloom_filter::membership_agent for each thread.
     */
    [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) & noexcept
    {
        assert(result_buffer_vector.size() >= 1);
        binning_bitvector & result_buffer = result_buffer_vector.front();

        assert(ibf_ptr != nullptr);
        assert(result_buffer.size() == ibf_ptr->bin_count());

        std::array<size_t, 5> bloom_filter_indices;
        std::memcpy(&bloom_filter_indices, &ibf_ptr->hash_seeds, sizeof(size_t) * ibf_ptr->hash_funs);

        for (size_t i = 0; i < ibf_ptr->hash_funs; ++i)
            bloom_filter_indices[i] = ibf_ptr->hash_and_fit(value, bloom_filter_indices[i]);

        for (size_t batch = 0; batch < ibf_ptr->bin_words; ++batch)
        {
            size_t tmp{-1ULL};
           for (size_t i = 0; i < ibf_ptr->hash_funs; ++i)
           {

                assert(bloom_filter_indices[i] < ibf_ptr->data.size());
                tmp &= ibf_ptr->data.get_int(bloom_filter_indices[i]);
                bloom_filter_indices[i] += 64;
           }

            result_buffer.set_int(batch << 6, tmp);
        }
        return result_buffer;  // Checks which bins are set per k-mer.
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    [[nodiscard]] binning_bitvector const & bulk_contains(size_t const & value) && noexcept = delete;

    //!\overload
    template <typename hash_range_t>
    [[nodiscard]] std::vector<binning_bitvector> const & bulk_contains(hash_range_t && hash_range) & noexcept
    {
        assert(ibf_ptr != nullptr);
        // assert(result_buffer.size() == ibf_ptr->bin_count());
        uint32_t const hash_range_size = std::ranges::distance(hash_range); // TODO: maybe not sized range?

        result_buffer_vector.resize(hash_range_size, binning_bitvector{ibf_ptr->bin_count()});
        auto result_buffer_it = result_buffer_vector.begin();

        // Step 1: Split hash value computation and ibf lookup.
        std::vector<uint64_t> bloom_filter_indices{};
        bloom_filter_indices.resize(hash_range_size * ibf_ptr->hash_funs);
        auto bloom_filter_indices_it = bloom_filter_indices.begin();
        for (size_t const hash_value : hash_range)
        { // Nice to vectorise!
            for (size_t i = 0; i < ibf_ptr->hash_funs; ++i, ++bloom_filter_indices_it)
                *bloom_filter_indices_it = ibf_ptr->hash_and_fit(hash_value, ibf_t::hash_seeds[i]);
        }

        // Step 2: Run ibf lookup and update bitvector.
        std::vector<size_t> tmp_subvector{};
        tmp_subvector.resize(ibf_ptr->bin_words);

        for (bloom_filter_indices_it = bloom_filter_indices.begin();
             bloom_filter_indices_it != bloom_filter_indices.end();)
        {
            std::ranges::fill(tmp_subvector, -1ULL);

            for (size_t i = 0; i < ibf_ptr->hash_funs; ++i, ++bloom_filter_indices_it)
            {
                for (size_t & tmp : tmp_subvector)
                {
                    assert(*bloom_filter_indices_it < ibf_ptr->data.size());
                    tmp &= ibf_ptr->data.get_int(*bloom_filter_indices_it);
                    *bloom_filter_indices_it += 64;
                }
            }

            // Store the result_buffer
            for (size_t batch = 0; batch < tmp_subvector.size(); ++batch)
                result_buffer_it->set_int(batch << 6, tmp_subvector[batch]);

            ++result_buffer_it;
        }

        return result_buffer_vector;  // Checks which bins are set per k-mer.
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <typename hash_range_t>
    [[nodiscard]] std::vector<binning_bitvector> const & bulk_contains(hash_range_t const &) && noexcept = delete;
    //!\}

};

//!\brief A bitvector representing the result of a call to `bulk_contains` of the seqan3::interleaved_bloom_filter.
template <data_layout data_layout_mode>
class interleaved_bloom_filter<data_layout_mode>::membership_agent::binning_bitvector : private sdsl::bit_vector
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_bitvector() = default; //!< Defaulted.
    binning_bitvector(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector(binning_bitvector &&) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector &&) = default; //!< Defaulted.
    ~binning_bitvector() = default; //!< Defaulted.

    binning_bitvector(size_t new_size) : sdsl::bit_vector(new_size, 0)
    {}
    //!\}

    using sdsl::bit_vector::begin;
    using sdsl::bit_vector::end;
    using sdsl::bit_vector::operator==;
    using sdsl::bit_vector::operator[];
    using sdsl::bit_vector::size;

#if SEQAN3_DOXYGEN_ONLY(1)0
    //!\brief The iterator type of the `binning_bitvector`;
    using iterator_t = IMPLEMENTATION_DEFINED;
    //!\brief The reference type of the `binning_bitvector`;
    using reference_t = IMPLEMENTATION_DEFINED;
    //!\brief The const_reference type of the `binning_bitvector`;
    using const_reference_t = IMPLEMENTATION_DEFINED;
    //!\brief Returns an iterator to the begin of the bitvector.
    iterator_t begin() noexcept;
    //!\brief Returns an iterator to the end of the bitvector.
    iterator_t end() noexcept;
    //!\brief Compares two bitvectors.
    bool operator==(bit_vector const & other) const noexcept;
    //!\brief Returns a reference to position `idx` of the bitvector.
    reference_t operator[](size_t const & idx) noexcept;
    //!\brief Returns a const_reference to position `idx` of the bitvector.
    const_reference_t operator[](size_t const & idx) const noexcept;
    //!\brief Returns the size of the bitvector.
    size_t size() noexcept;
#endif

private:
    friend class membership_agent;
    using sdsl::bit_vector::resize;
    using sdsl::bit_vector::set_int;
    template <std::integral value_t, bool>
    friend class counting_vector;
};

/*!\brief Manages counting ranges of values for the seqan3::interleaved_bloom_filter.
 * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <data_layout data_layout_mode>
template <std::integral value_t, bool simd_enabled>
class interleaved_bloom_filter<data_layout_mode>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::interleaved_bloom_filter.
    using ibf_t = interleaved_bloom_filter<data_layout_mode>;
    //!\brief The type of the counting vector depending on the simd mode.
    using counting_vector_type = counting_vector<value_t, simd_enabled>;

    //!\brief A pointer to the augmented seqan3::interleaved_bloom_filter.
    ibf_t const * ibf_ptr{nullptr};

    //!\brief Store a seqan3::interleaved_bloom_filter::membership_agent to call `bulk_contains`.
    typename ibf_t::membership_agent membership_agent;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    counting_agent_type(ibf_t const & ibf) : ibf_ptr(std::addressof(ibf)), membership_agent(ibf)
    {
        result_buffer.resize(ibf_ptr->bin_count());
    };
    //!\}

    //!\brief Stores the result of bulk_count().
    counting_vector_type result_buffer{};

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector_type const & bulk_count(value_range_t && values) & noexcept
    {
        assert(ibf_ptr != nullptr);
        assert(result_buffer.size() == ibf_ptr->bin_count());

        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, 0);

        // Step 1: get all hashes:
        auto const & membership_vector = membership_agent.bulk_contains(values);

        for (auto && bit_vector : membership_vector)
            result_buffer += bit_vector;

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector_type const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}

};

//!\}

} // namespace seqan3
