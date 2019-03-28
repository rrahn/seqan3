// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::to_simd view.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

template <std::ranges::View urng_t, Simd simd_t>
//!\cond
    requires std::ranges::ForwardRange<urng_t>
//!\endcond
class view_to_simd_chunk : public std::ranges::view_interface<view_to_simd_chunk<urng_t, simd_t>>
{
    using inner_range_t = value_type_t<urng_t>;

    static_assert(std::ranges::InputRange<inner_range_t>, "Expects a range of ranges as underlying range.");
    static_assert(Semialphabet<value_type_t<inner_range_t>>, "Expects semi-alphabet as value type.");

    static constexpr bool fast_load         = std::ranges::ContiguousRange<inner_range_t> &&
                                              std::SizedSentinel<std::ranges::iterator_t<inner_range_t>,
                                                                 std::ranges::sentinel_t<inner_range_t>> &&
                                              sizeof(alphabet_rank_t<value_type_t<inner_range_t>>) == 1 &&
                                              simd_traits<simd_t>::max_length < 64;  // Currently not supported on AVX512

    //!\brief The size of one chunk. Equals the number of elements in the vector.
    static constexpr int8_t chunk_size      = simd_traits<simd_t>::length;
    //!\brief The number of chunks that can be gathered with a single load.
    static constexpr int8_t chunks_per_load = simd_traits<simd_t>::max_length / chunk_size;
    //!\brief The total number of chunks that can be cached.
    static constexpr int8_t total_chunks    = (fast_load) ? (chunks_per_load * chunks_per_load) : 1;
    //!\brief Helper index sequence to
    static constexpr auto   idx_seq         = std::make_index_sequence<chunk_size>{};

    using chunk_t                           = std::array<simd_t, chunk_size>;
    using scalar_type                       = typename simd_traits<simd_t>::scalar_type;
    using max_simd_t                        = typename simd_type<int8_t>::type;

    class iterator_type
    {
    public:

        /*!\name Associated types
         * \{
         */
        //!\brief The reference type.
        using reference = chunk_t &;
        //!\brief The value type.
        using value_type = chunk_t;
        //!\brief The pointer type.
        using pointer = chunk_t *;
        //!\brief The difference type.
        using difference_type = ptrdiff_t;
        //!\brief The iterator category;
        using iterator_category = std::input_iterator_tag;
        //!\}

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr iterator_type()                                  = default; //!< Defaulted.
        constexpr iterator_type(iterator_type const &)             = default; //!< Defaulted.
        constexpr iterator_type(iterator_type &&)                  = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type const &) = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type &&)      = default; //!< Defaulted.
        ~iterator_type()                                           = default; //!< Defaulted.

        constexpr iterator_type(view_to_simd_chunk & _this_view) : this_view{&_this_view}, cached_chunk_pos{0}
        {
            // Initialise the iterator of the sub ranges.
            final_chunk = true;
            auto i_it = cached_iter.begin();
            auto i_sentinel = cached_sentinel.begin();
            for (auto & irng : this_view->urng)
            {
                *i_it = std::ranges::begin(irng);
                *i_sentinel = std::ranges::end(irng);
                final_chunk = (final_chunk) ? (*i_it == *i_sentinel) : final_chunk;
                ++i_it;
                ++i_sentinel;
            }

            load_next_chunk();
        }
        //!\}

        constexpr reference operator*() const noexcept
        {
            assert(this_view != nullptr);
            return this_view->cached_simd_chunks[cached_chunk_pos];
        }

        constexpr iterator_type & operator++(/*pre-increment*/)
        {
            if constexpr (fast_load)
            {
                if (cached_chunk_pos == final_chunk_pos)
                {
                    load_next_chunk();
                    cached_chunk_pos = 0;
                }
                else
                {
                    ++cached_chunk_pos;
                }
            }
            else
            {
                load_next_chunk();
            }

            return *this;
        }

        constexpr void operator++(int /*post-increment*/)
        {
            ++(*this);
        }

        /*!\name Comparison operators
         * \{
         */
        //!\brief Checks for equality with sentinel.
        constexpr bool operator==(std::ranges::default_sentinel_t const &) const noexcept
        {
            return at_end;
        }

        friend constexpr bool operator==(std::ranges::default_sentinel_t const & /*lhs*/,
                                         iterator_type const & rhs) noexcept
        {
            return rhs.at_end;
        }

        constexpr bool operator!=(std::ranges::default_sentinel_t const & /*rhs*/) const noexcept
        {
            return !at_end;
        }

        friend constexpr bool operator!=(std::ranges::default_sentinel_t const & /*lhs*/,
                                         iterator_type const & rhs) noexcept
        {
            return !rhs.at_end;
        }
        //!\}

    private:

        auto unpack(max_simd_t const & row)
        {
            if constexpr (chunk_size == simd_traits<max_simd_t>::length / 2)  // upcast into 2 vectors.
            {
                return std::array{simd::upcast_auto<simd_t>(row),                        // 1. part
                                  simd::upcast_auto<simd_t>(extract_halve(row, 0x01))};  // 2. part
            }
            else if constexpr (chunk_size == simd_traits<max_simd_t>::length / 4) // upcast into 4 vectors.
            {
                return std::array{simd::upcast_auto<simd_t>(row),                        // 1. part
                                  simd::upcast_auto<simd_t>(extract_quarter(row, 0x1)),  // 2. part
                                  simd::upcast_auto<simd_t>(extract_quarter(row, 0x2)),  // 3. part
                                  simd::upcast_auto<simd_t>(extract_quarter(row, 0x3))}; // 4. part
            }
            else if constexpr (chunk_size == simd_traits<max_simd_t>::length / 8) // upcast into 8 vectors.
            {
                return std::array{simd::upcast_auto<simd_t>(row),                         // 1. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x01)),   // 2. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x02)),   // 3. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x03)),   // 4. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x04)),   // 5. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x05)),   // 6. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x06)),   // 7. part
                                  simd::upcast_auto<simd_t>(extract_eighth(row, 0x07))};  // 8. part
            }
            else  // no cast necessary.
            {
                static_assert(chunk_size == simd_traits<max_simd_t>::length, "Not supported simd type.");
                return reinterpret_cast<simd_t>(row);
            }
        }

        constexpr void unpack_and_cache(std::array<max_simd_t, simd_traits<max_simd_t>::length> matrix)
        {
            // Iterate over the rows of the matrix
            for (int8_t i = 0; i < static_cast<int8_t>(matrix.size()); ++i)
            {
                // unpack depending on chunks_per_load.
                auto res = unpack(matrix[i]);
                // debug_stream << "DEBUG " << __LINE__ << ": "<<res[0] << " " << res[1] << "\n";

                if constexpr (std::Same<decltype(res), simd_t>)
                {
                    static_assert(simd_traits<max_simd_t>::length == chunk_size, "Expected byte packed simd type.");
                    this_view->cached_simd_chunks[0][i] = std::move(res);
                }
                else // We need to parse the tuple elements and store them in the cached simd chunks.
                {
                    static_assert(res.size() == chunks_per_load, "Expected chunks_per_load many simd vectors.");

                    // int8_t block = i / chunks_per_load;
                    for (int8_t j = 0; j < chunks_per_load; ++j)  // store chunks in respective cached entries.
                    {
                        size_t idx = (j * simd_traits<max_simd_t>::length + i) / simd_traits<simd_t>::length;
                        this_view->cached_simd_chunks[idx][i % simd_traits<simd_t>::length] =
                            std::move(res[j]);
                    }
                }
            }
        }

        template <size_t ... indices>
        constexpr bool
        all_iterators_reached_sentinel(std::index_sequence<indices...> const & SEQAN3_DOXYGEN_ONLY(idx_seq)) noexcept
        {
            return (true && ... && (cached_iter[indices] == cached_sentinel[indices]));
        }

        constexpr void load_next_chunk()
            requires fast_load
        {
            at_end = final_chunk;
            if (at_end)  // reached end of stream.
                return;
            // For the efficient load we assume at most one byte sized alphabets.
            // Hence we can load max_simd_t length many elements at once.
            // Depending on the packing of simd_t we can prefetch blocks and store them in the cached_simd_chunks.
            // E.g. assume simd_t with length 8 on SSE4 with max length 16.
            // To fill the 16x16 matrix we need four 8x8 matrices.
            // Thus, for the 8 sequences we need to load two times 16 consecutive bytes to fill the matrix.
            // This quadratic byte matrix can be transposed efficiently with simd instructions.

            constexpr int8_t max_size = simd_traits<max_simd_t>::length;
            constexpr int8_t num_chunks = max_size / chunk_size;
            std::array<max_simd_t, max_size> matrix{};
            final_chunk_pos = 0;  // reset the final chunk position, since this could be the last load.
            // Iterate over each sequence.
            for (int8_t i = 0; i < chunk_size; ++i)
            {  // Iterate over each block depending on the packing of the target simd vector.
                int8_t last_chunk_pos = -1;
                for (int8_t j = 0; j < num_chunks; ++j)
                {
                    int8_t pos = j * chunk_size + i; // matrix entry to fill
                    if (cached_sentinel[i] - cached_iter[i] >= max_size) // not in final block
                    {
                        matrix[pos] = simd::load<max_simd_t>(std::addressof(*cached_iter[i]));
                        std::advance(cached_iter[i], max_size);
                        last_chunk_pos += num_chunks; // We have read num_chunks at once.
                    }
                    else  // Loads the final block byte wise in order to not load from uninitialised memory.
                    {
                        auto & it = cached_iter[i];
                        auto & sent = cached_sentinel[i];
                        max_simd_t & tmp = matrix[pos];
                        tmp = this_view->cached_default;
                        for (int8_t idx = 0; it != sent; ++it, ++idx)
                        {
                            tmp[idx] = seqan3::to_rank(*it);
                            if (idx % num_chunks == 0)
                                ++last_chunk_pos; // increment whenever we store an element from the num_chunks boundary
                        }
                    }
                }
                final_chunk_pos = std::max(final_chunk_pos, last_chunk_pos);
            }
            final_chunk = all_iterators_reached_sentinel(idx_seq);
            simd::transpose(matrix);

            unpack_and_cache(std::move(matrix));
        }

        // Iterate for the current simd vector over the column.
        template <size_t ... indices>
        constexpr void fill_simd(chunk_t & target,
                                 size_t const pos,
                                 std::index_sequence<indices...> const & SEQAN3_DOXYGEN_ONLY(idx_seq))
        {
            auto get_and_update = [this, &pos] (size_t const idx)
            {
                if (cached_iter[idx] == cached_sentinel[idx])
                {
                    if (this_view->padding_value.has_value())
                        return static_cast<scalar_type>(this_view->padding_value.value());
                    else
                        return this_view->cached_simd_chunks[cached_chunk_pos][pos][idx];
                }
                else
                {  // only increment if not at end.
                    return static_cast<scalar_type>(to_rank(*(cached_iter[idx]++)));
                }
            };
            target[pos] = simd::set<simd_t>(get_and_update(indices)...);
        }

        constexpr void load_next_chunk()
            requires !fast_load
        {
            at_end = final_chunk;
            if (at_end)  // reached end of stream.
                return;

            for (size_t i = 0; i < chunk_size; ++i)
                fill_simd(this_view->cached_simd_chunks[0], i, idx_seq);

            final_chunk = all_iterators_reached_sentinel(idx_seq);
        }

        std::array<std::ranges::iterator_t<inner_range_t>, chunk_size> cached_iter{};
        std::array<std::ranges::sentinel_t<inner_range_t>, chunk_size> cached_sentinel{};
        view_to_simd_chunk *                                           this_view{nullptr};
        int8_t                                                         final_chunk_pos{0};
        int8_t                                                         cached_chunk_pos{0};
        bool                                                           final_chunk{false};
        bool                                                           at_end{false};
    };

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr view_to_simd_chunk()                                       = default; //!< Defaulted.
    constexpr view_to_simd_chunk(view_to_simd_chunk const &)             = default; //!< Defaulted.
    constexpr view_to_simd_chunk(view_to_simd_chunk &&)                  = default; //!< Defaulted.
    constexpr view_to_simd_chunk & operator=(view_to_simd_chunk const &) = default; //!< Defaulted.
    constexpr view_to_simd_chunk & operator=(view_to_simd_chunk &&)      = default; //!< Defaulted.
    ~view_to_simd_chunk()                                                = default; //!< Defaulted.

    //!\brief Construction from the underlying view.
    constexpr view_to_simd_chunk(urng_t _urng, std::optional<scalar_type> const _padding_value = std::nullopt) :
        urng{std::move(_urng)},
        padding_value{std::move(_padding_value)}
    {
        // Check if the size is less or equal the simd size.
        if (std::ranges::distance(urng) > chunk_size)
            throw std::invalid_argument{"The size of the underlying range must be less than or equal to the size of "
                                        "the given simd type!"};

        if constexpr (fast_load)
        {
            if (padding_value.has_value()) // Leave uninitialised if not set.
                cached_default = simd::fill<max_simd_t>(padding_value.value());
        }
    }

    //!\brief Construction from std::ranges::ViewableRange.
    template <typename other_urng_t>
    //!\cond
    requires !std::Same<remove_cvref_t<other_urng_t>, view_to_simd_chunk> &&
             std::ranges::ViewableRange<other_urng_t> &&  // Must come after self type check to avoid conflicts with the move constructor.
             std::Constructible<urng_t, ranges::ref_view<std::remove_reference_t<other_urng_t>>>
    //!\endcond
    constexpr view_to_simd_chunk(other_urng_t && _urng, std::optional<scalar_type> const _padding_value = std::nullopt) :
        view_to_simd_chunk{std::view::all(_urng), std::move(_padding_value)}
    {}
    //!\}

    constexpr iterator_type begin() noexcept
    {
        return {*this};
    }

    constexpr void begin() const noexcept = delete;
    constexpr void cbegin() const noexcept = delete;

    constexpr std::ranges::default_sentinel_t end() noexcept
    {
        return std::ranges::default_sentinel;
    }

    constexpr void end() const noexcept = delete;
    constexpr void cend() const noexcept = delete;

    constexpr size_t size() const noexcept
    //!\cond
        requires std::ranges::SizedRange<inner_range_t>
    //!\endcond
    {
        size_t max_size = 0;
        for (auto & rng : urng)
            max_size = std::max(max_size, (std::ranges::size(rng) + chunk_size - 1) / chunk_size);

        return max_size;
    }

private:

    urng_t                            urng;
    std::array<chunk_t, total_chunks> cached_simd_chunks;
    max_simd_t                        cached_default{};
    std::optional<scalar_type>        padding_value;
};

// /*!\name Type deduction guides
//  * \{
//  */
//
// /*!\brief Type deduction guide for View inputs.
// * \relates seqan3::detail::view_to_simd_chunk
// */
// template <std::ranges::View urng_t, Simd simd_t>
// view_to_simd_chunk(urng_t, simd_t) -> view_to_simd_chunk<urng_t, simd_t>;
//
// /*!\brief Type deduction guide for ViewableRange inputs.
//  * \relates seqan3::detail::view_to_simd_chunk
//  */
// template <std::ranges::ViewableRange urng_t, Simd simd_t>
// view_to_simd_chunk(urng_t &&, size_t, simd_t)
//     -> view_to_simd_chunk<decltype(view::all(std::declval<urng_t>())), simd_t>;
//
// //!\}

// ============================================================================
//  to_simd_chunk_fn (adaptor definition)
// ============================================================================

/*!\brief view::to_simd_chunk's range adaptor object type (non-closure).
 */
template <Simd simd_t>
struct to_simd_chunk_fn
{
    using padding_t = typename simd_traits<simd_t>::scalar_type;
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(padding_t const padding) const noexcept
    {
        return detail::adaptor_from_functor{*this, std::optional{padding}};
    }

    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()() const noexcept
    {
        return detail::adaptor_from_functor{*this, std::optional<padding_t>{std::nullopt}};
    }

    /*!\brief            Call the view's constructor with the underlying view as argument.
     * \param[in] urange The input range to process. Must model std::ranges::ForwardRange and std::ranges::ViewableRange.
     * \param[in] i      The inserted range to process. Must model std::ranges::ForwardRange.
     * \param[in] size   The step size for insertion into the input range.
     * \returns          A range of with the inserted range interleaved into the underlying range at the specified intervals.
     */
    template <std::ranges::Range urng_t>
    constexpr auto operator()(urng_t && urange, std::optional<padding_t> const padding = std::nullopt) const noexcept
    {
        static_assert(std::ranges::ForwardRange<urng_t>,
            "The underlying range parameter in view::interleave must model std::ranges::ForwardRange.");
        static_assert(std::ranges::ViewableRange<urng_t>,
            "The underlying range parameter in view::interleave must model std::ranges::ViewableRange.");

        return view_to_simd_chunk<std::ranges::all_view<urng_t>, simd_t>{std::forward<urng_t>(urange), padding};
    }

    //!\brief This adaptor is usable without setting a padding value.
    template <std::ranges::Range urng_t>
    constexpr friend auto operator|(urng_t && urange, to_simd_chunk_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};


} // namespace seqan3::detail

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief A view that interleaves a given range into another range at regular intervals.
 * \tparam urng_t The type of the range being processed.
 * \tparam inserted_rng_t The type of the range being inserted.
 * \param[in] urange The range being processed.
 * \param[in] inserted_range The range being inserted.
 * \param[in] step_size A value of size_type which indicates the interval to insert the inserted_range.
 * \returns A range with the second range inserted at regular intervals. See below for properties of said range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to insert one range into another range at regular intervals. It behaves essentially like
 * `| std::view::chunk(step_size) | std::view::join(inserted_range)` except that for input that models
 * std::ranges::RandomAccessRange and std::ranges::SizedRange a more efficient data structure is returned
 * (otherwise it returns exactly the above combination of views).
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/interleave.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |---------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                     |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                     |
 * | std::ranges::BidirectionalRange | *required*                            | *preserved*                     |
 * | std::ranges::RandomAccessRange  | *required*                            | *preserved*                     |
 * | std::ranges::ContiguousRange    |                                       | *lost*                          |
 * |                                 |                                       |                                 |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                    |
 * | std::ranges::View               |                                       | *guaranteed*                    |
 * | std::ranges::SizedRange         | *required*                            | *preserved*                     |
 * | std::ranges::CommonRange        |                                       | *preserved*                     |
 * | std::ranges::OutputRange        |                                       | *preserved*                     |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>     |
 *
 *
 * If above requirements are not met, this adaptor forwards to
 * `| ranges::view::chunk(step_size) | ranges::view::join(inserted_range)`
 * which returns a view with the following properties:
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |---------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                     |
 * | std::ranges::ForwardRange       | *required*                            | *lost*                          |
 * | std::ranges::BidirectionalRange |                                       | *lost*                          |
 * | std::ranges::RandomAccessRange  |                                       | *lost*                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                          |
 * |                                 |                                       |                                 |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                    |
 * | std::ranges::View               |                                       | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                       | *lost*                          |
 * | std::ranges::CommonRange        |                                       | *lost*                          |
 * | std::ranges::OutputRange        |                                       | *lost*                          |
 * | seqan3::ConstIterableRange      |                                       | *lost*                          |
 * |                                 |                                       |                                 |
 * | seqan3::reference_t             |                                       | seqan3::value_type_t<urng_t>    |
 *
 * * `urng_t` is the type of the range modified by this view (input).
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref view.
 *
 * ### Example
 *
 * \include test/snippet/range/view/interleave.cpp
 * \hideinitializer
 */

template <Simd simd_t>
inline constexpr auto to_simd_chunk = detail::to_simd_chunk_fn<simd_t>{};

//!\}
} // namespace seqan3::view
