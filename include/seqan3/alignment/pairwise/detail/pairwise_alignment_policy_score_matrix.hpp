// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_policy_score_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/affine_cell_proxy.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

template <typename score_t>
class single_column_score_matrix
{
private:
    using vertical_column_t = decltype(views::repeat_n(score_t{}, 1));

    std::vector<score_t> optimal_column{};
    std::vector<score_t> horizontal_column{};
    vertical_column_t vertical_column{};
    size_t columns_count{};

    class iterator
    {
    private:
        single_column_score_matrix * host_ptr{nullptr};
        size_t current_column{};
        size_t end_column{};

    public:

        /*!\name Associated types
         * \{
         */
        // using value_type = zip_column_t;
        // using reference = zip_column_t;
        using pointer = void;
        using difference = void;
        using iterator_tag = std::forward_iterator_tag;
        //!\}

        /*!\name Constructor, assignment and destructor
         * \{
         */
        iterator() noexcept = default; //!< Defaulted.
        iterator(iterator const &) noexcept = default; //!< Defaulted.
        iterator(iterator &&) noexcept = default; //!< Defaulted.
        iterator & operator=(iterator const &) noexcept = default; //!< Defaulted.
        iterator & operator=(iterator &&) noexcept = default; //!< Defaulted.
        ~iterator() = default; //!< Defaulted.

        explicit iterator(single_column_score_matrix & matrix_obj) noexcept :
            host_ptr{std::addressof(matrix_obj)},
            end_column{host_ptr->columns_count}
        {}
        //!\}

        auto operator*() const noexcept
        {
            return views::zip(host_ptr->optimal_column, host_ptr->horizontal_column, host_ptr->vertical_column)
                 | std::views::transform([] (auto && tpl)
                {
                    using fwd_tuple_t = decltype(tpl);
                    return affine_cell_proxy<remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
                });
        }

        iterator & operator++() noexcept
        {
            ++current_column;
            return *this;
        }

        iterator operator++(int) noexcept
        {
            iterator tmp{*this};
            ++(*this);
            return tmp;
        }

        friend bool operator==(iterator const & lhs, iterator const & rhs) noexcept
        {
            return lhs.current_column == rhs.current_column;
        }

        friend bool operator!=(iterator const & lhs, iterator const & rhs) noexcept
        {
            return !(lhs == rhs);
        }

        friend bool operator==(iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
        {
            return lhs.current_column == lhs.end_column;
        }

        friend bool operator==(std::ranges::default_sentinel_t const & lhs, iterator const & rhs) noexcept
        {
            return rhs == lhs;
        }

        friend bool operator!=(iterator const & lhs, std::ranges::default_sentinel_t const & rhs) noexcept
        {
            return !(lhs == rhs);
        }

        friend bool operator!=(std::ranges::default_sentinel_t const & lhs, iterator const & rhs) noexcept
        {
            return rhs != lhs;
        }
    };

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    single_column_score_matrix() = default; //!< Defaulted.
    single_column_score_matrix(single_column_score_matrix const &) = default; //!< Defaulted.
    single_column_score_matrix(single_column_score_matrix &&) = default; //!< Defaulted.
    single_column_score_matrix & operator=(single_column_score_matrix const &) = default; //!< Defaulted.
    single_column_score_matrix & operator=(single_column_score_matrix &&) = default; //!< Defaulted.
    ~single_column_score_matrix() = default; //!< Defaulted.
    //!\}

    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void reset_matrix(sequence1_t && seq1, sequence2_t && seq2)
    {
        size_t column_size = std::ranges::distance(seq2) + 1;
        optimal_column.clear();
        optimal_column.resize(column_size, 0);
        horizontal_column.resize(column_size, 0);
        vertical_column = views::repeat_n(score_t{}, column_size);
        columns_count = std::ranges::distance(seq1) + 1;
    }

    iterator begin()
    {
        return iterator{*this};
    }

    std::ranges::default_sentinel_t end()
    {
        return std::ranges::default_sentinel;
    }
};

template <typename alignment_config_t>
class pairwise_alignment_policy_score_matrix
{
private:
    using traits_type = alignment_configuration_traits<alignment_config_t>;
    using score_type = typename traits_type::score_t;
    using score_matrix_type = single_column_score_matrix<score_type>;

    score_matrix_type score_matrix{};

    struct alignment_matrix
    {
    private:
        score_matrix_type * score_matrix_ptr{};
        score_matrix_type local_score_matrix{};

    public:
        alignment_matrix() = default;
        alignment_matrix(alignment_matrix const &) = default;
        alignment_matrix(alignment_matrix &&) = default;
        alignment_matrix & operator=(alignment_matrix const &) = default;
        alignment_matrix & operator=(alignment_matrix &&) = default;

        alignment_matrix(score_matrix_type & global_matrix) : score_matrix_ptr{std::addressof(global_matrix)}
        { // Exchange the cached matrix with the local one.
            assert(score_matrix_ptr != nullptr);
            local_score_matrix = std::move(*score_matrix_ptr);
        }

        ~alignment_matrix()
        { // Restore global matrix to keep buffer reallocation to a minimum.
            assert(score_matrix_ptr != nullptr);
            *score_matrix_ptr = std::move(local_score_matrix);
        }

        auto begin()
        {
            return local_score_matrix.begin();
        }

        auto end()
        {
            return local_score_matrix.end();
        }
    };

public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    pairwise_alignment_policy_score_matrix() = default; //!< Defaulted.
    pairwise_alignment_policy_score_matrix(pairwise_alignment_policy_score_matrix const &) = default; //!< Defaulted.
    pairwise_alignment_policy_score_matrix(pairwise_alignment_policy_score_matrix &&) = default; //!< Defaulted.
    pairwise_alignment_policy_score_matrix & operator=(pairwise_alignment_policy_score_matrix const &) = default; //!< Defaulted.
    pairwise_alignment_policy_score_matrix & operator=(pairwise_alignment_policy_score_matrix &&) = default; //!< Defaulted.
    ~pairwise_alignment_policy_score_matrix() = default; //!< Defaulted.

    explicit pairwise_alignment_policy_score_matrix(alignment_config_t const & /*config*/)
    {
        // handle band parameters
    }
    //!\}

protected:
    /*!\brief Acquires a local alignment matrix to be used for the alignment algorithm.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] seq1 The first sequence.
     * \param[in] seq2 The second sequence.
     *
     * \details
     *
     * Resets the cached alignment matrix with the new dimensions given by the size of the passed sequences and returns
     * a seqan3::detail::pairwise_alignment_policy_score_matrix::alignment_matrix wrapping the actual matrix.
     * The returned matrix is a local variable for the algorithm but reuses the memory of the globally cached matrix
     * in order to avoid unecessary allocations in between calls to the same algorithm instance.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    alignment_matrix acquire_alignment_matrix(sequence1_t && seq1, sequence2_t && seq2)
    {
        score_matrix.reset_matrix(std::forward<sequence1_t>(seq1), std::forward<sequence2_t>(seq2));
        return alignment_matrix{score_matrix};
    }
};

} // namespace seqan3::detail
