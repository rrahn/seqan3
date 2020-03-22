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

#include <seqan3/core/platform.hpp>

#include <seqan3/std/span>

namespace seqan3::detail
{

template <typename score_t>
class single_column_score_matrix
{
private:
    using combined_score_t = std::pair<score_t, score_t>;

    std::vector<combined_score_t> optimal_column{};
    size_t number_columns{};

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
        using value_type = decltype(host_ptr->optimal_column);
        using reference = std::span<combined_score_t>;
        using const_reference = std::span<combined_score_t const>;
        using pointer = void;
        using differnce = void;
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
            end_column{host_ptr->number_columns}
        {}
        //!\}

        reference operator*() const noexcept
        {
            return std::span{host_ptr->optimal_column};
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
        optimal_column.clear();
        optimal_column.resize(std::ranges::distance(seq2) + 1, {0, 0});
        number_columns = std::ranges::distance(seq1) + 1;
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
     * The returned matrix is a local variable for the algorithm but reuses the memory of the globally cahed matrix
     * in order to avoid unecessary allocations inbetween calls to the same algorithm instance.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    alignment_matrix acquire_alignment_matrix(sequence1_t && seq1, sequence2_t && seq2)
    {
        score_matrix.reset_matrix(std::forward<sequence1_t>(seq1), std::forward<sequence2_t>(seq2));
        return alignment_matrix{score_matrix};
    }
};

} // namespace seqan3::detail
