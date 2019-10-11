// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::edit_distance_algorithm.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>

namespace seqan3::detail
{

/*!\brief This algorithm unifies different edit distance implementations and uses the appropriate one depending on the
 * given configuration.
 * \implements std::invocable
 * \tparam config_t The configuration type.
 * \tparam traits_t The traits type.
 *
 * \details
 *
 * This wrapper class is used to decouple the sequence types from the algorithm class type.
 * Within the alignment configuration a std::function object storing this wrapper is returned
 * if an edit distance should be computed. On invocation it delegates the call to the actual implementation
 * of the edit distance algorithm, while the interface is unified with the execution model of the pairwise alignment
 * algorithms.
 */
template <typename config_t, typename traits_t>
class edit_distance_algorithm
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr edit_distance_algorithm() = default;                                            //!< Defaulted
    constexpr edit_distance_algorithm(edit_distance_algorithm const &) = default;             //!< Defaulted
    constexpr edit_distance_algorithm(edit_distance_algorithm &&) = default;                  //!< Defaulted
    constexpr edit_distance_algorithm & operator=(edit_distance_algorithm const &) = default; //!< Defaulted
    constexpr edit_distance_algorithm & operator=(edit_distance_algorithm &&) = default;      //!< Defaulted
    ~edit_distance_algorithm() = default;                                                     //!< Defaulted

    /*!\brief Constructs the wrapper with the passed configuration.
     * \param cfg The configuration to be passed to the algorithm.
     *
     * \details
     *
     * The configuration is copied once to the heap during construction and maintained by a std::shared_ptr.
     * The configuration is not passed to the function-call-operator of this function object in order to avoid
     * incompatible configurations between the passed configuration and the one used during configuration of this
     * class. Further, the function object will be stored in a std::function which requires copyable objects and
     * in parallel executions the function object must be copied as well.
     */
    constexpr edit_distance_algorithm(config_t const & cfg) : cfg_ptr{new config_t(cfg)}
    {}
    //!}

    /*!\brief Invokes the actual alignment computation given two sequences.
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::forward_range.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::forward_range.
     * \param[in] idx            The index of the current sequence pair.
     * \param[in] first_range    The first sequence (or packed sequences).
     * \param[in] second_range   The second sequence (or packed sequences).
     */
    template <std::ranges::forward_range first_range_t, std::ranges::forward_range second_range_t>
    constexpr [[nodiscard]] auto operator()(size_t const idx,
                                            first_range_t && first_range,
                                            second_range_t && second_range)
    {
        using edit_traits = default_edit_distance_trait_type<first_range_t,
                                                             second_range_t,
                                                             config_t,
                                                             typename traits_t::is_semi_global_type>;
        edit_distance_unbanded algo{first_range, second_range, *cfg_ptr, edit_traits{}};
        return algo(idx);
    }

private:
    //!\brief The alignment configuration stored on the heap.
    std::shared_ptr<remove_cvref_t<config_t>> cfg_ptr{};
};

} // namespace seqan3::detail
