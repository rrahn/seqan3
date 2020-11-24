// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_element.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\interface seqan3::detail::alignment_matrix_element <>
 * \brief
 * \ingroup alignment_matrix
 *
 * This concept describes the requirements for a type that is used to determine the value and reference types of an
 * alignment matrix column.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 */

/*!\name Requirements for seqan3::detail::alignment_matrix_element
 * \brief You can expect these functions and type aliases on all types that model
 *        seqan3::detail::alignment_matrix_element.
 * \relates seqan3::detail::alignment_matrix_element
 *
 * \details
 *
 * In the following description:
 * * `t` is a type meeting the requirements of seqan3::detail::alignment_matrix_element
 * * `e` is a non-const value meeting the requirements of seqan3::detail::alignment_matrix_element
 *
 * \{
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT alignment_matrix_element = requires (t & e)
{
//!\endcond

    /*!\typedef typename t\:\:value_type;
     * \brief The underlying type of the element stored in the alignment matrix.
     */
    typename t::value_type;
    /*!\typedef typename t\:\:storage_value_type;
     * \brief The type of an entry in the matrix.
     */
    typename t::storage_value_type;
    /*!\typedef typename t\:\:element_type;
     * \brief The type of an entry in the matrix.
     */
    typename t::element_type;
    /*!\typedef typename t\:\:element_reference;
     * \brief The element reference type returned by the alignment matrix column iterator.
     *
     * \details
     *
     * This
     */
    typename t::element_reference;

    requires std::semiregular<typename t::value_type>;

    /*!\fn element_reference make_element(storage_value_type & value)
    * \brief Assign an ungapped sequence to a gapped sequence.
    *
    * \param[in] value The value of the underlying matrix storage that is transformed to the seqan3::detail::alignment_matrix_element::element_reference type.
    *
    * \details
    *
    * An aligned sequence has to be assignable from its unaligned counter part. For example a
    * std::vector<seqan3::gapped<seqan3::dna4>> as well as a seqan3::gap_decorator<std::vector<seqan3::dna4>>
    * can be assigned from s std::vector<seqan3::dna4> via seqan3::assign_unaligned.
    *
    * \attention This is a concept requirement, not an actual function (however types
    *            modelling this concept will provide an implementation).
    */
    //!\cond
    SEQAN3_RETURN_TYPE_CONSTRAINT(e.make_element(std::declval<typename t::storage_value_type &>()),
                                  std::same_as, typename t::element_reference);
    //!\endcond

    /*!\fn storage_value_type initialise(value_type & init)
    * \brief Assign an ungapped sequence to a gapped sequence.
    *
    * \param[in] init The value to initialise the matrix storage and the element itself.
    *
    * \details
    *
    * An aligned sequence has to be assignable from its unaligned counter part. For example a
    * std::vector<seqan3::gapped<seqan3::dna4>> as well as a seqan3::gap_decorator<std::vector<seqan3::dna4>>
    * can be assigned from s std::vector<seqan3::dna4> via seqan3::assign_unaligned.
    *
    * \attention This is a concept requirement, not an actual function (however types
    *            modelling this concept will provide an implementation).
    */
    //!\cond
    SEQAN3_RETURN_TYPE_CONSTRAINT(e.initialise(std::declval<typename t::value_type>()),
                                  std::same_as, typename t::storage_value_type);
    //!\endcond
//!\cond
};
//!\endcond
}  // namespace seqan3::detail
