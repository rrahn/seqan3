// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_matrix_storage.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/detail/alignment_matrix_element_concept.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>

namespace seqan3::detail
{
/*!\interface seqan3::detail::alignment_matrix_storage <>
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

/*!\name Requirements for seqan3::detail::alignment_matrix_storage
 * \brief You can expect these functions and type aliases on all types that model
 *        seqan3::detail::alignment_matrix_storage.
 * \relates seqan3::detail::alignment_matrix_storage
 *
 * \details
 *
 * In the following description:
 * * `t` is a type meeting the requirements of seqan3::detail::alignment_matrix_storage
 * * `s` is a non-const value meeting the requirements of seqan3::detail::alignment_matrix_storage
 *
 * \{
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT alignment_matrix_storage =
    requires (t & s, row_index_type<size_t> rdim, column_index_type<size_t> cdim, size_t i)
{
//!\endcond

    /*!\typedef typename t\:\:value_type;
     * \brief The underlying type of the element stored in the alignment matrix.
     */
    typename t::column_type;

    /*!\typedef typename t\:\:value_type;
     * \brief The underlying type of the element stored in the alignment matrix.
     */
    typename t::matrix_element_type;

    //!\brief The element type must meet the requirements of seqan3::detail::alignment_matrix_element
    requires alignment_matrix_element<typename t::matrix_element_type>;

    /*!\fn element_reference make_element(storage_value_type & value)
    * \brief Assign an ungapped sequence to a gapped sequence.
    *
    * \param[in] value The value of the underlying matrix storage that is transformed to the seqan3::detail::alignment_matrix_storage::element_reference type.
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
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.resize(rdim, cdim), std::same_as, void);
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.resize(rdim, cdim, std::declval<typename t::matrix_element_type::value_type>()),
                                  std::same_as, void);
    //!\endcond

    //!\cond
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.clear(), std::same_as, void);
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.cols(), std::same_as, size_t);
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.rows(), std::same_as, size_t);
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
    SEQAN3_RETURN_TYPE_CONSTRAINT(s.column_at(i), std::same_as, typename t::column_type);
    //!\endcond
//!\cond
};
//!\endcond
}  // namespace seqan3::detail
