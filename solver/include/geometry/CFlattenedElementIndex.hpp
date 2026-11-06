#pragma once

/*!
 * @file
 * @brief Zone-local address of one flattened element.
 */


#include <cstddef>


/*!
 * @brief The zone index and zone-local index of one flattened element.
 */
struct CFlattenedElementIndex
{
	std::size_t mIndexZone;  ///< Index of the containing zone.
	std::size_t mIndexElem;  ///< Element index local to the zone.
};
