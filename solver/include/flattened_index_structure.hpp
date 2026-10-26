#pragma once

/*!
 * @file
 * @brief Zone-local element addresses and category-local face addresses.
 */


#include <cstddef>

#include "option_structure.hpp"


/*!
 * @brief A face category and its category-local batch and face indices.
 */
struct CFlattenedFaceIndex
{
  std::size_t            mIndexBatch;  ///< Batch index local to the selected face category.
  std::size_t            mIndexFace;   ///< Face index local to its batch.
  ETypeFaceGeometry mFaceType;         ///< Face category selecting the interior, boundary or interface collection.
};


/*!
 * @brief The zone index and zone-local index of one flattened element.
 */
struct CFlattenedElementIndex
{
	std::size_t mIndexZone;          ///< Index of the containing zone.
	std::size_t         mIndexElem;  ///< Element index local to the zone.
};


