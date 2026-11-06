#pragma once

/*!
 * @file
 * @brief Category-local address of one flattened face.
 */


#include <cstddef>

#include "common/AS3Options.hpp"


/*!
 * @brief A face category and its category-local batch and face indices.
 */
struct CFlattenedFaceIndex
{
  std::size_t       mIndexBatch;  ///< Batch index local to the selected face category.
  std::size_t       mIndexFace;   ///< Face index local to its batch.
  ETypeFaceGeometry mFaceType;    ///< Face category selecting the interior, boundary or interface collection.
};
