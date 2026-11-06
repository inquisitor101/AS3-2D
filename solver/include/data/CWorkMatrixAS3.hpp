#pragma once

/*!
 * @file
 * @brief Borrowed matrix view with LIFO lifetime tracking.
 */


#include <cstddef>

#include "NError.hpp"
#include "common/AS3ScalarTypes.hpp"


/*!
 * @brief A non-owning row-major matrix borrowed from a LIFO workspace pool.
 *
 * @tparam T mutable built-in arithmetic scalar satisfying AS3MatrixScalar.
 *
 * @note Workspace matrices must be destroyed in reverse creation order. The backing pool must remain alive.
 */
template<AS3MatrixScalar T>
class CWorkMatrixAS3
{
	public:

		/*!
		 * @brief Constructor that wraps a borrowed pool block and records its LIFO identity.
		 *
		 * @param[in] data borrowed contiguous matrix storage; the pool retains ownership.
		 * @param[in] inuse address of the originating pool's entry counter; never a copy of its value.
		 * @param[in] ninst address of the originating pool's live-workspace counter.
		 * @param[in] r number of rows.
		 * @param[in] c number of columns.
		 */
		explicit CWorkMatrixAS3(T                *data,
		                       std::size_t      *inuse,
		                       std::size_t      *ninst,
		                       const std::size_t r,
		                       const std::size_t c)
			: mData(data), mInUse(inuse), mNInst(ninst), mNRow(r), mNCol(c), mID(*ninst) { }

		/*!
		 * @brief Destructor that returns this matrix block to its pool in LIFO order.
		 */
		~CWorkMatrixAS3(void) 
		{
			// These pointers refer to the pool's own counters, so every borrower updates the same state.
			// The pool must remain alive at its original address until all borrowers are destroyed.
			// This check ensures that the last-in first-out (LIFO) memory model is respected.
			if( mID != *mNInst ) ERROR("Object creation/deletion needs to abide by the LIFO convention.");
			*mInUse -= mNRow*mNCol; (*mNInst)--;
		}

		/*!
		 * @brief Return the total  size of the matrix.
		 *
		 * @return total number of scalar entries.
		 */
		std::size_t size(void) const {return mNRow*mNCol;}

		/*!
		 * @brief Return the row    size of the matrix.
		 *
		 * @return number of matrix rows.
		 */
		std::size_t row(void)  const {return mNRow;}

		/*!
		 * @brief Return the column size of the matrix.
		 *
		 * @return number of matrix columns.
		 */
		std::size_t col(void)  const {return mNCol;}

		/*!
		 * @brief Getter function which returns the pointer to the first stored entry.
		 *
		 * @return pointer to the first stored entry.
		 */
		T* data(void) 
		{
			return mData;
		}

		/*!
		 * @brief Getter function which returns a mutable reference to one flattened matrix entry.
		 *
		 * @param[in] index zero-based entry index.
		 *
		 * @return reference to the selected matrix entry.
		 */
		T& operator[](const std::size_t index)
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		/*!
		 * @brief Getter function which returns a const reference to one flattened matrix entry.
		 *
		 * @param[in] index zero-based entry index.
		 *
		 * @return reference to the selected matrix entry.
		 */
		const T& operator[](const std::size_t index) const
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		/*!
		 * @brief Getter function which returns a mutable reference to one matrix entry.
		 *
		 * @param[in] r zero-based row index.
		 * @param[in] c zero-based column index.
		 *
		 * @return reference to the selected matrix entry.
		 */
		T& operator()(const std::size_t r, const std::size_t c)
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}

		/*!
		 * @brief Getter function which returns a const reference to one matrix entry.
		 *
		 * @param[in] r zero-based row index.
		 * @param[in] c zero-based column index.
		 *
		 * @return reference to the selected matrix entry.
		 */
		const T& operator()(const std::size_t r, const std::size_t c) const
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}

	private:
		T           *mData  = nullptr; ///< Borrowed matrix entries; ownership remains with the pool.
		std::size_t *mInUse = nullptr; ///< Shared count of scalar entries currently borrowed from the pool.
		std::size_t *mNInst = nullptr; ///< Shared count of live workspace matrices.
		std::size_t  mNRow  = 0;       ///< Number of matrix rows.
		std::size_t  mNCol  = 0;       ///< Number of matrix columns.
		std::size_t  mID;              ///< LIFO instance identifier assigned when this matrix was borrowed.
		
		
		

		// Disable default constructor.
		CWorkMatrixAS3(void) = delete;
		// Disable default copy constructor.
		CWorkMatrixAS3(const CWorkMatrixAS3&) = delete;
		// Disable default copy operator.
		CWorkMatrixAS3& operator=(CWorkMatrixAS3&) = delete;
};
