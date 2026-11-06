#pragma once

/*!
 * @file
 * @brief Workspace pool for borrowed LIFO matrices.
 */


#include <cstddef>

#include "NError.hpp"
#include "data/CWorkMatrixAS3.hpp"


/*!
 * @brief A contiguous workspace pool that lends matrices in last-in, first-out order.
 *
 * @tparam T mutable built-in arithmetic scalar satisfying AS3MatrixScalar.
 *
 * @note All borrowed matrices must be released before the pool is destroyed. Each OpenMP thread owns its own
 *   pool. Pools stay at a fixed address because borrowed matrices refer to their counters.
 */
template<AS3MatrixScalar T>
class CPoolMatrixAS3
{
	public:

		/*!
		 * @brief Constructor that allocates and zero-initializes the workspace pool.
		 *
		 * @param[in] nsize total number of scalar entries to allocate.
		 */
		explicit CPoolMatrixAS3(std::size_t nsize) 
			: mNSize(nsize), mData(new T[nsize]{})
		{
			// Inline counters avoid two heap allocations and cannot leak if the buffer allocation throws.
		}

		/*!
		 * @brief Destructor that releases the pool after all borrowed matrices have been destroyed.
		 */
		~CPoolMatrixAS3(void)
		{
			if( mNInst ) ERROR("Some workspace matrices still exist.");
			delete[] mData;
		}

		/*!
		 * @brief Return a specific memory block, wrapped as a matrix.
		 *
		 * @param[in] r number of rows.
		 * @param[in] c number of columns.
		 *
		 * @return borrowed matrix, released automatically in LIFO order.
		 */
		CWorkMatrixAS3<T> GetWorkMatrixAS3(const std::size_t r, const std::size_t c) 
		{
			const std::size_t I0 = mInUse;
			// Check remaining storage before r*c or I0+n can wrap around.
			// A zero-row matrix borrows no entries, regardless of its column count.
			if( r && c > (mNSize - I0)/r ) ERROR("Requested size exceeds allocated memory.");
			const std::size_t n = r*c;
			mInUse = I0 + n;
			mNInst++;

			// Issue (via RVO) a work matrix, based on the specified dimensions.
			// Borrowers still store pointers; only the counters' ownership/storage changed.
			// Passing member addresses keeps their destruction updates visible to this pool.
			return CWorkMatrixAS3<T>( mData + I0, &mInUse, &mNInst, r, c );
		}

		/*!
		 * @brief Return the total size of the memory pool.
		 *
		 * @return total number of scalar entries.
		 */
		std::size_t size(void)      const {return  mNSize;}

		/*!
		 * @brief Return the number of object instances in use.
		 *
		 * @return number of live borrowed matrices.
		 */
		std::size_t ninst(void)     const {return mNInst;}

		/*!
		 * @brief Return the currently used memory of the memory pool.
		 *
		 * @return number of scalar entries currently borrowed.
		 */
		std::size_t inuse(void)     const {return mInUse;}

		/*!
		 * @brief Return the available memory for usage.
		 *
		 * @return number of scalar entries available for borrowing.
		 */
		std::size_t available(void) const {return (mNSize - mInUse);}


	private:
		std::size_t mNSize = 0;        ///< Total number of scalar entries allocated in the pool.
		std::size_t mInUse = 0;        ///< Borrowed entry count; work matrices point to this member.
		std::size_t mNInst = 0;        ///< Live workspace count; work matrices point to this member.
		T          *mData  = nullptr;  ///< Owned contiguous pool storage.

		// Disable default constructor.
		CPoolMatrixAS3(void) = delete;
		// Disable default copy constructor.
		CPoolMatrixAS3(const CPoolMatrixAS3&) = delete;
		// Disable default copy operator.
		CPoolMatrixAS3& operator=(const CPoolMatrixAS3&) = delete;
		// Moving would invalidate the counter pointers held by outstanding workspace matrices.
		CPoolMatrixAS3(CPoolMatrixAS3&&) = delete;
		CPoolMatrixAS3& operator=(CPoolMatrixAS3&&) = delete;
};
