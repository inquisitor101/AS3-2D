#pragma once

/*!
 * @file
 * @brief Contiguous owning matrices and LIFO workspace matrices.
 */


#include <algorithm>
#include <concepts>
#include <cstddef>

#include <cstring>
#include <limits>
#include <type_traits>
#include "error_structure.hpp"


/*!
 * @brief Mutable built-in arithmetic scalars supported by all AS3 matrix containers.
 *
 * These types have numeric zero value-initialization, non-throwing element lifetimes,
 * and bytewise copy semantics. Trivially copyable alone would also admit pointers
 * and custom classes with throwing constructors. Const/volatile entries cannot
 * satisfy the containers' ordinary mutable storage and memcpy contract.
 * This constraint is checked at compile time and adds no runtime work.
 */
template<typename T>
concept AS3MatrixScalar = (std::integral<T> || std::floating_point<T>) &&
                           std::same_as<T, std::remove_cv_t<T>>;


/*!
 * @brief An owning row-major matrix with contiguous entries and zero-initialized allocations.
 *
 * @tparam T mutable built-in arithmetic scalar satisfying AS3MatrixScalar.
 *
 * @note Resizing resets all active entries to zero and retains allocated capacity when shrinking.
 */
template<AS3MatrixScalar T>
class CMatrixAS3
{
	public:
		

		/*!
		 * @brief Default constructor that initializes an empty CMatrixAS3.
		 */
		CMatrixAS3(void) = default;

		/*!
		 * @brief Constructor that initializes CMatrixAS3 from the specified data.
		 *
		 * @param[in] nrow number of rows.
		 * @param[in] ncol number of columns; defaults to one.
		 */
		explicit CMatrixAS3(std::size_t nrow, std::size_t ncol = 1)
			: mNRow(nrow), mNCol(ncol), mCapacity(CheckedSize(nrow, ncol)),
			  mData(mCapacity ? new T[mCapacity]{} : nullptr)
		{
			// Value-initialize scalars in the allocation; throwing new[] never returns nullptr.
		}

		/*!
		 * @brief Copy constructor, does a deep copy of same class.
		 *
		 * @param[in] other source matrix or object.
		 */
		CMatrixAS3(const CMatrixAS3<T> &other)
			: mNRow(other.mNRow), mNCol(other.mNCol), mCapacity(other.size()),
			  mData(mCapacity ? new T[mCapacity] : nullptr)
		{
			// Empty matrices can have null data pointers, so skip even zero-byte memcpy calls.
			if( mCapacity ) std::memcpy(mData, other.mData, mCapacity*sizeof(T));
		}

		/*!
		 * @brief Move constructor.
		 *
		 * @param[in,out] other source matrix or object.
		 */
		CMatrixAS3(CMatrixAS3&& other) noexcept
		{
			this->mNRow     = other.mNRow;
			this->mNCol     = other.mNCol;
			this->mCapacity = other.mCapacity;
			this->mData     = other.mData;
			other.mNRow     = 0;
			other.mNCol     = 0;
			other.mCapacity = 0;
			other.mData     = nullptr;
		}


		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CMatrixAS3(void)
		{
			if( mData )
			{
				delete[] mData;
				mData = nullptr;	
			}
		}

		/*!
		 * @brief Copy assignment operator.
		 *
		 * @param[in] other source matrix or object.
		 *
		 * @return reference to this matrix.
		 */
		CMatrixAS3<T>& operator=(const CMatrixAS3<T> &other)
		{
			if( this != &other )
			{
				const std::size_t count = other.size();
				if( count > mCapacity )
				{
					// Allocate first: if new[] throws, this matrix still owns its original buffer.
					T* data = new T[count];
					delete[] mData;
					mData = data;
					mCapacity = count;
				}
				// Reuse storage whenever it fits, including assignments with a different shape.
				if( count ) std::memcpy(mData, other.mData, count*sizeof(T));
				mNRow = other.mNRow;
				mNCol = other.mNCol;
			}
			return *this;
		}

		/*!
		 * @brief Move assignment operator.
		 *
		 * @param[in,out] other source matrix or object.
		 *
		 * @return reference to this matrix.
		 */
		CMatrixAS3<T>& operator=(CMatrixAS3<T> &&other) noexcept
		{
			if( this != &other )
			{
				if( mData ) delete[] mData;
				this->mData     = other.mData;
				this->mNRow     = other.mNRow;
				this->mNCol     = other.mNCol;
				this->mCapacity = other.mCapacity;
				other.mData     = nullptr;
				other.mNRow     = 0;
				other.mNCol     = 0;
				other.mCapacity = 0;
			}
			return *this;
		}


		/*!
		 * @brief Function that resizes the matrix and resets every entry to zero.
		 *
		 * @note Allocated capacity is retained when shrinking; existing entries are not preserved.
		 *
		 * @param[in] nrow number of rows.
		 * @param[in] ncol number of columns; defaults to one.
		 */
		void resize(std::size_t nrow, std::size_t ncol = 1)
		{
			const std::size_t count = CheckedSize(nrow, ncol);
			// Logical size may have shrunk; compare against the allocation, not the old shape.
			if( count > mCapacity )
			{
				T* data = new T[count];
				delete[] mData;
				mData = data;
				mCapacity = count;
			}
			mNRow = nrow;
			mNCol = ncol;
			reset(); // Resets all values to zero.
		}

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
		 * @brief Getter function which returns the pointer to the first stored entry.
		 *
		 * @return pointer to the first stored entry.
		 */
		const T* data(void) const 
		{
			return mData;
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
		 * @brief Reset values to zero.
		 */
		void reset(void)
		{
			for(std::size_t i=0; i<mNRow*mNCol; i++) mData[i] = static_cast<T>(0);
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

		/*!
		 * @brief Begin a mutable traversal of the contiguous entries in row-major order.
		 * @return pointer to the first entry; equals end() when empty.
		 */
    T* begin() noexcept
    {
      return mData;
    }
    
		/*!
		 * @brief End a mutable traversal of the active entries.
		 * @return pointer past the last entry; equals begin() when empty.
		 */
    T* end() noexcept
    {
      const std::size_t n = size();
      return n ? mData + n : mData;
    }
    
		/*!
		 * @brief Begin a read-only traversal of the contiguous entries in row-major order.
		 * @return pointer to the first entry; equals end() when empty.
		 */
    const T* begin() const noexcept
    {
      return mData;
    }
    
		/*!
		 * @brief End a read-only traversal of the active entries.
		 * @return pointer past the last entry; equals begin() when empty.
		 */
    const T* end() const noexcept
    {
      const std::size_t n = size();
      return n ? mData + n : mData;
    }

	private:
		/*!
		 * @brief Validate the element and byte counts before multiplying dimensions.
		 */
		static std::size_t CheckedSize(std::size_t nrow, std::size_t ncol)
		{
			// Division cannot wrap. Zero-sized shapes remain valid with a large other dimension.
			constexpr std::size_t maxCount = std::numeric_limits<std::size_t>::max()/sizeof(T);
			if( ncol && nrow > maxCount/ncol ) ERROR("Matrix dimensions exceed addressable storage.");
			return nrow*ncol;
		}

		std::size_t mNRow     = 0;       ///< Number of matrix rows.
		std::size_t mNCol     = 0;       ///< Number of matrix columns.
		std::size_t mCapacity = 0;       ///< Allocated entry count, independent of the current shape.
		T          *mData     = nullptr; ///< Contiguous row-major matrix entries.
};


//-----------------------------------------------------------------------------------


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
		explicit CWorkMatrixAS3(T *data, std::size_t *inuse, std::size_t *ninst, const std::size_t r, const std::size_t c) 
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


//-----------------------------------------------------------------------------------


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





