#pragma once

/*!
 * @file
 * @brief Owning contiguous row-major matrix.
 */


#include <cstddef>
#include <cstring>
#include <limits>

#include "NError.hpp"
#include "common/AS3ScalarTypes.hpp"


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
