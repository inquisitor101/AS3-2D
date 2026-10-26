#pragma once

/*!
 * @file
 * @brief Contiguous owning matrices and LIFO workspace matrices.
 */


#include <algorithm>
#include <cstddef>

#include <cstring>
#include "error_structure.hpp"


/*!
 * @brief An owning row-major matrix with contiguous entries and zero-initialized allocations.
 *
 * @tparam T type of each stored scalar.
 */
template<typename T>
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
		explicit CMatrixAS3(std::size_t nrow, std::size_t ncol = 1) : mNRow(nrow), mNCol(ncol), mData(nullptr)
		{
			mData = new T[nrow*ncol];
			if( !mData ) ERROR("Allocation failed.");
			std::fill_n(mData, nrow*ncol, static_cast<T>(0.0));
		}

		/*!
		 * @brief Copy constructor, does a deep copy of same class.
		 *
		 * @param[in] other source matrix or object.
		 */
		CMatrixAS3(const CMatrixAS3<T> &other)
		{
			this->mNRow = other.mNRow;
			this->mNCol = other.mNCol;
			this->mData = new T[mNRow*mNCol];
			if( !mData ) ERROR("Allocation failed.");
			std::memcpy(mData, other.mData, mNRow*mNCol*sizeof(T));
		}

		/*!
		 * @brief Move constructor.
		 *
		 * @param[in,out] other source matrix or object.
		 */
		CMatrixAS3(CMatrixAS3&& other) noexcept
		{
			this->mNRow = other.mNRow;
			this->mNCol = other.mNCol;
			this->mData = other.mData;
			other.mNRow = 0;
			other.mNCol = 0;
			other.mData = nullptr;
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
				if( (mNRow != other.mNRow) || (mNCol != other.mNCol) )
				{
					if( mData ) delete[] mData;
					mData = new T[other.size()];
					if( !mData ) ERROR("Allocation failed.");
					this->mNRow = other.mNRow;
					this->mNCol = other.mNCol;
				}
				std::memcpy(mData, other.mData, other.size()*sizeof(T));
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
		CMatrixAS3<T>& operator=(CMatrixAS3<T> &&other)
		{
			if( this != &other )
			{
				if( mData ) delete[] mData;
				this->mData = other.mData;
				this->mNRow = other.mNRow;
				this->mNCol = other.mNCol;
				other.mData = nullptr;
				other.mNRow = 0;
				other.mNCol = 0;
			}
			return *this;
		}


		/*!
		 * @brief Function that resizes the matrix and resets every entry to zero.
		 *
		 * @param[in] nrow number of rows.
		 * @param[in] ncol number of columns; defaults to one.
		 */
		void resize(std::size_t nrow, std::size_t ncol = 1)
		{
			const std::size_t size  = nrow*ncol;
			const std::size_t mSize = this->size(); 
			if( size <= mSize )
			{
				this->mNRow = nrow;
				this->mNCol = ncol;
			}
			else
			{
				T* data = new T[size];
				if( !data ) ERROR("Resize failed.");
				delete[] mData;
				this->mData = data;
				this->mNRow = nrow;
				this->mNCol = ncol;
			}
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
		void reset(void) const
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


	private:
		std::size_t mNRow = 0;   ///< Number of matrix rows.
		std::size_t mNCol = 0;   ///< Number of matrix columns.
		T     *mData = nullptr;  ///< Contiguous row-major matrix entries.
};


//-----------------------------------------------------------------------------------


/*!
 * @brief A non-owning row-major matrix borrowed from a LIFO workspace pool.
 *
 * @tparam T type of each stored scalar.
 *
 * @note Workspace matrices must be destroyed in reverse creation order. The backing pool must remain alive.
 */
template<typename T>
class CWorkMatrixAS3
{
	public:

		/*!
		 * @brief Constructor that wraps a borrowed pool block and records its LIFO identity.
		 *
		 * @param[in] data borrowed contiguous matrix storage; the pool retains ownership.
		 * @param[in] inuse shared pool entry counter.
		 * @param[in] ninst shared live-workspace counter.
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
		T      *mData  = nullptr;       ///< Borrowed matrix entries; ownership remains with the pool.
		std::size_t *mInUse = nullptr;  ///< Shared count of scalar entries currently borrowed from the pool.
		std::size_t *mNInst = nullptr;  ///< Shared count of live workspace matrices.
		std::size_t  mNRow  = 0;        ///< Number of matrix rows.
		std::size_t  mNCol  = 0;        ///< Number of matrix columns.
		std::size_t  mID;               ///< LIFO instance identifier assigned when this matrix was borrowed.
		
		
		

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
 * @tparam T type of each stored scalar.
 *
 * @note All borrowed matrices must be released before the pool is destroyed. Each OpenMP thread owns its own
 *   pool.
 */
template<typename T>
class CPoolMatrixAS3
{
	public:

		/*!
		 * @brief Constructor that allocates and zero-initializes the workspace pool.
		 *
		 * @param[in] nsize total number of scalar entries to allocate.
		 */
		explicit CPoolMatrixAS3(std::size_t nsize) 
			: mNSize(nsize), mInUse(new std::size_t(0)), mNInst(new std::size_t(0)), mData(nullptr)
		{
			mData = new T[nsize];
			if( !mData ) ERROR("Allocation failed.");
			std::fill_n(mData, nsize, static_cast<T>(0.0));
		}

		/*!
		 * @brief Destructor that releases the pool after all borrowed matrices have been destroyed.
		 */
		~CPoolMatrixAS3(void)
		{
			if( mData )
			{
				if( *mNInst )            ERROR("Some shared pointer instances still exist.");
				if( !mInUse || !mNInst ) ERROR("Unexcpected behavior encountered.");	
				delete[] mData;  mData  = nullptr;
				delete   mInUse; mInUse = nullptr;
				delete   mNInst; mNInst = nullptr;
			}
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
			const std::size_t n  =  r*c;    // overall size.
			const std::size_t I0 = *mInUse; // address of starting index.
			const std::size_t I1 =  I0 + n; // address of last index.
			
			if( I1 > mNSize ) ERROR("Requested size exceeds allocated memory.");
			*mInUse = I1;              // Update used memory so far.
			(*mNInst)++;               // Update the number of instances in use.

			// Issue (via RVO) a work matrix, based on the specified dimensions.
			return CWorkMatrixAS3<T>( &mData[I0], mInUse, mNInst, r, c ); 
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
		std::size_t ninst(void)     const {return *mNInst;}

		/*!
		 * @brief Return the currently used memory of the memory pool.
		 *
		 * @return number of scalar entries currently borrowed.
		 */
		std::size_t inuse(void)     const {return *mInUse;}

		/*!
		 * @brief Return the available memory for usage.
		 *
		 * @return number of scalar entries available for borrowing.
		 */
		std::size_t available(void) const {return (mNSize - *mInUse);}


	private:
		std::size_t  mNSize = 0;        ///< Total number of scalar entries allocated in the pool.
		std::size_t *mInUse = nullptr;  ///< Shared count of scalar entries currently borrowed from the pool.
		std::size_t *mNInst = nullptr;  ///< Shared count of live workspace matrices.
		T      *mData  = nullptr;       ///< Owned contiguous pool storage; separate from the shared usage counters.

		// Disable default constructor.
		CPoolMatrixAS3(void) = delete;
		// Disable default copy constructor.
		CPoolMatrixAS3(const CPoolMatrixAS3&) = delete;
		// Disable default copy operator.
		CPoolMatrixAS3& operator=(CPoolMatrixAS3&) = delete;
};





