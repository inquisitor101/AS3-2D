#pragma once

#include <cstring>
#include "error_structure.hpp"


/*!
 * @brief A class for a matrix data structure. 
 */
template<typename T>
class CMatrixAS3
{
	public:
		
		// Default constructor.
		CMatrixAS3(void) = default;

		// Constructor, based on data size. Default is 1D array.
		explicit CMatrixAS3(size_t nrow, size_t ncol = 1) : mNRow(nrow), mNCol(ncol), mData(nullptr)
		{
			mData = new T[nrow*ncol];
			if( !mData ) ERROR("Allocation failed.");
			std::fill_n(mData, nrow*ncol, static_cast<T>(0.0));
		}

		// Copy constructor, does a deep copy of same class.
		CMatrixAS3(const CMatrixAS3<T> &other)
		{
			this->mNRow = other.mNRow;
			this->mNCol = other.mNCol;
			this->mData = new T[mNRow*mNCol];
			if( !mData ) ERROR("Allocation failed.");
			std::memcpy(mData, other.mData, mNRow*mNCol*sizeof(T));
		}

		// Move constructor.
		CMatrixAS3(CMatrixAS3&& other) noexcept
		{
			this->mNRow = other.mNRow;
			this->mNCol = other.mNCol;
			this->mData = other.mData;
			other.mNRow = 0;
			other.mNCol = 0;
			other.mData = nullptr;
		}


		// Destructor, which frees memory.
		~CMatrixAS3(void)
		{
			if( mData )
			{
				delete[] mData;
				mData = nullptr;	
			}
		}

		// Copy assignment operator.
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

		// Move assignment operator.
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


		// Resize the data and initializes all data to zero (always).
		void resize(size_t nrow, size_t ncol = 1)
		{
			const size_t size  = nrow*ncol;
			const size_t mSize = this->size(); 
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

		// Return the pointer to mData.
		T* data(void) 
		{
			return mData;
		}

		// Return a const pointer to mData.
		const T* data(void) const 
		{
			return mData;
		}

		// Return the total  size of the matrix.
		size_t size(void) const {return mNRow*mNCol;}
		// Return the row    size of the matrix.
		size_t row(void)  const {return mNRow;}
		// Return the column size of the matrix.
		size_t col(void)  const {return mNCol;}

		// Reset values to zero.
		void reset(void) const
		{
			for(size_t i=0; i<mNRow*mNCol; i++) mData[i] = static_cast<T>(0);
		}

		// Return value by 1D index [].
		T& operator[](const size_t index)
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		// Return constant value by 1D index [].
		const T& operator[](const size_t index) const
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		// Return value by 2D index ().
		T& operator()(const size_t r, const size_t c)
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}

		// Return constant value by 2D index ().
		const T& operator()(const size_t r, const size_t c) const
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}


	private:
		size_t mNRow = 0;
		size_t mNCol = 0;
		T     *mData = nullptr;
};


//-----------------------------------------------------------------------------------


/*!
 * @brief A class that shares existing contiguous memory, as a matrix. 
 */
template<typename T>
class CWorkMatrixAS3
{
	public:

		// Explicit constructor.
		explicit CWorkMatrixAS3(T *data, size_t *inuse, size_t *ninst, const size_t r, const size_t c) 
			: mData(data), mInUse(inuse), mNInst(ninst), mNRow(r), mNCol(c) { mID = *ninst; }

		// Default destructor.
		~CWorkMatrixAS3(void) 
		{
			// This check ensures that the last-in first-out (LIFO) memory model is respected.
			if( mID != *mNInst ) ERROR("Object creation/deletion needs to abide by the LIFO convention.");
			*mInUse -= mNRow*mNCol; (*mNInst)--;
		}

		// Return the pointer to mData.
		T* data(void) 
		{
			return mData;
		}

		// Return value by 1D index [].
		T& operator[](const size_t index)
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		// Return constant value by 1D index [].
		const T& operator[](const size_t index) const
		{
#if DEBUG
			if( index >= mNRow*mNCol ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[index];
		}

		// Return value by 2D index ().
		T& operator()(const size_t r, const size_t c)
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}

		// Return constant value by 2D index ().
		const T& operator()(const size_t r, const size_t c) const
		{
#if DEBUG
			if( (r >= mNRow) || (c >= mNCol) ) ERROR("Index exceeds maximum data size.");
#endif
			return mData[r*mNCol+c];
		}

	private:
		size_t  mID;
		size_t  mNRow  = 0;
		size_t  mNCol  = 0;	
		size_t *mInUse = nullptr;
		size_t *mNInst = nullptr;
		T      *mData  = nullptr;

		// Disable default constructor.
		CWorkMatrixAS3(void) = delete;
		// Disable default copy constructor.
		CWorkMatrixAS3(const CWorkMatrixAS3&) = delete;
		// Disable default copy operator.
		CWorkMatrixAS3& operator=(CWorkMatrixAS3&) = delete;
};


//-----------------------------------------------------------------------------------


/*!
 * @brief A class that creates one contiguous memory block, and issues chunks of this block as separate matrices. 
 */
template<typename T>
class CPoolMatrixAS3
{
	public:

		// Constructor, based on data size. Default is 1D array.
		explicit CPoolMatrixAS3(size_t nsize) 
			: mNSize(nsize), mInUse(new size_t(0)), mNInst(new size_t(0)), mData(nullptr)
		{
			mData = new T[nsize];
			if( !mData ) ERROR("Allocation failed.");
			std::fill_n(mData, nsize, static_cast<T>(0.0));
		}

		// Destructor, which frees memory.
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

		// Return a sepecific memory block, wrapped as a matrix.
		CWorkMatrixAS3<T> GetWorkMatrixAS3(const size_t r, const size_t c) 
		{
			const size_t n  =  r*c;    // overall size.
			const size_t I0 = *mInUse; // address of starting index.
			const size_t I1 =  I0 + n; // address of last index.
			
			if( I1 > mNSize ) ERROR("Requested size exceeds allocated memory.");
			*mInUse = I1;              // Update used memory so far.
			(*mNInst)++;               // Update the number of instances in use.

			// Issue (via RVO) a work matrix, based on the specified dimensions.
			return CWorkMatrixAS3<T>( &mData[I0], mInUse, mNInst, r, c ); 
		}

		// Return the total size of the memory pool.
		size_t size(void)      const {return  mNSize;}
		// Return the number of object instantances in use.
		size_t ninst(void)     const {return *mNInst;}
		// Return the currently used memory of the memory pool.
		size_t inuse(void)     const {return *mInUse;}
		// Return the available memory for usage.
		size_t available(void) const {return (mNSize - *mInUse);}


	private:
		size_t  mNSize = 0;
		size_t *mInUse = nullptr;
		size_t *mNInst = nullptr;
		T      *mData  = nullptr;

		// Disable default constructor.
		CPoolMatrixAS3(void) = delete;
		// Disable default copy constructor.
		CPoolMatrixAS3(const CPoolMatrixAS3&) = delete;
		// Disable default copy operator.
		CPoolMatrixAS3& operator=(CPoolMatrixAS3&) = delete;
};





