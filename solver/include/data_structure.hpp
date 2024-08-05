#pragma once

#include <cstring>
#include "error_structure.hpp"


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
		T*     mData = nullptr;
};


