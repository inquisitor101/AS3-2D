#pragma once 

#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


struct CIndexElement
{
	CIndexElement(unsigned short iZone, unsigned int iElem) : mZone(iZone), mElem(iElem) {}
	unsigned short mZone;
	unsigned int   mElem;
};


/*!
 * @brief A class for the OpenMP shared memory parallelization. 
 */
class COpenMP
{
	public:
		
		/*!
		 * @brief Constructor of COpenMP, which initializes the OpenMP class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 */
		COpenMP(CConfig                               *config_container,
				    CGeometry                             *geometry_container,
						as3vector1d<std::unique_ptr<ISolver>> &solver_container);
	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~COpenMP(void);

		/*!
		 * @brief Getter function which returns the total number of volume elements.
		 */
		size_t GetnElemTotal(void) const {return mIndexVolume.size();}

		/*!
		 * @brief Getter function which returns the total number of internal IMAX faces.
		 */
		size_t GetnInternIFace(void) const {return mInternIFace.size();}

		/*!
		 * @brief Getter function which returns the total number of internal JMAX faces.
		 */
		size_t GetnInternJFace(void) const {return mInternJFace.size();}

		/*!
		 * @brief Getter function that returns the value of mIndexVolume at the specified index.
		 *
		 * @param[in] index index of the vector.
		 *
		 * @return mIndexVolume[index]
		 */
		CIndexElement *GetIndexVolume(size_t index) const {return mIndexVolume[index].get();}

		/*!
		 * @brief Getter function that returns the value of mInternIFace at the specified index.
		 *
		 * @param[in] index index of the vector.
		 *
		 * @return mInternIFace[index]
		 */
		CIndexElement *GetInternIFace(size_t index) const {return mInternIFace[index].get();}

		/*!
		 * @brief Getter function that returns the value of mInternJFace at the specified index.
		 *
		 * @param[in] index index of the vector.
		 *
		 * @return mInternJFace[index]
		 */
		CIndexElement *GetInternJFace(size_t index) const {return mInternJFace[index].get();}

		/*!
		 * @brief Getter function that returns the IMin residual, based on the specified index.
		 *
		 * @return mResIMin[index].
		 */
		CMatrixAS3<as3double> &GetResIMin(size_t index) {return mResIMin[index];}

		/*!
		 * @brief Getter function that returns the JMin residual, based on the specified index.
		 *
		 * @return mResJMin[index].
		 */
		CMatrixAS3<as3double> &GetResJMin(size_t index) {return mResJMin[index];}

	protected:

	private:
		as3vector1d<std::unique_ptr<CIndexElement>> mIndexVolume; ///< Global map for the volume indices.
		as3vector1d<std::unique_ptr<CIndexElement>> mInternIFace; ///< Global map for the (IMAX) indices.
		as3vector1d<std::unique_ptr<CIndexElement>> mInternJFace; ///< Global map for the (JMAX) indices.

		as3vector1d<CMatrixAS3<as3double>> mResIMin;
		as3vector1d<CMatrixAS3<as3double>> mResJMin;


		/*!
		 * @brief Function that initializes the volume indices for each element.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitializeVolumeIndices(CConfig   *config_container,
				                         CGeometry *geometry_container);

		/*!
		 * @brief Function that initializes the internal face indices in the i-direction.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitializeSurfaceIFaces(CConfig   *config_container,
				                         CGeometry *geometry_container);

		/*!
		 * @brief Function that initializes the internal face indices in the j-direction.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitializeSurfaceJFaces(CConfig   *config_container,
				                         CGeometry *geometry_container);

		// Disable default constructor.
		COpenMP(void) = delete;
		// Disable default copy constructor.
		COpenMP(const COpenMP&) = delete;
		// Disable default copy operator.
		COpenMP& operator=(COpenMP&) = delete;
};

