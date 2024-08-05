#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "tensor_structure.hpp"


/*!
 * @brief A class used for initializing and definining a physical element. 
 */
class CPhysicalElement
{
	public:

		unsigned short        mNPoly;           ///< Solution polynomial order.
		CMatrixAS3<as3double> mSol2D;           ///< Solution on the physical element in 2D.
		CMatrixAS3<as3double> mInvMassMatrix;   ///< Inverse of the mass matrix.
		
		CMatrixAS3<as3double> mMetricSol2D;     ///< Metrics at the volume solution points in 2D. Rows: 
																				    ///< [0]: |J|, [1]: drdx, [2]: drdy, [3]: dsdx, [4]: dsdy.

		CMatrixAS3<as3double> mMetricInt2D;     ///< Metrics at the volume integration points in 2D. Rows: 
																				    ///< [0]: |J|, [1]: drdx, [2]: drdy, [3]: dsdx, [4]: dsdy.

		CMatrixAS3<as3double> mMetricIntIMin1D; ///< Metrics at the IMIN surface integration points in 1D. Rows:
																						///< [0]: ||n||, [1]: nx, [2]: ny, [3]: drdx, [4]: drdy, [5]: dsdx, [6]: dsdy.
		CMatrixAS3<as3double> mMetricIntIMax1D; ///< Metrics at the IMAX surface integration points in 1D. Rows:
																						///< [0]: ||n||, [1]: nx, [2]: ny, [3]: drdx, [4]: drdy, [5]: dsdx, [6]: dsdy.
		CMatrixAS3<as3double> mMetricIntJMin1D; ///< Metrics at the JMIN surface integration points in 1D. Rows:
																						///< [0]: ||n||, [1]: nx, [2]: ny, [3]: drdx, [4]: drdy, [5]: dsdx, [6]: dsdy.
		CMatrixAS3<as3double> mMetricIntJMax1D; ///< Metrics at the JMAX surface integration points in 1D. Rows:
																						///< [0]: ||n||, [1]: nx, [2]: ny, [3]: drdx, [4]: drdy, [5]: dsdx, [6]: dsdy.


		/*!
		 * @brief Default constructor of CPhysicalElement, which initializes a CPhysicalElement class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] element_geometry geometry of the physical element. 
		 * @param[in] iZone current zone ID.
		 * @param[in] nVar number of working variables.
		 */
		CPhysicalElement(CConfig          *config_container,
										 CStandardElement *standard_element,
										 ITensorProduct   *tensor_container,
										 CElementGeometry *element_geometry,
						         unsigned short    iZone,
										 unsigned short    nVar);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CPhysicalElement(void);


		/*!
		 * @brief Function that returns the number of working variables in the solution.
		 * 
		 * @return number of working variables.
		 */
		size_t GetnWorkingVar(void) const {return mSol2D.row();}

	protected:

	private:
	
		/*!
		 * @brief Function that computes the volume metrics at the solution points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsSolVolume(CStandardElement      *standard_element,
																 ITensorProduct        *tensor_container,
				                         CMatrixAS3<as3double> &coord);

		/*!
		 * @brief Function that computes the volume metrics at the integration points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsIntVolume(CStandardElement      *standard_element,
																 ITensorProduct        *tensor_container,
				                         CMatrixAS3<as3double> &coord);

		/*!
		 * @brief Function that computes the IMIN surface metrics at the integration points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsIntSurfIMIN(CStandardElement      *standard_element,
																   ITensorProduct        *tensor_container,
				                           CMatrixAS3<as3double> &coord);

		/*!
		 * @brief Function that computes the IMAX surface metrics at the integration points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsIntSurfIMAX(CStandardElement      *standard_element,
																   ITensorProduct        *tensor_container,
				                           CMatrixAS3<as3double> &coord);
	
		/*!
		 * @brief Function that computes the JMIN surface metrics at the integration points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsIntSurfJMIN(CStandardElement      *standard_element,
																   ITensorProduct        *tensor_container,
				                           CMatrixAS3<as3double> &coord);

		/*!
		 * @brief Function that computes the JMAX surface metrics at the integration points.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] coord physical coordinates over the element.
		 */
		void ComputeMetricsIntSurfJMAX(CStandardElement      *standard_element,
																   ITensorProduct        *tensor_container,
				                           CMatrixAS3<as3double> &coord);

	
		/*!
		 * @brief Function that computes the inverse mass matrix on the physical element.
		 * 
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 */
		void ComputeInverseMassMatrix(CStandardElement *standard_element,
																  ITensorProduct   *tensor_container);



		// Disable default constructor.
		CPhysicalElement(void) = delete;
		// Disable default copy constructor.
		CPhysicalElement(const CPhysicalElement&) = delete;
		// Disable default copy operator.
		CPhysicalElement& operator=(CPhysicalElement&) = delete;
};
