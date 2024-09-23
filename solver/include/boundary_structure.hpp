#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "marker_structure.hpp"


/*!
 * @brief An interface class used for the boundary condition.
 */
class IBoundary
{
	public:
	
		/*!
		 * @brief Constructor of IBoundary, which serves as an interface for the boundary condition.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] marker_container input marker container.
		 * @param[in] index index of the element.
		 */
		IBoundary(CConfig      *config_container,
				      CGeometry    *geometry_container,
						  CMarker      *marker_container,
							unsigned int  index);
	
		/*!
		 * @brief Function that returns the element index of this boundary.
		 */
		unsigned int GetIndex(void) const {return mIndex;}

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IBoundary(void);


		/*!
		 * @brief Pure virtual function that returns the type of BC. Must be overridden.
		 */
		virtual ETypeBC GetTypeBC(void) const = 0;

		/*!
		 * @brief Pure virtual function that computes a boundary state. Must be overridden.
		 */
		virtual void ComputeBoundaryState(CMatrixAS3<as3double>     &metric, 
				                              CWorkMatrixAS3<as3double> &varI,
																			CWorkMatrixAS3<as3double> &dVarDxI,
																			CWorkMatrixAS3<as3double> &dVarDyI,
																			CWorkMatrixAS3<as3double> &varJ,
																			CWorkMatrixAS3<as3double> &dVarDxJ,
																			CWorkMatrixAS3<as3double> &dVarDyJ) = 0;

	protected:
		unsigned int mIndex; ///< Index of the owner element.

	private:
		// Disable default constructor.
		IBoundary(void) = delete;
		// Disable default copy constructor.
		IBoundary(const IBoundary&) = delete;
		// Disable default copy operator.
		IBoundary& operator=(IBoundary&) = delete;
};

//-----------------------------------------------------------------------------------



