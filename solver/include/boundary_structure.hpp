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
		 */
		IBoundary(CConfig   *config_container,
				      CGeometry *geometry_container,
						  CMarker   *marker_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IBoundary(void);

		/*!
		 * @brief Function that returns a pointer to the current marker.
		 *
		 * @return pointer to mMarker
		 */
		CMarker *GetMarker(void) const {return mMarker;}

	protected:
		CMarker *mMarker;   ///< Pointer to the current marker container.

	private:
		// Disable default constructor.
		IBoundary(void) = delete;
		// Disable default copy constructor.
		IBoundary(const IBoundary&) = delete;
		// Disable default copy operator.
		IBoundary& operator=(IBoundary&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a boundary specification based on periodic BCs. 
 */
class CPeriodicBoundary : public IBoundary
{
	public:

		/*!
		 * @brief Constructor of CPeriodicBoundary, which initializes a periodic boundary condition.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] marker_container input marker container.
		 */
		CPeriodicBoundary(CConfig   *config_container,
				              CGeometry *geometry_container,
							        CMarker   *marker_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CPeriodicBoundary(void) override;

		/*!
		 * @brief Function that returns a pointer to the matching marker.
		 *
		 * @return pointer to mMatchingMarker
		 */
		CMarker *GetMatchingMarker(void) const {return mMatchingMarker;}

	protected:

	private:
		CMarker *mMatchingMarker;   ///< Pointer to the matching marker container.

		/*!
		 * @brief Function that finds the matching marker for this periodic boundary.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container. 
		 */
		void FindMatchingMarker(CConfig   *config_container,
				                    CGeometry *geometry_container);
};



