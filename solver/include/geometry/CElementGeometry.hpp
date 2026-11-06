#pragma once

/*!
 * @file
 * @brief Coordinates of one quadrilateral grid element.
 */


#include <cstddef>

#include "common/AS3Options.hpp"
#include "data/CMatrixAS3.hpp"


/*!
 * @brief Coordinates of the geometry nodes of one quadrilateral element.
 */
class CElementGeometry
{
	public:

    /*!
     * @brief Constructor that initializes CElementGeometry from the specified data.
     *
     * @param[in] nNode2D number of geometry nodes in the element.
     */
    explicit CElementGeometry(std::size_t nNode2D);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CElementGeometry(void);
	

		/*!
		 * @brief Getter function which returns the coordinates at the solution DOFs.
		 *
		 * @return mCoordSolDOFs.
		 */
		CMatrixAS3<as3double> &GetCoordSolDOFs(void) {return mCoordSolDOFs;}

		/*!
		 * @brief Function that compares matching coordinates across a translated interface.
		 *
		 * @param[in] xi coordinate on the side to be translated.
		 * @param[in] xj matching coordinate on the other side.
		 * @param[in] translation displacement from the i side to the j side.
		 * @return whether the coordinates agree within the geometric tolerance and roundoff.
		 */
		static bool CoordinatesMatch(as3double xi,
		                             as3double xj,
		                             as3double translation);

	protected:

	private:
		CMatrixAS3<as3double> mCoordSolDOFs;  ///< Coordinates at the solution DOFs.


		// Disable default constructor.
		CElementGeometry(void) = delete;
		// Disable default copy constructor.
		CElementGeometry(const CElementGeometry&) = delete;
		// Disable default copy operator.
		CElementGeometry& operator=(CElementGeometry&) = delete;	
};
