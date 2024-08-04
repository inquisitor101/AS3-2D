#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"
#include "marker_structure.hpp"

// Forward declaration to avoid compiler problems.
class CZoneGeometry;
class CElementGeometry;


/*!
 * @brief A class used for storing the entire (multi-zone) grid.
 */
class CGeometry
{
	public:
		// Disable default constructor.
		CGeometry(void) = delete;

		/*!
		 * @brief Constructor of CGeometry, which is responsible for the entire grid geometry.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 */
		CGeometry(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CGeometry(void);


		/*!
		 * @brief Getter function which returns a vector of grid zone containers.
		 *
		 * @return mZoneGeometry
		 */
		const as3vector1d<std::unique_ptr<CZoneGeometry>> &GetZoneGeometry(void) const {return mZoneGeometry;}

		/*!
		 * @brief Getter function which returns a specific grid zone container.
		 *
		 * @return mZoneGeometry[iZone]
		 */
		CZoneGeometry *GetZoneGeometry(unsigned short iZone) const {return mZoneGeometry[iZone].get();}

	protected:

	private:
		const unsigned short                        mNZone;            ///< Total number of zones.
		as3vector1d<std::unique_ptr<CZoneGeometry>> mZoneGeometry;     ///< Container with the zone geometry.
    as3vector1d<std::string>                    mZoneGridFilename; ///< Grid filename per each zone.

		/*!
		 * @brief Function that extracts the grid zone connectivity information.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 */
    void ExtractZoneConnectivity(CConfig *config_container);

		/*!
		 * @brief Function that maps each interface to its assigned neighbor.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 */
    void MatchInterfaceMarkers(CConfig *config_container);
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class used for storing a (single) zone grid.
 */
class CZoneGeometry
{
	public:
		// Disable default constructor.
		CZoneGeometry(void) = delete;

		/*!
		 * @brief Constructor of CZoneGeometry, which is responsible for a single grid geometry.
		 * 
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] gridfiles entire multi-zone grid filenames.
		 * @param[in] iZone current zone index.
		 */
		CZoneGeometry(CConfig        *config_container,
									std::string     gridfile,
									unsigned short  iZone);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CZoneGeometry(void);



		/*!
		 * @brief Function that initializes and defines all the element geometry in this zone.
		 */
		void InitializeElements(as3vector2d<double> &x, 
				                    as3vector2d<double> &y,
														unsigned int         nxElem,
														unsigned int         nyElem);

		/*!
		 * @brief Function that defines all the interface markers in this zone.
		 */
		void InitializeMarkers(as3vector2d<unsigned int> &mark,
				                   as3vector2d<EFaceElement> &face,
				                   as3vector1d<std::string>  &name);


		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		unsigned short GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the value of mGridFile.
		 *
		 * @return mGridFile.
		 */
		const std::string &GetGridFile(void) const {return mGridFile;}

		/*!
		 * @brief Getter function which returns the mMarkerInterface.
		 *
		 * @return mMarkerInterface.
		 */
		as3vector1d<std::shared_ptr<CInterfaceMarker>> &GetMarkerInterface(void) {return mMarkerInterface;}

		/*!
		 * @brief Getter function which returns the element of the specified index.
		 *
		 * @return mElementGeometry[index].
		 */
		CElementGeometry *GetElementGeometry(size_t index) const {return mElementGeometry[index].get();}

		/*!
		 * @brief Getter function which returns the total number of elements in this zone.
		 *
		 * @return mElementGeometry.size().
		 */
		size_t GetnElem(void) const {return mElementGeometry.size();}

		/*!
		 * @brief Getter function which returns the number of elements in the x-direction.
		 *
		 * @return mNxElem.
		 */
		size_t GetnxElem(void) const {return static_cast<size_t>(mNxElem);}

		/*!
		 * @brief Getter function which returns the number of elements in the y-direction.
		 *
		 * @return mNyElem.
		 */
		size_t GetnyElem(void) const {return static_cast<size_t>(mNyElem);}

		/*!
		 * @brief Getter function which returns the number of grid nodes in 2D.
		 *
		 * @return (mNPolyGrid+1)*(mNPolyGrid+1).
		 */
		unsigned int GetnNodeGrid2D(void) const {return static_cast<unsigned int>( (mNPolyGrid+1)*(mNPolyGrid+1) );}

		/*!
		 * @brief Getter function which returns the polynomial order of the grid element.
		 *
		 * @return mNPolyGrid.
		 */
		unsigned short GetnPolyGrid(void) const {return mNPolyGrid;}


	protected:

	private:
		const unsigned short                           mZoneID;          ///< Current zone index.
		const std::string															 mGridFile;				 ///< Associated grid file name of this zone.
		unsigned short                                 mNPolyGrid;       ///< Polynomial order of the grid element.
		unsigned int                                   mNxElem;          ///< Number of elements in x-direction.
		unsigned int                                   mNyElem;          ///< Number of elements in y-direction.	
		as3vector1d<std::shared_ptr<CInterfaceMarker>> mMarkerInterface; ///< Container with the marker data.
		as3vector1d<std::unique_ptr<CElementGeometry>> mElementGeometry; ///< Container with the element geometry.

		/*!
		 * @brief Function that reads the marker interface for this zone.
		 * 
		 * @param[in] config_container pointer to the configuration container.
		 */
		void ReadMarkerInterfaceZone(CConfig *config_container);
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class used for storing a (single) zone grid.
 */
class CElementGeometry
{
	public:
		// Disable default constructor.
		CElementGeometry(void) = delete;

		/*!
		 * @brief Constructor of CElementGeometry, which is responsible for a single element geometry.
		 */
		CElementGeometry(as3vector1d<double> &x, as3vector1d<double> &y);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CElementGeometry(void);
	

		/*!
		 * @brief Getter function which returns the coordinates at the solution DOFs.
		 *
		 * @return mCoordSolDOFs.
		 */
		CMatrixAS3<as3double> &GetCoordSolDOFs(void) {return mCoordSolDOFs;}


	protected:

	private:
		CMatrixAS3<as3double> mCoordSolDOFs;  ///< Coordinates at the solution DOFs.
};


