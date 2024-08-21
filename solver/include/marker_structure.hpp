#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"


/*!
 * @brief A class used for storing a (single) marker region.
 */
class CMarker
{
	public:

		/*!
		 * @brief Constructor of CMarker, which is responsible for a single marker geometry.
		 *
		 * @param[in] iZone current zone index of this marker.
		 * @param[in] tagname name of the marker tag.
		 */
		CMarker(unsigned short            iZone,
				    std::string               tagname,
						ETypeBC                   typemarker,
						as3vector1d<unsigned int> elements,
						as3vector1d<EFaceElement> faces);
	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CMarker(void);


		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		unsigned short GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the value of mMarkerBC.
		 *
		 * @return mMarkerBC.
		 */
		ETypeBC GetMarkerBC(void) const {return mMarkerBC;}

		/*!
		 * @brief Getter function which returns the value of mNameMarkerTag.
		 *
		 * @return mNameMarkerTag.
		 */
		const std::string &GetNameMarkerTag(void) const {return mNameMarkerTag;}

		/*!
		 * @brief Getter function which returns the value of mElementIndices.
		 *
		 * @return mElementIndices.
		 */
		const as3vector1d<unsigned int> &GetElementIndices(void) const {return mElementIndices;}

		/*!
		 * @brief Getter function which returns the value of mElementFaces.
		 *
		 * @return mElementFaces.
		 */
		const as3vector1d<EFaceElement> &GetElementFaces(void) const {return mElementFaces;}

	protected:

	private:
		unsigned short            mZoneID;         ///< Current zone index.
		std::string               mNameMarkerTag;  ///< Name of the marker tag.
		ETypeBC                   mMarkerBC;       ///< Type of marker region.
    as3vector1d<unsigned int> mElementIndices; ///< Element indices on this marker.
		as3vector1d<EFaceElement> mElementFaces;   ///< Face indices on each element of this marker.


		// Disable default constructor.
		CMarker(void) = delete;
		// Disable default copy constructor.
		CMarker(const CMarker&) = delete;
		// Disable default copy operator.
		CMarker& operator=(CMarker&) = delete;
};


