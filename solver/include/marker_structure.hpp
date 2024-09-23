#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"


/*!
 * @brief A class used for storing a single element face in a marker region.
 */
struct CFaceMarker
{
	unsigned int mIndex; ///< Element index containing this face marker.
	EFaceElement mFace;  ///< Face location for this marker.
};


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
		 * @param[in] type type of BC.
		 * @param[in] tagname name of the marker tag.
		 * @param[in] face face locations on each element of this marker.
		 * @param[in] element element indices on this marker. 
		 */
		CMarker(unsigned short            zone,
				    ETypeBC                   type,
				    std::string               name,	
						as3vector1d<EFaceElement> face,
						as3vector1d<unsigned int> mark);
	
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
		 * @brief Getter function which returns the value of mTypeBC.
		 *
		 * @return mTypeBC.
		 */
		ETypeBC GetTypeBC(void) const {return mTypeBC;}

		/*!
		 * @brief Getter function which returns the value of mNameMarker.
		 *
		 * @return mNameMarker.
		 */
		const std::string &GetNameMarker(void) const {return mNameMarker;}

		/*!
		 * @brief Getter function which returns the number of elements.
		 *
		 * @return mElementFaces.size().
		 */
		size_t GetnElem(void) const {return mElementFaces.size();}

		/*!
		 * @brief Getter function which returns mElementFaces.
		 *
		 * @return mElementFaces.
		 */
		const as3vector1d<CFaceMarker> &GetElementFaces(void) const {return mElementFaces;}

		/*!
		 * @brief Getter function which returns value of mElementFaces at a specific index.
		 *
		 * @return mElementFaces[index].
		 */
		const CFaceMarker &GetElementFaces(size_t index) const {return mElementFaces[index];}

	protected:

	private:
		unsigned short           mZoneID;       ///< Current zone index.
		ETypeBC                  mTypeBC;       ///< Type of boundary condition.
		std::string              mNameMarker;   ///< Name of the marker tag.
		as3vector1d<CFaceMarker> mElementFaces; ///< Elements and their faces on this marker.

		// Disable default constructor.
		CMarker(void) = delete;
		// Disable default copy constructor.
		CMarker(const CMarker&) = delete;
		// Disable default copy operator.
		CMarker& operator=(CMarker&) = delete;
};


