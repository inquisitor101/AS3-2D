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
		// Disable default constructor.
		CMarker(void) = delete;

		/*!
		 * @brief Constructor of CMarker, which is responsible for a single marker geometry.
		 *
		 * @param[in] iZone current zone index of this marker.
		 * @param[in] tagname name of the marker tag.
		 * @param[in] markertype type of the marker BC.
		 */
		CMarker(unsigned short iZone,
				    std::string    tagname,
						ETypeBCs       markertype);
	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~CMarker(void);


		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		unsigned short GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the value of mTypeMarker.
		 *
		 * @return mTypeMarker.
		 */
		ETypeBCs GetTypeMarker(void) const {return mTypeMarker;}

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


		/*!
		 * @brief Setter function which assigns the vector mElementIndices.
		 *
		 * @param[in] index vector containing the element indices on this marker.
		 */
		void SetElementIndices(as3vector1d<unsigned int> index) {mElementIndices = index;}

		/*!
		 * @brief Setter function which assigns the vector mElementFaces.
		 *
		 * @param[in] index vector containing the face indice of each element on this marker.
		 */
		void SetElementFaces(as3vector1d<EFaceElement> index) {mElementFaces = index;}



	protected:

	private:
		unsigned short            mZoneID;         ///< Current zone index.
		std::string               mNameMarkerTag;  ///< Name of the marker tag.
		ETypeBCs                  mTypeMarker;     ///< Type of marker region.
    as3vector1d<unsigned int> mElementIndices; ///< Element indices on this marker.
		as3vector1d<EFaceElement> mElementFaces;   ///< Face indices on each element of this marker.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class used for storing a (single) marker region.
 */
class CInterfaceMarker : public CMarker
{
	public:
		// Disable default constructor.
		CInterfaceMarker(void) = delete;

		/*!
		 * @brief Constructor of CInterfaceMarker, which is responsible for an interface marker.
		 *
		 * @param[in] iZone current zone index of this marker.
		 * @param[in] tagname name of the marker tag.
		 * @param[in] markertype type of the marker BC.
		 */
		explicit CInterfaceMarker(unsigned short iZone,
				                      std::string    tagnameI,
															unsigned short jZone,
															std::string    tagnameJ,
									            ETypeBCs       markertype);
	

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CInterfaceMarker(void) override;
	
		/*!
		 * @brief Getter function which returns the value of mMatchingZoneID.
		 *
		 * @return mMatchingZoneID.
		 */
		unsigned short GetMatchingZoneID(void) const {return mMatchingZoneID;}

		/*!
		 * @brief Getter function which returns the value of mMatchingNameMarkerTag.
		 *
		 * @return mMatchingNameMarkerTag.
		 */
		const std::string &GetMatchingNameMarkerTag(void) const {return mMatchingNameMarkerTag;}

		/*!
		 * @brief Getter function which returns mMatchingMarker.
		 *
		 * @return mMatchingMarker.
		 */
		std::weak_ptr<CMarker> GetMatchingMarker(void) const {return mMatchingMarker;}



		/*!
		 * @brief Setter function which assigns a value for the pointer mMatchingMarker.
		 *
		 * @param[in] mark weak pointer to the matching marker.
		 */
		void SetMatchingMarker(std::weak_ptr<CMarker> mark) {mMatchingMarker = mark;}


	
	protected:

	private:
		unsigned short         mMatchingZoneID;        ///< Matching interface zone index.
		std::string            mMatchingNameMarkerTag; ///< Name of the matching interface marker tag.
		std::weak_ptr<CMarker> mMatchingMarker;        ///< Matching interface marker.
		
};


