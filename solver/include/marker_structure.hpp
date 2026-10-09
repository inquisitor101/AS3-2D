#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"

// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CExternalFamilyMarker;


// NOTE, markers are divided into three:
//   *) external, which are boundary markers.
//   *) internal, which are interface markers.
//   *) periodic, which.. are periodic markers.
// ... all are defined on each zone's entire edge.


// Marker for a family of periodic zone-boundaries.
class CPeriodicFamilyMarker
{
  private:
    
    // Marker for a periodic zone-boundary.
    struct CPeriodicMarker
    {
      size_t mIndexZoneI;
      size_t mIndexZoneJ;

      EFaceLocation mFaceLocationI;
      EFaceLocation mFaceLocationJ;

      std::string mNameI;
      std::string mNameJ;
      
      bool mIsReversed;

      std::array<as3double, 2> mTranslationVector; // From I to J.
    };

  public:
    
    CPeriodicFamilyMarker(const CMultizoneGeometry       *multizone_geometry_container,
                          const CExternalFamilyMarker    &iexternal_marker,
                          const CExternalFamilyMarker    &jexternal_marker,
                          const std::array<as3double, 2> &translation);

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CPeriodicMarker> mMarkers;
};



// Marker for a family of internal zone-boundaries (i.e. interfaces).
class CInternalFamilyMarker
{
  private:
    
    // Marker for an internal zone-boundary (i.e. interface).
    struct CInternalMarker
    {
      size_t mIndexZoneI;
      size_t mIndexZoneJ;

      EFaceLocation mFaceLocationI;
      EFaceLocation mFaceLocationJ;

      bool mIsReversed;
    };

  public:
    
    CInternalFamilyMarker(std::size_t nMarker)
    {
      mMarkers.reserve( nMarker );
    }

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CInternalMarker> mMarkers;
};


// Marker for a family of external zone-boundaries.
class CExternalFamilyMarker
{
  private:

    // Marker for an external zone-boundary.
    struct CExternalMarker
    {
      std::size_t   mIndexZone;
      EFaceLocation mFaceLocation;
    };

  public:
    
    CExternalFamilyMarker(std::string name,
                          std::size_t nMarker)
      : mName(std::move(name))
    {
      mMarkers.reserve( nMarker );
    }

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const std::string& GetName(void) const {return mName;}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}

  private:
    
    std::string mName;
    as3vector1d<CExternalMarker> mMarkers;
};




/*!
 * @brief A class used for storing a single element face in a marker region.
 */
struct CFaceMarker
{
	unsigned int  mIndex; ///< Element index containing this face marker.
	EFaceLocation mFace;  ///< Face location for this marker.
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
		CMarker(unsigned short             zone,
				    ETypeBC                    type,
				    std::string                name,	
						as3vector1d<EFaceLocation> face,
						as3vector1d<unsigned int>  mark);
	
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
		unsigned int GetnElem(void) const {return static_cast<unsigned int>( mElementFaces.size() );}

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

    EFaceLocation GetFaceLocation(void) const {return mFaceLocation;}
    ETypeFace     GetTypeFace(void)     const {return mFaceType;}

	protected:

	private:
		unsigned short           mZoneID;       ///< Current zone index.
		ETypeBC                  mTypeBC;       ///< Type of boundary condition.
		std::string              mNameMarker;   ///< Name of the marker tag.
		EFaceLocation            mFaceLocation; 
    ETypeFace                mFaceType;
    as3vector1d<CFaceMarker> mElementFaces; ///< Elements and their faces on this marker.

		// Disable default constructor.
		CMarker(void) = delete;
		// Disable default copy constructor.
		CMarker(const CMarker&) = delete;
		// Disable default copy operator.
		CMarker& operator=(CMarker&) = delete;
};


