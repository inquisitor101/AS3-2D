#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "input_structure.hpp"
#include "marker_structure.hpp"



/*!
 * @brief A namespace used for storing specific file import utility functions.
 */
namespace NImportFile
{

	/*!
	 * @brief Function that imports an AS3 grid file.
	 */
	void ImportAS3Grid(CConfig            *config_container, 
			               CMultizoneGeometry *multizone_geometry_container);

	/*!
	 * @brief Function that imports an AS3 grid file in binary format.
	 */
	void ImportAS3GridBinary(CConfig            *config_container, 
			                     CMultizoneGeometry *multizone_geometry_container);


  // TODO: move this into NAS3BinaryFile namespace.
  /*!
	 * @brief Function that checks whether byte-swapping is required.
	 */
	template<typename T>
	bool CheckByteSwapping(const T expect, T test)
	{
		bool swap = ( expect != test ) ? true : false;

		// Check if byte-swapping works or not.
		if( swap )
		{
			T tmp = test;
			NInputUtility::SwapBytes(&tmp, sizeof(T), 1);
			if( expect != tmp )
				ERROR("Issue in file, could not read it.");
		}

		return swap;
	}


  // CHANGED: the below is based on the new grid format.

  namespace NAS3BinaryFile
  {
    using as3fileuint = std::uint32_t;

    struct CAS3ZoneMetadata
    {
      size_t mZoneID;
      size_t mNPoly;
      size_t mNiElem;
      size_t mNjElem;
    };

    struct CAS3Mappings
    {
      std::map<as3fileuint, EFaceLocation>   mFaces;
      std::map<as3fileuint, ETypeZoneMarker> mBoundaries;
      std::map<as3fileuint, ETypeDOF>        mNodalDistributions;
    };

    struct CAS3BinaryMetadata
    {
      std::string                   mFilename;
      std::size_t                   mNDim     = 0;
      bool                          mByteSwap = false;
      std::vector<CAS3ZoneMetadata> mZones;
      CAS3Mappings                  mMappings;
    };

    struct CBinaryFile
    {
      using uptr_file = std::unique_ptr<std::FILE, decltype(&std::fclose)>;
      uptr_file   mHandle{nullptr, &std::fclose};
      std::string mFilename;
      bool        mSwap = false;
    };


    CAS3BinaryMetadata ReadAS3BinaryMetadata(const std::string &filename);
    std::unique_ptr<CSinglezoneGeometry> ReadAS3BinarySinglezoneGrid(const std::string        &filename,
                                                                     const CAS3BinaryMetadata &metadata,
                                                                     std::size_t               iZone);
 
    as3vector1d<CExternalFamilyMarker> ReadAS3ExternalMarkers(const std::string        &filename,
                                                              const CAS3BinaryMetadata &metadata);

    as3vector1d<CInternalFamilyMarker> ReadAS3InternalMarkers(const std::string        &filename,
                                                              const CAS3BinaryMetadata &metadata);



 
    template<class InputPrecision>
    void ReadElementCoordinates(const CBinaryFile   &file,
                                CSinglezoneGeometry &grid);

 
    // Helper functions for reading entries.
    CBinaryFile OpenAS3BinaryFile(const std::string &filename);
    size_t DetermineCoordinatePrecisionBytes(const CBinaryFile &file,
                                             size_t             nCoor);

    void ReportReadError(const CBinaryFile &file);
    as3fileuint ReadAS3UInt(const CBinaryFile &file);
    std::string ReadString(const CBinaryFile &file);
    bool ReadBoolean(const CBinaryFile &file);

  } // End of NAS3BinaryFile
}

#include "import_structure.inl"
