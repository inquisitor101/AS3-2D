#pragma once

/*!
 * @file
 * @brief AS3 binary grid metadata, marker records and coordinate import.
 */


#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <vector>

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
	 * @brief Function that determines whether swapping the read value reproduces the expected value.
	 *
	 * @tparam T matrix scalar type.
	 *
	 * @param[in] expect expected value in host byte order.
	 * @param[in] test value read from the file.
	 *
	 * @return true if byte swapping is required; terminates if neither byte order matches.
	 *
	 * @todo move this into NAS3BinaryFile namespace.
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


  /*!
   * @brief Reading and decoding the AS3 binary grid format.
   */

  namespace NAS3BinaryFile
  {

    /*!
     * @brief Unsigned 32-bit integer type stored in the AS3 binary file.
     */
    using as3fileuint = std::uint32_t;

    /*!
     * @brief Polynomial order and structured element counts for one AS3 grid zone.
     */
    struct CAS3ZoneMetadata
    {
      std::size_t mZoneID;  ///< Zero-based grid zone identifier.
      std::size_t mNPoly;   ///< Polynomial order of the geometry nodes.
      std::size_t mNiElem;  ///< Number of elements in the logical i-direction.
      std::size_t mNjElem;  ///< Number of elements in the logical j-direction.
    };

    /*!
     * @brief AS3 file identifiers mapped to solver face, marker and nodal-distribution enums.
     */
    struct CAS3Mappings
    {
      std::map<std::size_t, EFaceLocation>   mFaces;               ///< File face identifiers mapped to local edge locations.
      std::map<std::size_t, ETypeZoneMarker> mBoundaries;          ///< File marker identifiers mapped to internal or external marker types.
      std::map<std::size_t, ETypeDOF>        mNodalDistributions;  ///< File nodal-distribution identifiers mapped to solver enums.
    };

    /*!
     * @brief The AS3 grid header, byte order and zone metadata read from one file.
     */
    struct CAS3BinaryMetadata
    {
      std::string                   mFilename;          ///< Associated AS3 binary grid filename.
      std::size_t                   mNDim     = 0;      ///< Number of spatial dimensions recorded in the file.
      bool                          mByteSwap = false;  ///< Whether multibyte input values require byte swapping.
      std::vector<CAS3ZoneMetadata> mZones;             ///< Metadata indexed by zone.
      CAS3Mappings                  mMappings;          ///< File-to-solver enumeration mappings.
    };

    /*!
     * @brief An owning binary file handle together with its filename and byte-swap flag.
     */
    struct CBinaryFile
    {

      /*!
       * @brief Owning FILE pointer that closes its handle through std::fclose.
       */
      using uptr_file = std::unique_ptr<std::FILE, decltype(&std::fclose)>;
      uptr_file   mHandle{nullptr, &std::fclose};  ///< Owned FILE handle, closed automatically by its deleter.
      std::string mFilename;                       ///< Associated AS3 binary grid filename.
      bool        mSwap = false;                   ///< Whether values read through this handle require byte swapping.
    };


    /*!
     * @brief Function that reads and validates the AS3 binary grid header and enumeration mappings.
     *
     * @param[in] filename input filename.
     *
     * @return grid header, zone metadata and file-to-solver identifier mappings.
     */
    CAS3BinaryMetadata ReadAS3BinaryMetadata(const std::string &filename);

    /*!
     * @brief Function that imports one zone and its element coordinates using previously read metadata.
     *
     * @param[in] filename input filename.
     * @param[in] metadata previously read AS3 header and zone metadata.
     * @param[in] iZone zero-based zone index.
     *
     * @return owning pointer to the imported zone geometry.
     */
    std::unique_ptr<CSinglezoneGeometry> ReadAS3BinarySinglezoneGrid(const std::string        &filename,
                                                                     const CAS3BinaryMetadata &metadata,
                                                                     std::size_t               iZone);
 

    /*!
     * @brief Function that reads named physical boundary marker families from the AS3 file.
     *
     * @param[in] filename input filename.
     * @param[in] metadata previously read AS3 header and zone metadata.
     *
     * @return named physical boundary marker families.
     */
    as3vector1d<CExternalFamilyMarker> ReadAS3ExternalMarkers(const std::string        &filename,
                                                              const CAS3BinaryMetadata &metadata);

    /*!
     * @brief Function that reads internal paired-zone marker families from the AS3 file.
     *
     * @param[in] filename input filename.
     * @param[in] metadata previously read AS3 header and zone metadata.
     *
     * @return internal paired-zone marker families.
     */
    as3vector1d<CInternalFamilyMarker> ReadAS3InternalMarkers(const std::string        &filename,
                                                              const CAS3BinaryMetadata &metadata);



 

    /*!
     * @brief Function that reads element coordinates directly or converts them through a reusable precision
     *   buffer.
     *
     * @tparam InputPrecision floating-point precision stored in the coordinate payload.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     * @param[in,out] grid zone geometry whose element coordinates are filled.
     */
    template<class InputPrecision>
    void ReadElementCoordinates(const CBinaryFile   &file,
                                CSinglezoneGeometry &grid);

 

    /*!
     * @brief Function that opens an AS3 file, validates its magic number and determines byte order.
     *
     * @param[in] filename input filename.
     *
     * @return owning file handle with filename and byte-swap information.
     */
    CBinaryFile OpenAS3BinaryFile(const std::string &filename);

    /*!
     * @brief Function that identifies the coordinate precision from the remaining coordinate payload.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     * @param[in] nCoor number of coordinate scalars in the zone.
     *
     * @return number of bytes per coordinate scalar.
     */
    std::size_t DetermineCoordinatePrecisionBytes(const CBinaryFile &file,
                                             std::size_t             nCoor);

    /*!
     * @brief Function that reports a truncated or invalid AS3 input record and terminates execution.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     */
    void ReportReadError(const CBinaryFile &file);

    /*!
     * @brief Function that reads a 32-bit AS3 unsigned integer, swaps bytes if needed and widens it to
     *   std::size_t.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     *
     * @return decoded unsigned integer widened to std::size_t.
     */
    std::size_t ReadAS3UInt(const CBinaryFile &file);

    /*!
     * @brief Function that reads an AS3 length-prefixed string.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     *
     * @return decoded string.
     */
    std::string ReadString(const CBinaryFile &file);

    /*!
     * @brief Function that reads an AS3 boolean record.
     *
     * @param[in] file open AS3 binary file and its byte-order information.
     *
     * @return decoded boolean value.
     */
    bool ReadBoolean(const CBinaryFile &file);

  } // End of NAS3BinaryFile
}

#include "import_structure.inl"
