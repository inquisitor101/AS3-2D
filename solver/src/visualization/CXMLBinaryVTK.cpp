#include "visualization/CXMLBinaryVTK.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string_view>
#include <type_traits>

#include "CConfig.hpp"
#include "geometry/CElementGeometry.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"


namespace
{
	/*!
	 * @brief Function that computes the byte count stored in an appended array header.
	 */
	std::uint64_t CountBytes(std::size_t nValue, std::size_t valueSize)
	{
		if( nValue > std::numeric_limits<std::uint64_t>::max()/valueSize )
			ERROR("VTK XML array exceeds the 64-bit byte count range.");
		return static_cast<std::uint64_t>(nValue)*valueSize;
	}

	/*!
	 * @brief Function that writes a contiguous range in native byte order.
	 */
	template<typename T>
	void WriteValues(std::FILE *file, const T *values, std::size_t nValue)
	{
		if( nValue && std::fwrite(values, sizeof(T), nValue, file) != nValue )
			ERROR("Could not write VTK XML visualization data.");
	}

	/*!
	 * @brief Function that writes the byte count and payload of one appended array.
	 */
	template<typename T>
	void WriteArray(std::FILE *file, const as3vector1d<T> &values)
	{
		const std::uint64_t nByte = CountBytes(values.size(), sizeof(T));
		WriteValues(file, &nByte, 1);
		WriteValues(file, values.data(), values.size());
	}
}


//-----------------------------------------------------------------------------------
// CXMLBinaryVTK member functions.
//-----------------------------------------------------------------------------------


CXMLBinaryVTK::CXMLBinaryVTK
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
	:
		IFileVTK(config_container, multizone_geometry_container)
 /*
	* Constructor that initializes the VTK XML writer and selects its buffer precision.
	*/
{
	if( config_container->GetOutputVisPrecision() == EFloatingPointPrecision::DOUBLE )
		mRealData.emplace<CRealData<double>>();
}

//-----------------------------------------------------------------------------------

template<typename T>
void CXMLBinaryVTK::InitializeMesh
(
 const CMultizoneGeometry *multizone_geometry_container,
 CRealData<T>             &realData
)
 /*
	* Function that caches the grid as linear quadrilaterals with one point per solution DOF.
	*/
{
	constexpr std::size_t  nCorner  = 4;
	constexpr std::uint8_t quadType = 9;
	constexpr std::size_t  maxIndex = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

	const std::size_t nZone = multizone_geometry_container->GetnZone();
	mZonePointCounts.resize(nZone);
	realData.mPointData.resize(nZone);

	// Count the mesh entries before allocating. Cell offsets count individual point indices.
	std::size_t nPointTotal = 0;
	std::size_t nCellTotal  = 0;
	for(std::size_t iZone=0; iZone<nZone; iZone++)
	{
		const auto* zone = multizone_geometry_container->GetSinglezoneGeometry(iZone);
		const std::size_t nElem = zone->GetnElem();
		const std::size_t nNode = zone->GetnNodeGrid2D();
		const std::size_t nPoly = zone->GetnPolyGrid();

		if( (nElem && nNode > (maxIndex-nPointTotal)/nElem) ||
		    (nPoly && nPoly > maxIndex/nPoly) )
			ERROR("Mesh exceeds the VTK XML 32-bit point or cell index range.");

		const std::size_t nCellElem = nPoly*nPoly;
		if( nElem && nCellElem > (maxIndex/nCorner-nCellTotal)/nElem )
			ERROR("Mesh exceeds the VTK XML 32-bit cell offset range.");

		const std::size_t nPoint = nElem*nNode;
		if( mNPointComponents && nPoint > realData.mPointData[iZone].max_size()/mNPointComponents )
			ERROR("VTK XML point data exceeds the available buffer size.");

		mZonePointCounts[iZone] = nPoint;
		nPointTotal += nPoint;
		nCellTotal  += nElem*nCellElem;
	}
	if( nPointTotal > realData.mPoints.max_size()/3 )
		ERROR("VTK XML coordinates exceed the available buffer size.");

	realData.mPoints.resize(3*nPointTotal);
	mConnectivity.resize(nCorner*nCellTotal);
	mOffsets.resize(nCellTotal);
	mCellTypes.resize(nCellTotal, quadType);

	std::size_t iPoint = 0;
	std::size_t iCell  = 0;
	for(std::size_t iZone=0; iZone<nZone; iZone++)
	{
		const auto* zone = multizone_geometry_container->GetSinglezoneGeometry(iZone);
		const std::size_t nNode   = zone->GetnNodeGrid2D();
		const std::size_t nPoly   = zone->GetnPolyGrid();
		const std::size_t nNode1D = nPoly + 1;

		realData.mPointData[iZone].resize(mZonePointCounts[iZone]*mNPointComponents);
		for(std::size_t iElem=0; iElem<zone->GetnElem(); iElem++)
		{
			// VTK points are XYZ tuples, including for planar quadrilaterals.
			const auto& coordinates = zone->GetElementGeometry(iElem)->GetCoordSolDOFs();
			for(std::size_t iNode=0; iNode<nNode; iNode++)
			{
				const std::size_t i = 3*(iPoint+iNode);
				realData.mPoints[i  ] = static_cast<T>(coordinates(0,iNode));
				realData.mPoints[i+1] = static_cast<T>(coordinates(1,iNode));
				realData.mPoints[i+2] = static_cast<T>(0);
			}

			// Subdivide this element into linear quadrilaterals, as in the legacy writer.
			for(std::size_t j=0; j<nPoly; j++)
			{
				for(std::size_t i=0; i<nPoly; i++)
				{
					const std::size_t n0 = iPoint + j*nNode1D + i;
					const std::size_t k  = nCorner*iCell;
					mConnectivity[k  ] = static_cast<std::int32_t>(n0);
					mConnectivity[k+1] = static_cast<std::int32_t>(n0+1);
					mConnectivity[k+2] = static_cast<std::int32_t>(n0+1+nNode1D);
					mConnectivity[k+3] = static_cast<std::int32_t>(n0+nNode1D);
					mOffsets[iCell++]  = static_cast<std::int32_t>(k+nCorner);
				}
			}
			iPoint += nNode;
		}
	}
}

//-----------------------------------------------------------------------------------

template<typename T>
void CXMLBinaryVTK::InitializeFileHeader
(
 const CRealData<T> &realData
)
 /*
	* Function that declares the arrays and their byte offsets relative to the appended-data marker.
	*/
{
	static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big);

	constexpr const char *vtkType = std::is_same_v<T,float> ? "Float32" : "Float64";
	const char* byteOrder = std::endian::native == std::endian::little ? "LittleEndian" : "BigEndian";
	const std::size_t nPoint = realData.mPoints.size()/3;

	std::ostringstream header;
	header.imbue(std::locale::classic());
	header << "<?xml version=\"1.0\"?>\n"
	       << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"" << byteOrder
	       << "\" header_type=\"UInt64\">\n"
	       << "  <UnstructuredGrid>\n"
	       << "    <FieldData>\n"
	       << "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\""
	       << " format=\"appended\" offset=\"0\"/>\n"
	       << "    </FieldData>\n"
	       << "    <Piece NumberOfPoints=\"" << nPoint << "\" NumberOfCells=\"" << mCellTypes.size() << "\">\n";

	// Every array starts with an unsigned 64-bit payload length followed by the raw bytes.
	// TimeValue is the first appended array, so the point fields start after it.
	std::uint64_t offset = sizeof(std::uint64_t) + sizeof(double);
	auto lDeclareArray = [&](const char    *type,
	                         const char    *name,
	                         std::size_t    nComponent,
	                         std::uint64_t  nByte)
	{
		header << "        <DataArray type=\"" << type << "\" Name=\"" << name
		       << "\" NumberOfComponents=\"" << nComponent
		       << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";

		constexpr std::uint64_t maxOffset = std::numeric_limits<std::uint64_t>::max();
		if( offset > maxOffset-sizeof(std::uint64_t) ||
		    nByte > maxOffset-offset-sizeof(std::uint64_t) )
			ERROR("VTK XML appended data exceeds the 64-bit offset range.");
		offset += sizeof(std::uint64_t) + nByte;
	};

	header << "      <PointData>\n";
	for(const auto& field : mPointFields)
		lDeclareArray(vtkType, field.mName, field.mNComponent, CountBytes(nPoint, field.mNComponent*sizeof(T)));
	header << "      </PointData>\n"
	       << "      <Points>\n";
	lDeclareArray(vtkType, "Points", 3, CountBytes(realData.mPoints.size(), sizeof(T)));
	header << "      </Points>\n"
	       << "      <Cells>\n";
	lDeclareArray("Int32", "connectivity", 1, CountBytes(mConnectivity.size(), sizeof(std::int32_t)));
	lDeclareArray("Int32", "offsets",      1, CountBytes(mOffsets.size(), sizeof(std::int32_t)));
	lDeclareArray("UInt8", "types",         1, CountBytes(mCellTypes.size(), sizeof(std::uint8_t)));
	header << "      </Cells>\n"
	       << "    </Piece>\n"
	       << "  </UnstructuredGrid>\n"
	       << "  <AppendedData encoding=\"raw\">\n"
	       << "    _";

	mFileHeader = header.str();
}

//-----------------------------------------------------------------------------------

void CXMLBinaryVTK::WriteFileVTK
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const COpenMP            *openmp_container [[maybe_unused]],
 const CMultizoneSolver   *multizone_solver_container,
 as3double                 time
)
 /*
	* Function that dispatches XML binary output to the configured precision.
	*/
{
	std::visit([&](auto &realData)
	{
		WriteSnapshot(config_container, multizone_geometry_container, multizone_solver_container, time, realData);
	}, mRealData);
}

//-----------------------------------------------------------------------------------

template<typename T>
void CXMLBinaryVTK::WriteSnapshot
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container,
 as3double                 time,
 CRealData<T>             &realData
)
 /*
	* Function that writes a VTK XML snapshot in native byte order.
	*/
{
	std::cout << "--------------------------------------------------------------------------------------------\n"
	          << "Writing solution in VTK XML binary format... ";

	// Geometry and the array layout are fixed for the lifetime of the writer.
	if( GetFileNumber() == 0 )
	{
		InitializeMesh(multizone_geometry_container, realData);
		InitializeFileHeader(realData);
	}
	DetermineVisualizationData(multizone_geometry_container, multizone_solver_container, realData.mPointData);

	const std::string filename = std::string(config_container->GetOutputVisFilename()) + "_"
	                           + std::to_string(GetFileNumber()) + ".vtu";
	std::FILE *file = std::fopen(filename.c_str(), "wb");
	if( !file ) ERROR("Could not open VTK XML visualization file: " + filename);

	WriteValues(file, mFileHeader.data(), mFileHeader.size());

	// TimeValue is a dataset field and always uses double precision.
	const double timeValue = static_cast<double>(time);
	const std::uint64_t timeBytes = sizeof(timeValue);
	WriteValues(file, &timeBytes, 1);
	WriteValues(file, &timeValue, 1);

	// Each field is one appended array. Write its zone segments directly from the reusable buffers.
	const std::size_t nPoint = realData.mPoints.size()/3;
	for(const auto& field : mPointFields)
	{
		const std::uint64_t nByte = CountBytes(nPoint, field.mNComponent*sizeof(T));
		WriteValues(file, &nByte, 1);
		for(std::size_t iZone=0; iZone<mZonePointCounts.size(); iZone++)
		{
			const std::size_t nZonePoint = mZonePointCounts[iZone];
			if( nZonePoint == 0 ) continue;
			const T *data = realData.mPointData[iZone].data() + field.mComponentOffset*nZonePoint;
			WriteValues(file, data, field.mNComponent*nZonePoint);
		}
	}

	WriteArray(file, realData.mPoints);
	WriteArray(file, mConnectivity);
	WriteArray(file, mOffsets);
	WriteArray(file, mCellTypes);

	constexpr std::string_view footer = "\n  </AppendedData>\n</VTKFile>\n";
	WriteValues(file, footer.data(), footer.size());
	if( std::fclose(file) != 0 ) ERROR("Could not close VTK XML visualization file: " + filename);

	WriteSeriesIndex(config_container, filename, time);
	std::cout << "Done.\n"
	          << "--------------------------------------------------------------------------------------------" << std::endl;
}
