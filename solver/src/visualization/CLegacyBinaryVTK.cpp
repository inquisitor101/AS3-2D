#include "visualization/CLegacyBinaryVTK.hpp"
#include "data/CMatrixAS3.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "CConfig.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "geometry/CElementGeometry.hpp"
#include "io/NInputUtility.hpp"


//-----------------------------------------------------------------------------------
// CLegacyBinaryVTK member functions.
//-----------------------------------------------------------------------------------


CLegacyBinaryVTK::CLegacyBinaryVTK
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
	:
		IFileVTK(config_container, multizone_geometry_container)
 /*
	* Constructor that initializes legacy binary VTK output.
	*/
{
	// Check for big endian. If not, the bytes must be swapped when
  // writing a binary VTK file.
  union {std::uint32_t i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) mBigEndian = false;
  else                  mBigEndian = true;

}

//-----------------------------------------------------------------------------------

CLegacyBinaryVTK::~CLegacyBinaryVTK
(
 void
)
 /*
	* Destructor for CLegacyBinaryVTK.
	*/
{

}

//-----------------------------------------------------------------------------------

void CLegacyBinaryVTK::WriteFileVTK
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const COpenMP            *openmp_container [[maybe_unused]],
 const CMultizoneSolver   *multizone_solver_container,
 as3double                 time
)
 /*
	* Function that dispatches legacy binary output to the configured precision.
	*/
{
	if( config_container->GetOutputVisPrecision() == EFloatingPointPrecision::DOUBLE )
		WriteSnapshot<double>(config_container, multizone_geometry_container, multizone_solver_container, time);
	else
		WriteSnapshot<float>(config_container, multizone_geometry_container, multizone_solver_container, time);
}

//-----------------------------------------------------------------------------------

template<typename T>
void CLegacyBinaryVTK::WriteSnapshot
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container,
 as3double                 time
)
 /*
	* Function that writes a legacy binary VTK file for the current multizone solution.
	*/
{
 	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in legacy binary VTK format... ";

	// Determine number of zones.
	const std::size_t nZone = config_container->GetnZone();

	// Create VTK filename.
	const std::string filename = std::string(config_container->GetOutputVisFilename()) + "_"
	                           + std::to_string(GetFileNumber()) + ".vtk";


	// Define the maximum string length for the writing.
  constexpr std::size_t MAX_STRING_LENGTH = 255;
	constexpr const char *vtkType = std::is_same_v<T,float> ? "float" : "double";

	// Initialize the buffer data.
	as3vector2d<T>         vars_buf(nZone);
	as3vector2d<T>         coor_buf(nZone);
	as3vector2d<std::int32_t>  conn_buf(nZone);
	as3vector2d<std::int32_t>  type_buf(nZone);

	// AS3 uses quadrilateral elements only. Based on VTK syntax, these have the below properties.
	const std::int32_t quad_elem = 9;
	const std::size_t quad_npts = 4;

	// Accumulate the total number of DOFs written per zone.
	std::size_t nDOFsWritten    = 0;
	// Accumulate the total number of P1 sub-elements per zone.
	std::size_t nSubElemWritten = 0;


	// Assemble the data to write.
	for(std::size_t iZone=0; iZone<nZone; iZone++)
	{
		// Extract the relevant zone.
		auto* zone = multizone_geometry_container->GetSinglezoneGeometry(iZone);

		// Extract the current zone information.	
		const std::size_t nElem = zone->GetnElem();
		const std::size_t nNode = zone->GetnNodeGrid2D();
		const std::size_t nPoly = zone->GetnPolyGrid();
		
    // Validate products before forming them; legacy VTK uses 32-bit connectivity.
    const std::size_t vtkMax = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    if((nElem != 0 && nNode > vtkMax/nElem) ||
       (nPoly != 0 && nPoly > vtkMax/nPoly))
      ERROR("Mesh exceeds the legacy VTK 32-bit connectivity range.");

		// Total DOFs and P1 sub-elements in this zone.
		const std::size_t nDOFsTot = nNode*nElem;
		const std::size_t nSubElem = nPoly*nPoly;
    if(nDOFsWritten > vtkMax || nDOFsTot > vtkMax-nDOFsWritten ||
       nSubElemWritten > vtkMax/(quad_npts+1) ||
       (nSubElem != 0 && nElem > (vtkMax/(quad_npts+1)-nSubElemWritten)/nSubElem))
      ERROR("Mesh exceeds the legacy VTK 32-bit connectivity range.");

		// Extract number of elements in x and y in this zone.
		const std::size_t niElem = zone->GetniElem();
		const std::size_t njElem = zone->GetnjElem();

		// Deduce number of nodes for the DOFs in 1D.
		const std::size_t nDOFsSol1D = nPoly + 1;

		// Ensure the number of elements is correct.
		if( niElem*njElem != nElem ) ERROR("Inconsistency in the number of elements being written.");

		// Initialize and specify the coordinates. Note, ParaView expects them in 3D.
		coor_buf[iZone].resize( 3*nDOFsTot );
		// Initialize the global connectivity matrix.
		// Note, the addition of one is for the number of points per element.
		conn_buf[iZone].resize( nElem*nSubElem*(quad_npts+1) );
		// Initialize the element types. For this code, just specify a quadrilateral.
		type_buf[iZone].resize( nElem*nSubElem, quad_elem );

		// Initialize the buffer data.
  	for(std::size_t ijElem=0; ijElem<nElem; ++ijElem)
  	{
			// Deduce the local element indices in x and y.
			//const int jElem = ijElem/niElem;
			//const int iElem = ijElem - jElem*niElem;

			// Deduce the current element coordinate array location.
			T            *coorelem = coor_buf[iZone].data() + ijElem*3*nNode;  
			// Deduce the current element connectivity array location.
			std::int32_t *connelem = conn_buf[iZone].data() + ijElem*(quad_npts+1)*nSubElem;

			// Extract the current element coordinates.
			auto& xyelem = zone->GetElementGeometry(ijElem)->GetCoordSolDOFs();

			// Loop over every DOF on this element and store the coordinate.
			std::size_t idx = 0;
			for(std::size_t l=0; l<nNode; l++)
			{
				coorelem[idx++] = static_cast<T>( xyelem(0,l) );
				coorelem[idx++] = static_cast<T>( xyelem(1,l) );
				coorelem[idx++] = static_cast<T>(0);
			}

			// Deduce the start of the current element's connectivity array.
			const std::size_t offset = ijElem*nNode + nDOFsWritten;
			
			// Connectivity array.
			idx = 0;
			for(std::size_t j=0; j<nPoly; j++)
			{
				const std::size_t joff = offset + j*nDOFsSol1D;
				for(std::size_t i=0; i<nPoly; i++)
				{
					// Starting index of each nPoly=1 sub-element.
					const std::size_t n0 = joff + i; 
					// Number of points in this element.
					connelem[idx++] = static_cast<std::int32_t>(quad_npts);
					// Actual connectivity array.
					connelem[idx++] = static_cast<std::int32_t>(n0);
					connelem[idx++] = static_cast<std::int32_t>(n0 + 1);
					connelem[idx++] = static_cast<std::int32_t>(n0 + 1 + nDOFsSol1D);
					connelem[idx++] = static_cast<std::int32_t>(n0 + nDOFsSol1D);
				}
			}
		}
	
		// Update the total number of DOFs written, thus far.
		nDOFsWritten    += nDOFsTot;
		// Update the total number of P1 sub-elements written, thus far.
		nSubElemWritten += nElem*nSubElem; 

		// Allocate the size of the written data variables.
		vars_buf[iZone].resize(nDOFsTot*mNPointComponents, static_cast<T>(0));
	}

	// Compute and store the required data for visualization. 
	// Note, this is done over all elements in all zones simulataneoously.
	DetermineVisualizationData(multizone_geometry_container, multizone_solver_container, vars_buf); 
	

	// Check if there need be any swapping, since ParaView expects data in big endian format.
	if( !mBigEndian )
	{
		for(std::size_t iZone=0; iZone<nZone; iZone++)
		{
			NInputUtility::SwapBytes( coor_buf[iZone].data(), sizeof(T),        coor_buf[iZone].size() );
			NInputUtility::SwapBytes( vars_buf[iZone].data(), sizeof(T),        vars_buf[iZone].size() );
			NInputUtility::SwapBytes( conn_buf[iZone].data(), sizeof(std::int32_t), conn_buf[iZone].size() );
			NInputUtility::SwapBytes( type_buf[iZone].data(), sizeof(std::int32_t), type_buf[iZone].size() );
		}
	}


	// Save number of total DOFs per zone, as it is needed for writing scalar/vector data.
	std::vector<std::size_t> nDOFsZone(nZone);
	for(std::size_t iZone=0; iZone<nZone; iZone++)
	{
		const std::size_t nElem = multizone_geometry_container->GetSinglezoneGeometry(iZone)->GetnElem();
		const std::size_t nNode = multizone_geometry_container->GetSinglezoneGeometry(iZone)->GetnNodeGrid2D();
		nDOFsZone[iZone] = nElem*nNode; 
	} 


	/*
	 * Write the data to file.
	 */


	// Open the visualization file for binary writing.
	std::FILE *fh = std::fopen(filename.c_str(), "wb");
	if( !fh ) ERROR("Visualization file in legacy + binary format could not be opened.");


	// Write the header of the visualization file.
  char str_buf[MAX_STRING_LENGTH];
  std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
	std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);

  std::strcpy(str_buf, "vtk output\n");
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);

  std::strcpy(str_buf, "BINARY\n");
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);

  std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);


	// Write the snapshot time as one double in dataset-level field data.
	double timeValue = static_cast<double>(time);
	if( !mBigEndian ) NInputUtility::SwapBytes(&timeValue, sizeof(timeValue), 1);
	if( std::fputs("FIELD FieldData 1\nTimeValue 1 1 double\n", fh) == EOF ||
	    std::fwrite(&timeValue, sizeof(timeValue), 1, fh) != 1 ||
	    std::fputc('\n', fh) == EOF )
		ERROR("Could not write the legacy VTK snapshot time.");

  // Write the coordinates.
  std::snprintf(str_buf, sizeof(str_buf), "POINTS %zu %s\n", nDOFsWritten, vtkType);
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
	for( auto& coor: coor_buf ) std::fwrite(coor.data(), sizeof(T), coor.size(), fh);

  // Write the connectivity data.
  std::snprintf(str_buf, sizeof(str_buf), "\nCELLS %zu %zu\n", nSubElemWritten, (quad_npts+1)*nSubElemWritten);
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
  for( auto& conn: conn_buf ) std::fwrite(conn.data(), sizeof(std::int32_t), conn.size(), fh);

  // Write the element type data.
  std::snprintf(str_buf, sizeof(str_buf), "\nCELL_TYPES %zu\n", nSubElemWritten);
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
  for( auto& type: type_buf ) std::fwrite(type.data(), sizeof(std::int32_t), type.size(), fh);

  // Write the ASCII line for the point data.
  std::snprintf(str_buf, sizeof(str_buf), "\nPOINT_DATA %zu\n", nDOFsWritten);
  std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);


	// Write each scalar or interleaved vector field across all zones.
	for(const auto& field : mPointFields)
	{
		if( field.mNComponent == 3 )
		{
			std::snprintf(str_buf, sizeof(str_buf), "\nVECTORS %s %s\n", field.mName, vtkType);
			std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
		}
		else
		{
			std::snprintf(str_buf, sizeof(str_buf), "\nSCALARS %s %s 1\n", field.mName, vtkType);
			std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
			std::strcpy(str_buf, "LOOKUP_TABLE default\n");
			std::fwrite(str_buf, sizeof(char), std::strlen(str_buf), fh);
		}

		for(std::size_t iZone=0; iZone<nZone; iZone++)
		{
			const std::size_t nPoint = nDOFsZone[iZone];
			const T *data = vars_buf[iZone].data() + field.mComponentOffset*nPoint;
			std::fwrite(data, sizeof(T), field.mNComponent*nPoint, fh);
		}
	}


	// Close the snapshot and add it to the time-series index.
	const bool writeFailed = std::ferror(fh) != 0;
	if( std::fclose(fh) != 0 || writeFailed ) ERROR("Could not write legacy VTK snapshot: " + filename);
	WriteSeriesIndex(config_container, filename, time);

	// Report: all is complete!
	std::cout << "Done." << std::endl;
	std::cout << "----------------------------------------------"
							 "----------------------------------------------" << std::endl;
}
