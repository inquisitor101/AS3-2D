#include "log_structure.hpp"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>



//-----------------------------------------------------------------------------------
// NLogger namespace functions.
//-----------------------------------------------------------------------------------


void NLogger::PrintInitSolver
(
 CConfig *config_container
)
 /*
	* Function that prints the equation type used in each zone.
	*/
{
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	
	// Report specific information about solver type and buffer layer, per zone.
	std::cout << "Initiating simulation... " << std::endl;

	// Extract total number of zones.
	auto nZone = config_container->GetnZone();
	
	// Define temporary lambda that prints the type of buffer layer in a zone.
	auto LPrintTypeBufferLayer = [&](ETypeBufferLayer buf) -> std::string
	{
		// String containing the relevant message.
		std::string out;
		
		// Check which type of message to issue, depending on the type of buffer layer.
		switch(buf)
		{
			case(ETypeBufferLayer::NONE): { out = " no buffer layer"; break; }
			default: ERROR("Unkown buffer layer.");
		}

		// return the name of the buffer layer, if any.
		return out;
	};


	// Loop over each zone and print the corresponding solver specs.
	for(std::size_t iZone=0; iZone<nZone; iZone++)
	{
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
	      std::cout << "  Beginning EE Solver in iZone(" << iZone << "): "
					        << LPrintTypeBufferLayer(config_container->GetTypeBufferLayer(iZone));
	    	break;
			}

			default: ERROR("Unknown solver type.");
		}

		std::cout << std::endl;
	}

	std::cout << "Done." << std::endl;
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayBoundaryConditions
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container
)
 /*
	* Function that prints the configured physical and periodic boundary conditions.
	*/
{
  struct CBoundarySummary
  {
    ETypeBC mType;
    std::size_t mFaces;
  };
  std::map<std::string, CBoundarySummary> physical;

  // Aggregate initialized boundary batches across zones and both directions.
  // This reports the conditions actually used by the solver, once per name.
  for(std::size_t direction=0; direction<2; ++direction)
  {
    const auto& batches = direction == 0
      ? multizone_geometry_container->GetBoundaryBatchesIFace()
      : multizone_geometry_container->GetBoundaryBatchesJFace();
    for(std::size_t batch=0; batch<batches.size(); ++batch)
    {
      const auto index = direction == 0
        ? multizone_solver_container->GetIndexBoundaryIFace(batch)
        : multizone_solver_container->GetIndexBoundaryJFace(batch);
      const auto type = multizone_solver_container->GetBoundary(index)->GetTypeBC();
      const auto& faces = batches[batch];
      auto [entry, inserted] = physical.emplace(faces.GetName(), CBoundarySummary{type, 0});
      if( !inserted && entry->second.mType != type )
        ERROR("Inconsistent boundary types for marker: " + faces.GetName());
      entry->second.mFaces += faces.GetnFace();
    }
  }

  struct CReportRow
  {
    std::string mName, mCondition;
    std::size_t mFaces;
    std::string mParameters;
  };
  as3vector1d<CReportRow> rows;
  for(const auto& [name, boundary] : physical)
  {
    const auto* parameters = config_container->GetBoundaryParamMarker(name);
    if( parameters->GetTypeBC() != boundary.mType )
      ERROR("Boundary parameters do not match solver for marker: " + name);
    std::string condition;
    std::ostringstream brief;
    brief << std::setprecision(8);
    switch( boundary.mType )
    {
      case ETypeBC::SYMMETRY:
        condition = "Symmetry";
        brief << '-';
        break;
      case ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN:
      {
        condition = "Static mass-flow inlet (Riemann)";
        const auto* inlet = static_cast<const CStaticMassFlowInletRiemannParamMarker*>(parameters);
        brief << "density = "   << inlet->mDensity     << ", " 
              << "velocity = (" << inlet->mVelocity[0] << ", " << inlet->mVelocity[1] << ")";
        break;
      }
      case ETypeBC::PRESSURE_OUTLET_RIEMANN:
      {
        condition = "Pressure outlet (Riemann)";
        const auto* outlet = static_cast<const CPressureOutletRiemannParamMarker*>(parameters);
        brief << "pressure = " << outlet->mPressure;
        break;
      }
      case ETypeBC::FARFIELD:
      {
        condition = "Farfield";
        const auto* farfield = static_cast<const CFarfieldParamMarker*>(parameters);
        brief << "density = "   << farfield->mDensity     << ", "
              << "velocity = (" << farfield->mVelocity[0] << ", " << farfield->mVelocity[1] << "), "
              << "pressure = "  << farfield->mPressure;
        break;
      }
      default: ERROR("Unsupported physical boundary type for marker: " + name);
    }
    rows.push_back({name, condition, boundary.mFaces, brief.str()});
  }

  // Periodic markers are interfaces rather than physical boundary evaluators.
  // Count each paired face once and retain the configured translation direction.
  const auto& periodic = config_container->GetPeriodicParamMarker();
  for(const auto& parameters : periodic)
  {
    std::size_t pairs = 0;
    for(const auto& family : multizone_geometry_container->GetPeriodicFamilyMarkers())
      for(const auto& marker : family.GetMarkers())
        if( marker.mNameI == parameters->mName && marker.mNameJ == parameters->mNameMatching )
        {
          const auto* zone = multizone_geometry_container->GetSinglezoneGeometry(marker.mIndexZoneI);
          pairs += marker.GetiTypeFace() == ETypeFace::IFACE ? zone->GetnjElem() : zone->GetniElem();
        }
    std::ostringstream brief;
    brief << std::setprecision(8) << "translation=(" << parameters->mVectorTrans[0]
          << ", " << parameters->mVectorTrans[1] << ')';
    rows.push_back({parameters->mName + " -> " + parameters->mNameMatching,
                    "Periodic", pairs, brief.str()});
  }

  std::size_t nameWidth = std::string("Marker").size();
  std::size_t conditionWidth = std::string("Condition").size();
  std::size_t faceWidth = std::string("Faces").size();
  std::size_t parameterWidth = std::string("Parameters").size();
  for(const auto& row : rows)
  {
    nameWidth = std::max(nameWidth, row.mName.size());
    conditionWidth = std::max(conditionWidth, row.mCondition.size());
    faceWidth = std::max(faceWidth, std::to_string(row.mFaces).size());
    parameterWidth = std::max(parameterWidth, row.mParameters.size());
  }

  // A local stream keeps this report independent of earlier scientific output
  // and preserves the caller's precision, fill character and alignment.
  std::ostringstream report;
  report << std::string(92, '-') << "\nBoundary conditions:\n";
  if( rows.empty() ) report << "  None.\n";
  else
  {
    const auto row = [&](const auto& name, const auto& condition,
                         const auto& faces, const auto& parameters)
    {
      report << "  " << std::left << std::setw(static_cast<int>(nameWidth + 2)) << name
             << std::setw(static_cast<int>(conditionWidth + 2)) << condition
             << std::right << std::setw(static_cast<int>(faceWidth)) << faces
             << "  " << parameters << '\n';
    };
    row("Marker", "Condition", "Faces", "Parameters");
    report << "  " << std::string(nameWidth + conditionWidth + faceWidth + parameterWidth + 6, '-') << '\n';
    for(const auto& entry : rows)
      row(entry.mName, entry.mCondition, entry.mFaces, entry.mParameters);
    if( !periodic.empty() )
      report << "  Periodic rows count face pairs; translation maps the first marker to the second.\n";
  }
  std::cout << report.str();
}

//-----------------------------------------------------------------------------------

void NLogger::MonitorOutput
(
 CConfig      *config_container,
 CMonitorData *monitor_container,
 std::size_t   iter,
 as3double     time,
 as3double     step
)
 /*
	* Function that prints the monitoring header and current simulation values.
	*/
{
	// Extract the max number of iterations.
	std::size_t nMaxIter = config_container->GetMaxIterTime();
	// Number of output reports for monitoring progress.
	std::size_t nOutput  = std::max(std::size_t{1}, nMaxIter/100);
	// Compute number of max digits needed for output.
	std::size_t nDigits  = std::to_string(nMaxIter).size();

	// Display header.
	if( iter%(50*nOutput) == 0 )
	{
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
		std::cout << " Iteration "     << "\t"
			        << " Time "          << "\t" 
							<< " Sync time "     << "\t"
							<< " Steps "         << "\t"
							<< " min(dt) "       << "\t"
							<< " max(dt) "       << "\t"
							<< " Max(Mach) "     << "\n";
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
	}

  // Extract the maximum Mach number.
  const as3double   Mmax = monitor_container->mMachMax;
	// Extract the number of substeps per sync.
	const std::size_t nsub = monitor_container->mNSyncSubStep;
	// Extract the min and max time steps, per sync step.
	const as3double dtmin = monitor_container->mMinTimeStep;
	const as3double dtmax = monitor_container->mMaxTimeStep;

	// Display progress.
	std::cout << std::scientific 
						<< std::setprecision(6)
		        << " " 
						<< std::setw(static_cast<int>(nDigits)) << iter
						<< "\t"  << time
						<< "\t"  << step
						<< "\t"  << nsub
						<< "\t"  << dtmin
						<< "\t"  << dtmax
            << "\t"  << Mmax
						<< std::endl;
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayOpenMPInfo
(
 COpenMP                  *openmp_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container
)
 /*
  * Function that prints thread counts and estimated static work distribution.
  */
{
  (void)openmp_container;
  // Report output.
  std::cout << "----------------------------------------------"
               "----------------------------------------------\n";
#ifdef HAVE_OPENMP
  // Allocate for the maximum team; report the actual team observed below.
  const std::size_t nThreads = static_cast<std::size_t>(omp_get_max_threads());
  std::size_t actualThreads = 1;

  // Get the total number of elements in all zones.
  const std::size_t nElemTotal = multizone_geometry_container->GetnElemTotal();
  // Get the total number of i-faces in all zones.
  const std::size_t nIFace     = multizone_geometry_container->GetnIFace();
  // Get the total number of j-faces in all zones.
  const std::size_t nJFace     = multizone_geometry_container->GetnJFace();

  // Estimate computational work load of each thread.
  as3vector1d<std::size_t> workloadDOFs(nThreads, 0);
  as3vector1d<std::size_t> workloadElem(nThreads, 0);
  as3vector1d<std::size_t> workloadIFace(nThreads, 0);
  as3vector1d<std::size_t> workloadJFace(nThreads, 0);

  as3vector1d<std::size_t> workloadIFaceInternal_standard(nThreads,  0);
  as3vector1d<std::size_t> workloadIFaceBoundary_standard(nThreads,  0);
  as3vector1d<std::size_t> workloadIFaceInterface_standard(nThreads, 0);

  as3vector1d<std::size_t> workloadIFaceInternal_balanced(nThreads,  0);
  as3vector1d<std::size_t> workloadIFaceBoundary_balanced(nThreads,  0);
  as3vector1d<std::size_t> workloadIFaceInterface_balanced(nThreads, 0);

  as3vector1d<std::size_t> workloadJFaceInternal_standard(nThreads,  0);
  as3vector1d<std::size_t> workloadJFaceBoundary_standard(nThreads,  0);
  as3vector1d<std::size_t> workloadJFaceInterface_standard(nThreads, 0);

  as3vector1d<std::size_t> workloadJFaceInternal_balanced(nThreads,  0);
  as3vector1d<std::size_t> workloadJFaceBoundary_balanced(nThreads,  0);
  as3vector1d<std::size_t> workloadJFaceInterface_balanced(nThreads, 0);

  // Use one actual team for all three static loops and the before/after counts.
#pragma omp parallel
  {
#pragma omp single
    actualThreads = static_cast<std::size_t>(omp_get_num_threads());

    // Estimate the i-surface workload.
#pragma omp for schedule(static)
    for(std::size_t i=0; i<nIFace; i++)
    {
      // Thread index.
      const std::size_t iThread = static_cast<std::size_t>(omp_get_thread_num());

      // Accumulate the number of elements per thread.
      workloadIFace[iThread]++;

      // Get the relevant face information for the standard and load-balanced partitions.
      const auto face_info_standard = multizone_geometry_container->GetFlattenedIndexIFace(i);
      const auto face_info_balanced = multizone_geometry_container->GetFlattenedIndexIFaceLoadBalanced(i);

      // Estimate the work load for the standard approach.
      switch( face_info_standard.mFaceType )
      {
        case( ETypeFaceGeometry::INTERNAL ):  { workloadIFaceInternal_standard[iThread]++;  break; }
        case( ETypeFaceGeometry::BOUNDARY ):  { workloadIFaceBoundary_standard[iThread]++;  break; }
        case( ETypeFaceGeometry::INTERFACE ): { workloadIFaceInterface_standard[iThread]++; break; }
        default: ERROR("Unknown face type.");
      }

      // Then, estimate the workload for the load-balanced approach.
      switch( face_info_balanced.mFaceType )
      {
        case( ETypeFaceGeometry::INTERNAL ):  { workloadIFaceInternal_balanced[iThread]++;  break; }
        case( ETypeFaceGeometry::BOUNDARY ):  { workloadIFaceBoundary_balanced[iThread]++;  break; }
        case( ETypeFaceGeometry::INTERFACE ): { workloadIFaceInterface_balanced[iThread]++; break; }
        default: ERROR("Unknown face type.");
      }
    }

    // Estimate the j-surface workload.
#pragma omp for schedule(static)
    for(std::size_t i=0; i<nJFace; i++)
    {
      // Thread index.
      const std::size_t iThread = static_cast<std::size_t>(omp_get_thread_num());

      // Accumulate the number of elements per thread.
      workloadJFace[iThread]++;
  
      // Get the relevant face information for the standard and load-balanced partitions.
      const auto face_info_standard = multizone_geometry_container->GetFlattenedIndexJFace(i);
      const auto face_info_balanced = multizone_geometry_container->GetFlattenedIndexJFaceLoadBalanced(i);

      // Estimate the work load for the standard approach.
      switch( face_info_standard.mFaceType )
      {
        case( ETypeFaceGeometry::INTERNAL ):  { workloadJFaceInternal_standard[iThread]++;  break; }
        case( ETypeFaceGeometry::BOUNDARY ):  { workloadJFaceBoundary_standard[iThread]++;  break; }
        case( ETypeFaceGeometry::INTERFACE ): { workloadJFaceInterface_standard[iThread]++; break; }
        default: ERROR("Unknown face type.");
      }

      // Then, estimate the workload for the load-balanced approach.
      switch( face_info_balanced.mFaceType )
      {
        case( ETypeFaceGeometry::INTERNAL ):  { workloadJFaceInternal_balanced[iThread]++;  break; }
        case( ETypeFaceGeometry::BOUNDARY ):  { workloadJFaceBoundary_balanced[iThread]++;  break; }
        case( ETypeFaceGeometry::INTERFACE ): { workloadJFaceInterface_balanced[iThread]++; break; }
        default: ERROR("Unknown face type.");
      }
    }

    // Estimate the elements' workload.
#pragma omp for schedule(static)
    for(std::size_t i=0; i<nElemTotal; i++)
    {
      // Extract the element indices.
      const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

      // Deduce the current element's zone and index.
      const auto iZone = elem_info.mIndexZone; 

      // Thread index.
      const std::size_t iThread = static_cast<std::size_t>(omp_get_thread_num());

      // Accumulate the number of elements per thread.
      workloadElem[iThread]++;

      // Accumulate the number of solution DOFs per thread.
      workloadDOFs[iThread] += multizone_solver_container->GetSinglezoneSolver(iZone)->GetStandardElement()->GetnSol2D(); 
    }


  } // End of the parallel region.

  const auto& iPlan = multizone_geometry_container->GetIFaceLoadBalancedPermutation();
  const auto& jPlan = multizone_geometry_container->GetJFaceLoadBalancedPermutation();
  const auto load = [](const CLoadBalancedFacePermutation& plan,
                       std::size_t internal, std::size_t boundary, std::size_t interface)
  {
    return static_cast<double>(internal)*plan.GetFaceWeight(ETypeFaceGeometry::INTERNAL)
         + static_cast<double>(boundary)*plan.GetFaceWeight(ETypeFaceGeometry::BOUNDARY)
         + static_cast<double>(interface)*plan.GetFaceWeight(ETypeFaceGeometry::INTERFACE);
  };
  const auto number = [](double value)
  {
    std::ostringstream text;
    text << value;
    return text.str();
  };

  // Build the table independently of std::cout's current precision, alignment,
  // fill character, or scientific formatting. No tabs: all rows share widths.
  std::ostringstream report;
  report << "OpenMP threads: " << actualThreads;
  if( actualThreads != nThreads ) report << " (configured maximum: " << nThreads << ')';
  report << "\nFace load balancing: ";
  if( iPlan.GetStrategy() == jPlan.GetStrategy() ) report << iPlan.GetStrategyName();
  else report << "I-faces: " << iPlan.GetStrategyName() << "; J-faces: " << jPlan.GetStrategyName();
  report << "\nEstimated work distribution:\n";

  constexpr int labelWidth = 22;
  for(std::size_t i=0; i<actualThreads; ++i)
  {
    const auto iStandard = number(load(iPlan, workloadIFaceInternal_standard[i],
                         workloadIFaceBoundary_standard[i], workloadIFaceInterface_standard[i]));
    const auto iBalanced = number(load(iPlan, workloadIFaceInternal_balanced[i],
                         workloadIFaceBoundary_balanced[i], workloadIFaceInterface_balanced[i]));
    const auto jStandard = number(load(jPlan, workloadJFaceInternal_standard[i],
                         workloadJFaceBoundary_standard[i], workloadJFaceInterface_standard[i]));
    const auto jBalanced = number(load(jPlan, workloadJFaceInternal_balanced[i],
                         workloadJFaceBoundary_balanced[i], workloadJFaceInterface_balanced[i]));
    const auto width = std::max({std::size_t{8}, std::to_string(workloadIFace[i]).size(),
                 std::to_string(workloadJFace[i]).size(), iStandard.size(), iBalanced.size(),
                 jStandard.size(), jBalanced.size()});
    const int columnWidth = static_cast<int>(width) + 2;

    const auto row = [&](const char* label, const auto& standard, const auto& balanced)
    {
      report << "  " << std::left << std::setw(labelWidth) << label
             << std::right << std::setw(columnWidth) << standard
             << std::setw(columnWidth) << balanced << '\n';
    };
    
    report << std::string(2 + labelWidth + 2 * columnWidth, '-') << '\n';
    report << "(*) iThread: " << i << '\n';

    row("", "Standard", "Balanced");
    row("", "--------", "--------"); 
    
    row("(->) IFaces", workloadIFace[i], workloadIFace[i]);
    row("       Internal", workloadIFaceInternal_standard[i], workloadIFaceInternal_balanced[i]);
    row("       Boundary", workloadIFaceBoundary_standard[i], workloadIFaceBoundary_balanced[i]);
    row("       Interface", workloadIFaceInterface_standard[i], workloadIFaceInterface_balanced[i]);
    row("       Weighted load", iStandard, iBalanced);
    
    row("(->) JFaces", workloadJFace[i], workloadJFace[i]);
    row("       Internal", workloadJFaceInternal_standard[i], workloadJFaceInternal_balanced[i]);
    row("       Boundary", workloadJFaceBoundary_standard[i], workloadJFaceBoundary_balanced[i]);
    row("       Interface", workloadJFaceInterface_standard[i], workloadJFaceInterface_balanced[i]);
    row("       Weighted load", jStandard, jBalanced);
    
    report << "  (->) nElem: " << workloadElem[i] << '\n'
           << "  (->) nDOFs: " << workloadDOFs[i] << '\n';

  }
  std::cout << report.str();
#else
  std::cout << "This is a serial implementation." << std::endl;
#endif
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayAS3BinaryMetadata
(
 const NImportFile::NAS3BinaryFile::CAS3BinaryMetadata &metadata,
 std::ostream                                          &out
)
 /*
	* Function that prints AS3 grid metadata and file identifier mappings.
	*/
{
  const auto savedFlags = out.flags();
  const auto savedFill  = out.fill();

  out << std::dec << std::setfill(' ')
      << "\nAS3 mesh metadata\n"
      << "File          : " << metadata.mFilename << '\n'
      << "Dimensions    : " << metadata.mNDim << '\n'
      << "Zones         : " << metadata.mZones.size() << '\n'
      << "Byte swapping : "
      << (metadata.mByteSwap ? "yes" : "no") << "\n\n";

  out << std::right
      << std::setw(12) << "Zone"
      << std::setw(12) << "Order"
      << std::setw(12) << "i elements"
      << std::setw(12) << "j elements" << '\n'
      << std::string(48, '-') << '\n';

  for (const auto& zone : metadata.mZones)
  {
    out << std::right
        << std::setw(12) << zone.mZoneID
        << std::setw(12) << zone.mNPoly
        << std::setw(12) << zone.mNiElem
        << std::setw(12) << zone.mNjElem << '\n';
  }

  auto lPrintMapping = [&](const char* title, const auto& mapping, auto enumName)
  {
    out << '\n' << title << '\n'
        << std::right << std::setw(10) << "File ID"
        << "  " << std::left << std::setw(28) << "C++ value"
        << '\n' << std::string(40, '-') << '\n';

    for(const auto& entry : mapping)
    {
      out << std::right << std::setw(10) << entry.first
          << "  " << std::left << std::setw(28)
          << enumName(entry.second) << '\n';
    }
  };

  const auto& mappings = metadata.mMappings;

  lPrintMapping("Face mapping", mappings.mFaces, [](EFaceLocation value) -> const char*
  {
    switch (value)
    {
      case EFaceLocation::IMIN: return "EFaceLocation::IMIN";
      case EFaceLocation::IMAX: return "EFaceLocation::IMAX";
      case EFaceLocation::JMIN: return "EFaceLocation::JMIN";
      case EFaceLocation::JMAX: return "EFaceLocation::JMAX";
    }
    return "UNKNOWN";
  });

  lPrintMapping("Boundary mapping", mappings.mBoundaries, [](ETypeZoneMarker value) -> const char*
  {
    switch (value)
    {
      case ETypeZoneMarker::INTERNAL: return "ETypeZoneMarker::INTERNAL";
      case ETypeZoneMarker::EXTERNAL: return "ETypeZoneMarker::EXTERNAL";
    }
    return "UNKNOWN";
  });

  lPrintMapping("Nodal distribution mapping", mappings.mNodalDistributions, [](ETypeDOF value) -> const char*
  {
    switch (value)
    {
      case ETypeDOF::EQD: return "ETypeDOF::EQD";
      case ETypeDOF::LGL: return "ETypeDOF::LGL";
    }
    return "UNKNOWN";
  });

  out.flags(savedFlags);
  out.fill(savedFill);
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayInternalMarkers
(
 const as3vector1d<CInternalFamilyMarker> &families
)
 /*
  * Function that prints internal zone-edge pairings and traversal orientation.
  */
{
  // Return a readable face location.
  auto lFaceName = [](EFaceLocation face) -> const char*
  {
    switch( face )
    {
      case( EFaceLocation::IMIN ): return "imin";
      case( EFaceLocation::IMAX ): return "imax";
      case( EFaceLocation::JMIN ): return "jmin";
      case( EFaceLocation::JMAX ): return "jmax";
      default: return "unknown";
    }
  };

  // Print a separator matching the table width.
  auto lPrintSeparator = [](std::size_t width)
  {
    std::cout << std::string(width, '-') << '\n';
  };

  // Preserve the original stream formatting.
  const auto flags = std::cout.flags();
  const auto fill = std::cout.fill();

  std::cout << std::dec << std::right << std::setfill(' ');

  // Total width of the table columns.
  const std::size_t tableWidth = 62;

  std::cout << "\nInternal marker families: " << families.size() << '\n';

  for(std::size_t iFamily=0; iFamily<families.size(); iFamily++)
  {
    const auto& family = families[iFamily];

    std::cout << "\nFamily[" << iFamily << "]"
              << ", markers: " << family.GetnMarkers() << '\n';

    // Display the table headings.
    lPrintSeparator(tableWidth);

    std::cout << std::setw(10) << "Marker"
              << std::setw(10) << "Zone I"
              << std::setw(10) << "Face I"
              << std::setw(10) << "Zone J"
              << std::setw(10) << "Face J"
              << std::setw(12) << "Reversed" << '\n';

    lPrintSeparator(tableWidth);

    // Display each marker in this family.
    for(std::size_t iMarker=0; iMarker<family.GetnMarkers(); iMarker++)
    {
      const auto& marker = family.GetMarker(iMarker);

      std::cout << std::setw(10) << iMarker
                << std::setw(10) << marker.mIndexZoneI
                << std::setw(10) << lFaceName(marker.mFaceLocationI)
                << std::setw(10) << marker.mIndexZoneJ
                << std::setw(10) << lFaceName(marker.mFaceLocationJ)
                << std::setw(12) << (marker.mIsReversed ? "yes" : "no")
                << '\n';
    }

    lPrintSeparator(tableWidth);
  }

  std::cout << std::endl;

  // Restore the original stream formatting.
  std::cout.flags(flags);
  std::cout.fill(fill);
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayExternalMarkers
(
 const as3vector1d<CExternalFamilyMarker> &families
)
 /*
  * Function that prints named physical boundary markers and their zone edges.
  */
{
  // Return a readable face location.
  auto lFaceName = [](EFaceLocation face) -> const char*
  {
    switch( face )
    {
      case( EFaceLocation::IMIN ): return "imin";
      case( EFaceLocation::IMAX ): return "imax";
      case( EFaceLocation::JMIN ): return "jmin";
      case( EFaceLocation::JMAX ): return "jmax";
      default: return "unknown";
    }
  };

  // Print a separator matching the table width.
  auto lPrintSeparator = [](std::size_t width)
  {
    std::cout << std::string(width, '-') << '\n';
  };

  // Preserve the original stream formatting.
  const auto flags = std::cout.flags();
  const auto fill = std::cout.fill();

  std::cout << std::dec << std::right << std::setfill(' ');

  // Total width of the table columns.
  const std::size_t tableWidth = 30;

  std::cout << "\nExternal marker families: " << families.size() << '\n';

  for(std::size_t iFamily=0; iFamily<families.size(); iFamily++)
  {
    const auto& family = families[iFamily];

    std::cout << "\nFamily[" << iFamily << "]"
              << ", name: \"" << family.GetName() << "\""
              << ", markers: " << family.GetnMarkers() << '\n';

    // Display the table headings.
    lPrintSeparator(tableWidth);

    std::cout << std::setw(10) << "Marker"
              << std::setw(10) << "Zone"
              << std::setw(10) << "Face" << '\n';

    lPrintSeparator(tableWidth);

    // Display each marker in this family.
    for(std::size_t iMarker=0; iMarker<family.GetnMarkers(); iMarker++)
    {
      const auto& marker = family.GetMarker(iMarker);

      std::cout << std::setw(10) << iMarker
                << std::setw(10) << marker.mIndexZone
                << std::setw(10) << lFaceName(marker.mFaceLocation)
                << '\n';
    }

    lPrintSeparator(tableWidth);
  }

  std::cout << std::endl;

  // Restore the original stream formatting.
  std::cout.flags(flags);
  std::cout.fill(fill);
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayPeriodicMarkers
(
 const as3vector1d<CPeriodicFamilyMarker> &families
)
 /*
  * Function that prints periodic zone-edge pairings, orientation and translations.
  */
{
  // Return a readable face location.
  auto lFaceName = [](EFaceLocation face) -> const char*
  {
    switch( face )
    {
      case( EFaceLocation::IMIN ): return "imin";
      case( EFaceLocation::IMAX ): return "imax";
      case( EFaceLocation::JMIN ): return "jmin";
      case( EFaceLocation::JMAX ): return "jmax";
      default: return "unknown";
    }
  };

  // Print a separator matching the table width.
  auto lPrintSeparator = [](std::size_t width)
  {
    std::cout << std::string(width, '-') << '\n';
  };

  // Preserve the original stream formatting.
  const auto flags = std::cout.flags();
  const auto precision = std::cout.precision();
  const auto fill = std::cout.fill();

  std::cout << std::dec << std::scientific << std::setprecision(12)
            << std::right << std::setfill(' ');

  // Determine the column width needed for the boundary names.
  std::size_t nameWidth = 8;
  for( const auto& family : families )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      nameWidth = std::max(nameWidth, marker.mNameI.size());
      nameWidth = std::max(nameWidth, marker.mNameJ.size());
    }
  }

  const int nameColumnWidth = static_cast<int>(nameWidth + 2);

  // Total width of the fixed columns and the two boundary-name columns.
  const std::size_t tableWidth = 110 + 2 * (nameWidth + 2);

  std::cout << "\nPeriodic marker families: " << families.size() << '\n';

  for(std::size_t iFamily=0; iFamily<families.size(); iFamily++)
  {
    const auto& family = families[iFamily];

    std::cout << "\nFamily[" << iFamily << "]"
              << ", markers: " << family.GetnMarkers() << '\n';

    // Display the table headings.
    lPrintSeparator(tableWidth);

    std::cout << std::setw(10) << "Marker"
              << std::setw(nameColumnWidth) << "Name I"
              << std::setw(10) << "Zone I"
              << std::setw(10) << "Face I"
              << std::setw(nameColumnWidth) << "Name J"
              << std::setw(10) << "Zone J"
              << std::setw(10) << "Face J"
              << std::setw(12) << "Reversed"
              << std::setw(24) << "Translation X (I->J)"
              << std::setw(24) << "Translation Y (I->J)" << '\n';

    lPrintSeparator(tableWidth);

    // Display each marker in this family.
    for(std::size_t iMarker=0; iMarker<family.GetnMarkers(); iMarker++)
    {
      const auto& marker = family.GetMarker(iMarker);

      std::cout << std::setw(10) << iMarker
                << std::setw(nameColumnWidth) << marker.mNameI
                << std::setw(10) << marker.mIndexZoneI
                << std::setw(10) << lFaceName(marker.mFaceLocationI)
                << std::setw(nameColumnWidth) << marker.mNameJ
                << std::setw(10) << marker.mIndexZoneJ
                << std::setw(10) << lFaceName(marker.mFaceLocationJ)
                << std::setw(12) << (marker.mIsReversed ? "yes" : "no")
                << std::setw(24) << marker.mTranslationVector[0]
                << std::setw(24) << marker.mTranslationVector[1]
                << '\n';
    }

    lPrintSeparator(tableWidth);
  }

  std::cout << std::endl;

  // Restore the original stream formatting.
  std::cout.flags(flags);
  std::cout.precision(precision);
  std::cout.fill(fill);
}








