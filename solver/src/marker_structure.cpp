#include "marker_structure.hpp"

#include "geometry_structure.hpp"


//-----------------------------------------------------------------------------------
// CMarker member functions.
//-----------------------------------------------------------------------------------


CMarker::CMarker
(
 unsigned short             zone,
 ETypeBC                    type,
 std::string                name,
 as3vector1d<EFaceLocation> face,
 as3vector1d<unsigned int>  mark
)
	:
		mZoneID(zone),
		mTypeBC(type),
		mNameMarker(name)
 /*
	* Constructor for the marker geometry, which contains a single marker region.
	*/
{
	// Ensure the faces and index markes are of the same size.
	if( face.size() != mark.size() ) ERROR("Mismatch in the size of face and element indices.");

	// Reserve memory for the element face properties.
	mElementFaces.reserve( face.size() );
	
	// Initialize the member variables.
	for(size_t i=0; i<face.size(); i++)
	{
		mElementFaces.emplace_back( mark[i], face[i] );
	}

  // Ensure the element markers have the same face location, as assumed is the case.
  mFaceLocation = mElementFaces[0].mFace;
  for(const auto& face_marker : mElementFaces)
  {
    if(face_marker.mFace != mFaceLocation) ERROR(mNameMarker + " must have the same face location on all its elements.");
  }

  // Deduce the face type, based on its location.
  switch( mFaceLocation )
  {
    case( EFaceLocation::IMIN ) : case( EFaceLocation::IMAX ):  mFaceType = ETypeFace::IFACE; break;
    case( EFaceLocation::JMIN ) : case( EFaceLocation::JMAX ):  mFaceType = ETypeFace::JFACE; break;
    default: ERROR("Cannot deduce face type from face location.");
  }
}

//-----------------------------------------------------------------------------------

CMarker::~CMarker
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}

//-----------------------------------------------------------------------------------
// CPeriodicFamilyMarker member functions.
//-----------------------------------------------------------------------------------

CPeriodicFamilyMarker::CPeriodicFamilyMarker
(
 const CMultizoneGeometry       *multizone_geometry_container,
 const CExternalFamilyMarker    &iexternal_marker,
 const CExternalFamilyMarker    &jexternal_marker,
 const std::array<as3double, 2> &translation
)
 /*
  *
  */
{
  // Consistency check.
  if( iexternal_marker.GetnMarkers() != jexternal_marker.GetnMarkers() )
  {
    ERROR("Periodic marker must have the same number of markers on both sides.");
  }
 
  // Extract the boundary names.
  const auto iName = iexternal_marker.GetName();
  const auto jName = jexternal_marker.GetName();

  // Deduce the number of markers.
  const size_t nMarkers = iexternal_marker.GetnMarkers();

  // Reserve the necessary number of markers.
  mMarkers.reserve( nMarkers );

  // Relative tolerance value.
  const as3double rtol = static_cast<as3double>( 1.0e-8 );

  // Absolute tolerance value, in grid coordinate units.
  const as3double atol = static_cast<as3double>( 1.0e-10 );

  // Compare a translated i coordinate with its matching j coordinate.
  auto lCoordinatesMatch = [rtol, atol](as3double xi, as3double xj)
  {
    const as3double difference = std::abs( xi - xj );
    const as3double tolerance = atol + rtol * std::max( std::abs(xi), std::abs(xj) );

    return difference <= tolerance;
  };

  // Track jmarkers already matched with a previous imarker on this boundary.
  as3vector1d<bool> jmatched( nMarkers, false );

  for( const auto& imarker : iexternal_marker.GetMarkers() )
  {
    const auto iZone = imarker.mIndexZone;
    const auto iFace = imarker.mFaceLocation;

    const auto& igrid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

    const auto ivec_ielem = igrid->ComputeSurfaceElementIndices( iFace );
    const auto ivec_inode = igrid->ComputeSurfaceNodeIndices( iFace );

    bool found = false;

    for(std::size_t jMarker=0; jMarker<nMarkers; jMarker++)
    {
      // Avoid repetitive checks on the same jmarker, in case it was
      // already matched with a previous imarker on this boundary.
      if( jmatched[jMarker] ) continue;

      // Obtain the current jmarker.
      const auto& jmarker = jexternal_marker.GetMarker(jMarker);

      const auto jZone = jmarker.mIndexZone;
      const auto jFace = jmarker.mFaceLocation;
    
      const auto& jgrid = multizone_geometry_container->GetSinglezoneGeometry(jZone);

      const auto jvec_ielem = jgrid->ComputeSurfaceElementIndices( jFace );
      const auto jvec_inode = jgrid->ComputeSurfaceNodeIndices( jFace );

      if( ivec_ielem.size() != jvec_ielem.size() ) continue;
      
      // Prevent accepting an empty surface or accessing an empty nodal vector.
      if( ivec_ielem.empty() || ivec_inode.empty() || jvec_inode.empty() )
      {
        ERROR("Cannot match periodic boundaries with empty surface indices.");
      }
      
      // Check every element on this marker in the specified orientation.
      // The i indexing remains unchanged; only the j indexing is reversed.
      auto lSurfaceMatches = [&](bool reversed)
      {
        // Abbreviations for the nodal indices.
        const std::size_t in0 = ivec_inode.front();
        const std::size_t in1 = ivec_inode.back();

        // Swap the j endpoints when checking the reversed orientation.
        const std::size_t jn0 = reversed ? jvec_inode.back()  : jvec_inode.front();
        const std::size_t jn1 = reversed ? jvec_inode.front() : jvec_inode.back();

        // Loop over every element on this marker.
        for(size_t k=0; k<ivec_ielem.size(); k++)
        {
          // Abbreviations for the element indices.
          const std::size_t ie = ivec_ielem[k];
          const std::size_t je =
            jvec_ielem[ reversed ? jvec_ielem.size()-1-k : k ];

          // Extract the entire element coordinates.
          const auto& icoor = igrid->GetElementGeometry( ie )->GetCoordSolDOFs();
          const auto& jcoor = jgrid->GetElementGeometry( je )->GetCoordSolDOFs();

          // Extract the endpoint coordinates for the current element edge, apply the translation.
          const as3double xi0 = icoor( 0, in0 ) + translation[0], xi1 = icoor( 0, in1 ) + translation[0];
          const as3double yi0 = icoor( 1, in0 ) + translation[1], yi1 = icoor( 1, in1 ) + translation[1];

          // Extract coordinates for the matching marker.
          // For the reversed case, jn0 and jn1 already select the opposite endpoints.
          const as3double xj0 = jcoor( 0, jn0 ), xj1 = jcoor( 0, jn1 );
          const as3double yj0 = jcoor( 1, jn0 ), yj1 = jcoor( 1, jn1 );

          // Reject this candidate if either endpoint does not match.
          if( !lCoordinatesMatch(xi0, xj0) || !lCoordinatesMatch(yi0, yj0) ||
              !lCoordinatesMatch(xi1, xj1) || !lCoordinatesMatch(yi1, yj1) )
          {
            return false;
          }
        }

        // Every element passed in this orientation.
        return true;
      };

      // First, check the forward case.
      const bool forward_match = lSurfaceMatches(false);

      // If the forward detection failed, check the entire marker in reverse.
      const bool reversed_match = !forward_match && lSurfaceMatches(true);

      // Neither orientation matched: check the next available jmarker.
      if( !forward_match && !reversed_match ) continue;

      // Every element passed in one orientation: this jmarker matches the current imarker.
      jmatched[jMarker] = true; found = true;

      // Initialize the periodic marker with the detected orientation.
      mMarkers.push_back( CPeriodicMarker{iZone, jZone,
                                          iFace, jFace,
                                          iName, jName,
                                          reversed_match, translation} );

      // Stop searching jmarkers for this imarker.
      break;
    }

    if( !found ) 
    {
      ERROR("Periodic boundaries: " + iName + ", " + jName + " do not match.");
    }
  }





  // DEBUGGING

  //// Consistency check.
  //if( iexternal_marker.GetnMarkers() != jexternal_marker.GetnMarkers() )
  //{
  //  ERROR("Periodic marker must have the same number of markers on both sides.");
  //}
 
  //// Extract the boundary names.
  //const auto iName = iexternal_marker.GetName();
  //const auto jName = jexternal_marker.GetName();

  //// Deduce the number of markers.
  //const size_t nMarkers = iexternal_marker.GetnMarkers();

  //// Reserve the necessary number of markers.
  //mMarkers.reserve( nMarkers );

  //// Relative tolerance value.
  //const as3double rtol = static_cast<as3double>( 1.0e-8 );

  //// Absolute tolerance value, in grid coordinate units.
  //const as3double atol = static_cast<as3double>( 1.0e-10 );

  //// Compare a translated i coordinate with its matching j coordinate.
  //auto lCoordinatesMatch = [rtol, atol](as3double xi, as3double xj)
  //{
  //  const as3double difference = std::abs( xi - xj );
  //  const as3double tolerance = atol + rtol * std::max( std::abs(xi), std::abs(xj) );

  //  return difference <= tolerance;
  //};

  //// Return a readable face location for diagnostic messages.
  //auto lFaceName = [](EFaceLocation face) -> const char*
  //{
  //  switch( face )
  //  {
  //    case( EFaceLocation::IMIN ): return "imin";
  //    case( EFaceLocation::IMAX ): return "imax";
  //    case( EFaceLocation::JMIN ): return "jmin";
  //    case( EFaceLocation::JMAX ): return "jmax";
  //    default: return "unknown";
  //  }
  //};

  //// Display the values used in one coordinate comparison.
  //auto lPrintCoordinate = [rtol, atol](const char *coordinate,
  //                                   as3double   xi_raw,
  //                                   as3double   xi,
  //                                   as3double   xj)
  //{
  //  const as3double difference = std::abs( xi - xj );
  //  const as3double tolerance =
  //    atol + rtol * std::max( std::abs(xi), std::abs(xj) );

  //  // Preserve the original stream formatting.
  //  const auto flags = std::cout.flags();
  //  const auto precision = std::cout.precision();
  //  const auto fill = std::cout.fill();

  //  std::cout << std::scientific << std::setprecision(16)
  //            << std::right << std::setfill(' ');

  //  std::cout << std::setw(6)  << coordinate
  //            << std::setw(26) << xi_raw
  //            << std::setw(26) << xi
  //            << std::setw(26) << xj
  //            << std::setw(26) << difference
  //            << std::setw(26) << tolerance
  //            << "  " << (difference <= tolerance ? "PASS" : "FAIL")
  //            << std::endl;

  //  // Restore the original stream formatting.
  //  std::cout.flags(flags);
  //  std::cout.precision(precision);
  //  std::cout.fill(fill);
  //};

  //// Track jmarkers already matched with a previous imarker on this boundary.
  //as3vector1d<bool> jmatched( nMarkers, false );

  //for( const auto& imarker : iexternal_marker.GetMarkers() )
  //{
  //  const auto iZone = imarker.mIndexZone;
  //  const auto iFace = imarker.mFaceLocation;

  //  if( iZone == 9 ) continue;

  //  const auto& igrid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

  //  const auto ivec_ielem = igrid->ComputeSurfaceElementIndices( iFace );
  //  const auto ivec_inode = igrid->ComputeSurfaceNodeIndices( iFace );

  //  bool found = false;

  //  for(std::size_t jMarker=0; jMarker<nMarkers; jMarker++)
  //  {
  //    std::cout << "... searching in iZone: " << iZone << std::endl;

  //    // Avoid repetitive checks on the same jmarker, in case it was
  //    // already matched with a previous imarker on this boundary.
  //    if( jmatched[jMarker] )
  //    {
  //      std::cout << "SKIPPED: J marker " << jMarker
  //                << " already matched." << std::endl;
  //      continue;
  //    }

  //    // Obtain the current jmarker.
  //    const auto& jmarker = jexternal_marker.GetMarker(jMarker);

  //    const auto jZone = jmarker.mIndexZone;
  //    const auto jFace = jmarker.mFaceLocation;
  //  
  //    const auto& jgrid = multizone_geometry_container->GetSinglezoneGeometry(jZone);

  //    const auto jvec_ielem = jgrid->ComputeSurfaceElementIndices( jFace );
  //    const auto jvec_inode = jgrid->ComputeSurfaceNodeIndices( jFace );

  //    if( ivec_ielem.size() != jvec_ielem.size() )
  //    {
  //      std::cout << "SKIPPED: I zone " << iZone
  //                << ", J zone " << jZone
  //                << ", surface element counts: "
  //                << ivec_ielem.size() << " versus "
  //                << jvec_ielem.size() << std::endl;
  //      continue;
  //    }

  //    // Prevent accepting an empty surface or accessing an empty nodal vector.
  //    if( ivec_ielem.empty() || ivec_inode.empty() || jvec_inode.empty() )
  //    {
  //      ERROR("Cannot match periodic boundaries with empty surface indices.");
  //    }
  //    
  //    // Check every element on this marker in the specified orientation.
  //    // The i indexing remains unchanged; only the j indexing is reversed.
  //    auto lSurfaceMatches = [&](bool reversed)
  //    {
  //      // Abbreviations for the nodal indices.
  //      const std::size_t in0 = ivec_inode.front();
  //      const std::size_t in1 = ivec_inode.back();

  //      // Swap the j endpoints when checking the reversed orientation.
  //      const std::size_t jn0 = reversed ? jvec_inode.back()  : jvec_inode.front();
  //      const std::size_t jn1 = reversed ? jvec_inode.front() : jvec_inode.back();

  //      // Loop over every element on this marker.
  //      for(size_t k=0; k<ivec_ielem.size(); k++)
  //      {
  //        // Abbreviations for the element indices.
  //        const std::size_t ie = ivec_ielem[k];
  //        const std::size_t je =
  //          jvec_ielem[ reversed ? jvec_ielem.size()-1-k : k ];

  //        // Extract the entire element coordinates.
  //        const auto& icoor = igrid->GetElementGeometry( ie )->GetCoordSolDOFs();
  //        const auto& jcoor = jgrid->GetElementGeometry( je )->GetCoordSolDOFs();

  //        // Extract the endpoint coordinates for the current element edge, apply the translation.
  //        const as3double xi0 = icoor( 0, in0 ) + translation[0], xi1 = icoor( 0, in1 ) + translation[0];
  //        const as3double yi0 = icoor( 1, in0 ) + translation[1], yi1 = icoor( 1, in1 ) + translation[1];

  //        // Extract coordinates for the matching marker.
  //        // For the reversed case, jn0 and jn1 already select the opposite endpoints.
  //        const as3double xj0 = jcoor( 0, jn0 ), xj1 = jcoor( 0, jn1 );
  //        const as3double yj0 = jcoor( 1, jn0 ), yj1 = jcoor( 1, jn1 );

  //        if( iZone == 9 && jZone == 4 )
  //        {
  //          std::cout << std::scientific << std::setprecision(6) << std::showpos;
  //          std::cout << "  "   << k << ", is_reversed: " << reversed << "\n"
  //                    << "    Pi0: " << xi0 << ", " << yi0 << "\n"
  //                    << "    Pi1: " << xi1 << ", " << yi1 << "\n" 
  //                    << "    Pj0: " << xj0 << ", " << yj0 << "\n"
  //                    << "    Pj1: " << xj1 << ", " << yj1 << std::endl; 
  //        }


  //        // Reject this candidate if either endpoint does not match.
  //        if( !lCoordinatesMatch(xi0, xj0) || !lCoordinatesMatch(yi0, yj0) ||
  //            !lCoordinatesMatch(xi1, xj1) || !lCoordinatesMatch(yi1, yj1) )
  //        {
  //          //std::cout << "\nREJECTED: "
  //          //          << (reversed ? "reversed" : "forward")
  //          //          << ", I zone " << iZone << " (" << lFaceName(iFace) << ")"
  //          //          << ", J zone " << jZone << " (" << lFaceName(jFace) << ")"
  //          //          << "\n  k=" << k
  //          //          << ", iElem=" << ie << ", jElem=" << je
  //          //          << "\n  iNodes=(" << in0 << ", " << in1 << ")"
  //          //          << ", jNodes=(" << jn0 << ", " << jn1 << ")"
  //          //          << "\n  iMatrix=(" << icoor.row() << ", " << icoor.col() << ")"
  //          //          << ", jMatrix=(" << jcoor.row() << ", " << jcoor.col() << ")"
  //          //          << '\n';

  //          //std::cout << std::setw(6)  << "Coord"
  //          //          << std::setw(26) << "I original"
  //          //          << std::setw(26) << "I translated"
  //          //          << std::setw(26) << "J matching"
  //          //          << std::setw(26) << "Difference"
  //          //          << std::setw(26) << "Tolerance"
  //          //          << "  Result\n";

  //          //// Display both endpoints, including comparisons that passed.
  //          //lPrintCoordinate("x0", icoor(0, in0), xi0, xj0);
  //          //lPrintCoordinate("y0", icoor(1, in0), yi0, yj0);
  //          //lPrintCoordinate("x1", icoor(0, in1), xi1, xj1);
  //          //lPrintCoordinate("y1", icoor(1, in1), yi1, yj1);

  //          //if( iZone == 9 && jZone == 4 ) continue;
  //          return false;
  //        }
  //      }

  //      // Every element passed in this orientation.
  //      return true;
  //    };

  //    // First, check the forward case.
  //    const bool forward_match = lSurfaceMatches(false);

  //    // If the forward detection failed, check the entire marker in reverse.
  //    const bool reversed_match = !forward_match && lSurfaceMatches(true);

  //    // Neither orientation matched: check the next available jmarker.
  //    if( !forward_match && !reversed_match ) continue;


  //    std::cout << "ACCEPTED: I zone " << iZone
  //              << " (" << lFaceName(iFace) << ")"
  //              << " -> J zone " << jZone
  //              << " (" << lFaceName(jFace) << ")"
  //              << ", orientation="
  //              << (reversed_match ? "reversed" : "forward")
  //              << std::endl;

  //    // Every element passed in one orientation: this jmarker matches the current imarker.
  //    jmatched[jMarker] = true; found = true;
  // 
  //    // Initialize the periodic marker with the detected orientation.
  //    mMarkers.push_back( CPeriodicMarker{iZone, jZone,
  //                                        iFace, jFace,
  //                                        iName, jName,
  //                                        reversed_match, translation} );

  //    // Stop searching jmarkers for this imarker.
  //    break;
  //  }

  //  if( !found ) 
  //  {
  //    ERROR("Periodic boundaries: " + iName + ", " + jName +
  //          " do not match for I zone " + std::to_string(iZone) +
  //          " (" + lFaceName(iFace) + ").");
  //  }
  //}



  // EXPECTATION (periodic zone matches):
  //  iZone:  0 = jZone: 2, (forward)
  //  iZone: 10 = jZone: 3, (??) -- print coordinates and see?
  //  iZone:  9 = jZone: 4  (??) -- print coordinates and see?


}





