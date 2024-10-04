#include "iteration_structure.hpp"


//-----------------------------------------------------------------------------------
// CIteration member functions.
//-----------------------------------------------------------------------------------


CIteration::CIteration
(
 CConfig                               *config_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Constructor for the iteration class.
	*/
{
	// Check the relevant number of data entries required in the work array.
	for(unsigned short iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
		// For now, take the maximum number of items.
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
				const size_t nItem2D = 3;  // volume  terms needed.
				const size_t nItem1D = 2;  // surface terms needed.
				
				const size_t nVar    = solver_container[iZone]->GetnVar();
				const size_t nInt1D  = solver_container[iZone]->GetStandardElement()->GetnInt1D();
				const size_t nInt2D  = solver_container[iZone]->GetStandardElement()->GetnInt2D();

				// Compute the total number of required volume and surface terms in the work array.
				const size_t nVol  = nItem2D*nInt2D*nVar;
				const size_t nSurf = nItem1D*nInt1D*nVar;

				// Take the maximum storage between the volume and surface terms.
				const size_t nData = std::max( nVol, nSurf );

				// Take whichever is the max possible storage across all zones (can be inefficient).
				mNWorkData = std::max( mNWorkData, nData );
				
				break;
			}

			default: ERROR("Unknown solver type.");
		}
	}
}

//-----------------------------------------------------------------------------------

CIteration::~CIteration
(
 void
)
 /*
	* Destructor, which cleans up after the iteration class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PreProcessIteration
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime 
)
 /*
	* Function that preprocesses the solution, before sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PostProcessIteration
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime 
)
 /*
	* Function that postprocesses the solution, after sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::ComputeResidual
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime
)
 /*
	* Function that computes the residual in all zones.
	*/
{
	// Get the total number of elements in all zones.
	const size_t nElemTotal = openmp_container->GetnElemTotal();
	// Get the total number of internal i-faces in all zones.
	const size_t nInternalIFace = openmp_container->GetnInternIFace();
	// Get the total number of internal j-faces in all zones.
	const size_t nInternalJFace = openmp_container->GetnInternJFace();

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Deduce the current element's zone and index.
		const unsigned short iZone = openmp_container->GetIndexVolume(i)->mZone;
		const unsigned int   iElem = openmp_container->GetIndexVolume(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant grid.
		auto* grid    = geometry_container->GetZoneGeometry(iZone);

		// Compute the volume terms on this element. Note, this step also initializes the residual.
		solver->ComputeVolumeResidual(grid, workarray, localtime, iElem); 
	}


#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nInternalIFace; i++)
	{
		// Deduce the right element's zone and index.
		const unsigned short iZone = openmp_container->GetInternIFace(i)->mZone;
		const unsigned int   iElem = openmp_container->GetInternIFace(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant grid.
		auto* grid    = geometry_container->GetZoneGeometry(iZone);
		// Extract the left residual.
		auto& resL    = openmp_container->GetResIMin(i);

		// Reset the left residual.
		for(size_t l=0; l<resL.size(); l++) resL[l] = C_ZERO;

		// Compute all surface terms in the i-direction.
		solver->ComputeSurfaceResidualIDir(grid, workarray, localtime, iElem, resL);
	}


#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nInternalJFace; i++)
	{
		// Deduce the top element's zone and index.
		const unsigned short iZone = openmp_container->GetInternJFace(i)->mZone;
		const unsigned int   iElem = openmp_container->GetInternJFace(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant grid.
		auto* grid    = geometry_container->GetZoneGeometry(iZone);
		// Extract the bottom residual.
		auto& resB    = openmp_container->GetResJMin(i);

		// Reset the bottom residual.
		for(size_t l=0; l<resB.size(); l++) resB[l] = C_ZERO;

		// Compute all surface terms in the j-direction.
		solver->ComputeSurfaceResidualJDir(grid, workarray, localtime, iElem, resB);
	}


#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nInternalIFace; i++)
	{
		// Deduce the right element's zone and index.
		const unsigned short iZone = openmp_container->GetInternIFace(i)->mZone;
		const unsigned int   IR    = openmp_container->GetInternIFace(i)->mElem;

		// Deduce the left element's index.
		const unsigned int IL = IR-1;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the temporary left residual.
		auto& tmpL    = openmp_container->GetResIMin(i);
		// Extract the actual left residual.
		auto& resL    = solver->GetPhysicalElement(IL)->mRes2D;

		// Accumulate the left residual.
		for(size_t l=0; l<resL.size(); l++) resL[l] += tmpL[l];
	}


#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nInternalJFace; i++)
	{
		// Deduce the top element's zone and index.
		const unsigned short iZone = openmp_container->GetInternJFace(i)->mZone;
		const unsigned int   IT    = openmp_container->GetInternJFace(i)->mElem;

		// Deduce the bottom element's index.
		const size_t IB = IT - geometry_container->GetZoneGeometry(iZone)->GetnxElem();

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the temporary bottom residual.
		auto& tmpB    = openmp_container->GetResJMin(i);
		// Extract the actual bottom residual.
		auto& resB    = solver->GetPhysicalElement(IB)->mRes2D;

		// Accumulate the bottom residual.
		for(size_t l=0; l<resB.size(); l++) resB[l] += tmpB[l];
	}


	// Loop over the interface boundaries, if need be.
	for( auto& interface: interface_container )
	{
		interface->ComputeInterfaceResidual(solver_container, workarray);
	}



	// Multiply by the inverse mass matrix.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Deduce the current element's zone and index.
		const unsigned short iZone = openmp_container->GetIndexVolume(i)->mZone;
		const unsigned int   iElem = openmp_container->GetIndexVolume(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant physical element.
		auto* element = solver->GetPhysicalElement(iElem);

		// Get the inverse of the mass matrix.
		auto& m = element->mInvMassMatrix;
		// Get the residual on this element.
		auto& r = element->mRes2D;

		// Perform a matrix-matrix multiplication to obtain the residual.
		NLinearAlgebra::MatrixVectorTransMult(m, r, r);
	}
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container, 
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 as3double                                 localtime 
)
 /*
	* Function that performs a grid sweep over all the zones. 
	*/
{
	// Initialize a work array, to avoid multiple memory allocations.
	// Note, during parallelization, this needs to be allocated inside 
	// the (shared memory) parallel region -- not here.
	CPoolMatrixAS3<as3double> workarray(mNWorkData);


	// Check for any preprocessing steps.
	PreProcessIteration(config_container,
			                geometry_container,
											openmp_container,
											solver_container,
											interface_container,
											workarray,
											localtime);


	// Compute the residual over all zones.
	ComputeResidual(config_container,
			            geometry_container,
									openmp_container,
									solver_container,
									interface_container,
									workarray,
									localtime);


	// Check for any postprocessing steps.
	PostProcessIteration(config_container,
			                 geometry_container,
											 openmp_container,
											 solver_container,
											 interface_container,
											 workarray,
											 localtime);
}


