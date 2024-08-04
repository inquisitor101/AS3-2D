#include "vtk_structure.hpp"


//-----------------------------------------------------------------------------------
// IFileVTK member functions.
//-----------------------------------------------------------------------------------


IFileVTK::IFileVTK
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor for the interface VTK class.
	*/
{
	// Initialize the variables to write to false.
	mWriteDensity     = false;
	mWriteMomentum    = false;
	mWriteTotalEnergy = false;
	mWritePressure    = false;
	mWriteVelocity    = false;
	mWriteVorticity   = false;
	mWriteMach        = false;
	mWriteTemperature = false;
	mWriteEntropy     = false;

	// Determine what needs to be written.
	for( auto& var: config_container->GetWriteVisVar() )
	{
		switch(var)
		{
			case(EWriteVariable::DENSITY):      { mWriteDensity     = true; break; }
			case(EWriteVariable::MOMENTUM):     { mWriteMomentum    = true; break; }
			case(EWriteVariable::TOTAL_ENERGY): { mWriteTotalEnergy = true; break; }
			case(EWriteVariable::PRESSURE):     { mWritePressure    = true; break; }
			case(EWriteVariable::VELOCITY):     { mWriteVelocity    = true; break; }
			case(EWriteVariable::VORTICITY):    { mWriteVorticity   = true; break; }
			case(EWriteVariable::MACH):         { mWriteMach        = true; break; }
			case(EWriteVariable::TEMPERATURE):  { mWriteTemperature = true; break; }
			case(EWriteVariable::ENTROPY):      { mWriteEntropy     = true; break; }

			default: ERROR("Unknown variable to write.");
		}
	}

	// Reset the file number to zero.
	mFileNumber = 0;
}

//-----------------------------------------------------------------------------------

IFileVTK::~IFileVTK
(
 void
)
 /*
	* Destructor, which cleans up after the VTK class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CLegacyBinaryVTK member functions.
//-----------------------------------------------------------------------------------


CLegacyBinaryVTK::CLegacyBinaryVTK
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		IFileVTK(config_container, geometry_container)
 /*
	* Constructor for the legacy and binary VTK class.
	*/
{
	// Check for big endian. If not, the bytes must be swapped when
  // writing a binary VTK file.
  union {int i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) mBigEndian = false;
  else                  mBigEndian = true;


	// Determine the variable names to write. Note, this loop has to be similar to the one
	// in DetermineVisualizationData. Otherwise, different variables are written under 
	// different names.
	for( auto& var: config_container->GetWriteVisVar() )
	{
		// Initialize the relevant variable name, according to the correct order.
		switch(var)
		{
			case(EWriteVariable::DENSITY):  
			{ 
				mVariableNames.push_back("Density"); 
				break; 
			}

			case(EWriteVariable::MOMENTUM): 
			{
				mVariableNames.push_back("Momentum_x");
				mVariableNames.push_back("Momentum_y");
				mVariableNames.push_back("Momentum_z");
				break;
			}

			case(EWriteVariable::TOTAL_ENERGY):
			{
				mVariableNames.push_back("Energy");
				break;
			}

			case(EWriteVariable::PRESSURE):
			{
				mVariableNames.push_back("Pressure");
				break;
			}

			case(EWriteVariable::VELOCITY):
			{
				mVariableNames.push_back("Velocity_x");
				mVariableNames.push_back("Velocity_y");
				mVariableNames.push_back("Velocity_z");
				break;
			}

			case(EWriteVariable::TEMPERATURE):
			{
				mVariableNames.push_back("Temperature");
				break;
			}

			case(EWriteVariable::MACH):
			{
				mVariableNames.push_back("Mach");
				break;
			}

			case(EWriteVariable::VORTICITY):
			{
				mVariableNames.push_back("Vorticity");
				break;
			}

			case(EWriteVariable::ENTROPY):
			{
				mVariableNames.push_back("Entropy");
				break;
			}

			default: ERROR("Unknown variable name.");
		}	
	}
}

//-----------------------------------------------------------------------------------

CLegacyBinaryVTK::~CLegacyBinaryVTK
(
 void
)
 /*
	* Destructor, which cleans up after the legacy and binary VTK class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CLegacyBinaryVTK::WriteFileVTK
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Function that writes VTK data to a file using a binary and legacy format.
	*/
{
 	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Writing solution in legacy binary VTK format... ";

	// Determine number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// Create VTK filename.
	std::ostringstream fn;
	fn << config_container->GetOutputVisFilename() << "_" << mFileNumber << ".vtk";


	// Define the maximum string length for the writing.
  const int MAX_STRING_LENGTH = 255;

	// Initialize the buffer data.
	as3vector2d<float> vars_buf(nZone);
	as3vector2d<float> coor_buf(nZone);
	as3vector2d<int>   conn_buf(nZone);
	as3vector2d<int>   type_buf(nZone);

	// AS3 uses quadrilateral elements only. Based on VTK syntax, these have the below properties.
	const unsigned int quad_elem = 9;
	const unsigned int quad_npts = 4;

	// Accumulate the total number of DOFs written per zone.
	int nDOFsWritten    = 0;
	// Accumulate the total number of P1 sub-elements per zone.
	int nSubElemWritten = 0;


	// Assemble the data to write.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		// Extract the relevant zone.
		auto* zone = geometry_container->GetZoneGeometry(iZone);

		// Extract the current zone information.	
		const int nElem = static_cast<int>( zone->GetnElem() );
		const int nNode = static_cast<int>( zone->GetnNodeGrid2D() );
		const int nPoly = static_cast<int>( zone->GetnPolyGrid() );
		
		// Total DOFs and P1 sub-elements in this zone.
		const int nDOFsTot = nNode*nElem;
		const int nSubElem = nPoly*nPoly;

		// Extract number of elements in x and y in this zone.
		const int nxElem = static_cast<int>( zone->GetnxElem() );
		const int nyElem = static_cast<int>( zone->GetnyElem() );

		// Deduce number of nodes for the DOFs in 1D.
		const int nDOFsSol1D = nPoly + 1;

		// Ensure the number of elements is correct.
		if( nxElem*nyElem != nElem ) ERROR("Inconsistency in the number of elements being written.");

		// Initialize and specify the coordinates. Note, ParaView expects them in 3D.
		coor_buf[iZone].resize( 3*nDOFsTot );
		// Initialize the global connectivity matrix.
		// Note, the addition of one is for the number of points per element.
		conn_buf[iZone].resize( nElem*nSubElem*(quad_npts+1) );
		// Initialize the element types. For this code, just specify a quadrilateral.
		type_buf[iZone].resize( nElem*nSubElem, quad_elem );

		// Initialize the buffer data.
  	for(int ijElem=0; ijElem<nElem; ++ijElem)
  	{
			// Deduce the local element indices in x and y.
			//const int jElem = ijElem/nxElem;
			//const int iElem = ijElem - jElem*nxElem;

			// Deduce the current element coordinate array location.
			float *coorelem = coor_buf[iZone].data() + ijElem*3*nNode;  
			// Deduce the current element connectivity array location.
			int   *connelem = conn_buf[iZone].data() + ijElem*(quad_npts+1)*nSubElem;

			// Extract the current element coordinates.
			auto& xyelem = zone->GetElementGeometry(ijElem)->GetCoordSolDOFs();

			// Loop over every DOF on this element and store the coordinate.
			int idx = 0;
			for(int l=0; l<nNode; l++){
				coorelem[idx++] = static_cast<float>( xyelem(0,l) );
				coorelem[idx++] = static_cast<float>( xyelem(1,l) );
				coorelem[idx++] = 0.0f;
			}

			// Deduce the start of the current element's connectivity array.
			const int offset = ijElem*nNode + nDOFsWritten;
			
			// Connectivity array.
			idx = 0;
			for(int j=0; j<nPoly; j++){
				const int joff = offset + j*nDOFsSol1D;
				for(int i=0; i<nPoly; i++){

					// Starting index of each nPoly=1 sub-element.
					const int n0 = joff + i; 
					// Number of points in this element.
					connelem[idx++] = quad_npts;
					// Actual connectivity array.
					connelem[idx++] = n0;
					connelem[idx++] = n0 + 1;
					connelem[idx++] = n0 + 1 + nDOFsSol1D;
					connelem[idx++] = n0 + nDOFsSol1D;
				}
			}
		}
	
		// Update the total number of DOFs written, thus far.
		nDOFsWritten    += nDOFsTot;
		// Update the total number of P1 sub-elements written, thus far.
		nSubElemWritten += nElem*nSubElem; 

		// Allocate the size of the written data variables.
		vars_buf[iZone].resize(nDOFsTot*mVariableNames.size(), 0.0);
	}

	// Compute and store the required data for visualization. 
	// Note, this is done over all elements in all zones simulataneoously.
	DetermineVisualizationData(config_container,
			                       geometry_container, 
														 solver_container, vars_buf); 
	

	// Check if there need be any swapping, since ParaView expects data in big endian format.
	if( !mBigEndian ){
		for(unsigned short iZone=0; iZone<nZone; iZone++){
			NInputUtility::SwapBytes( coor_buf[iZone].data(), sizeof(float), coor_buf[iZone].size() );
			NInputUtility::SwapBytes( vars_buf[iZone].data(), sizeof(float), vars_buf[iZone].size() );
			NInputUtility::SwapBytes( conn_buf[iZone].data(), sizeof(int),   conn_buf[iZone].size() );
			NInputUtility::SwapBytes( type_buf[iZone].data(), sizeof(int),   type_buf[iZone].size() );
		}
	}


	// Save number of total DOFs per zone, as it is needed for writing scalar/vector data.
	int nDOFsZone[nZone];
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		const int nElem  = static_cast<int>( geometry_container->GetZoneGeometry(iZone)->GetnElem() );
		const int nNode  = static_cast<int>( geometry_container->GetZoneGeometry(iZone)->GetnNodeGrid2D() );
		nDOFsZone[iZone] = nElem*nNode; 
	} 


	/*
	 * Write the data to file.
	 */


	// Open the visualization file for binary writing.
	FILE *fh = std::fopen( fn.str().c_str(), "wb" );
	if( !fh ) ERROR("Visualization file in legacy + binary format could not be opened.");


	// Write the header of the visualization file.
  char str_buf[MAX_STRING_LENGTH];
  std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
	std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "vtk output\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "BINARY\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);


  // Write the coordinates.
  std::sprintf(str_buf, "POINTS %i float\n", nDOFsWritten);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++) 
		std::fwrite(coor_buf[iZone].data(), sizeof(float), coor_buf[iZone].size(), fh);


  // Write the connectivity data.
  std::sprintf(str_buf, "\nCELLS %i %i\n", nSubElemWritten, (quad_npts+1)*nSubElemWritten);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++)
		std::fwrite(conn_buf[iZone].data(), sizeof(int), conn_buf[iZone].size(), fh);


  // Write the element type data.
  std::sprintf(str_buf, "\nCELL_TYPES %i\n", nSubElemWritten);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  for(int iZone=0; iZone<nZone; iZone++)
		std::fwrite(type_buf[iZone].data(), sizeof(int), type_buf[iZone].size(), fh);


  // Write the ASCII line for the point data.
  std::sprintf(str_buf, "\nPOINT_DATA %i\n", nDOFsWritten);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);


	// Loop over each variables and write the output to a file.
	for(size_t iVar=0; iVar<mVariableNames.size(); iVar++)
	{
		// Copy variable name, since it will get modified.
		std::string varname = mVariableNames[iVar];

		// Check whether this is a scalar or vector variable.
    // Note that the y- and z-components of a vector are skipped.
    bool writevar = true, isvector = false;
    size_t found = varname.find("_x");
    if(found != std::string::npos) isvector = true;

    found = varname.find("_y");
    if(found != std::string::npos) writevar = false;

    found = varname.find("_z");
    if(found != std::string::npos) writevar = false;

    // Check for a vector field.
    if( isvector )
    {
      // Vector variable. Remove the trailing "_x".
      varname.erase(varname.end()-2, varname.end());

      // Write the ASCII line with the information.
      std::sprintf(str_buf, "\nVECTORS %s float\n", varname.c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the vector data.
			for(int iZone=0; iZone<nZone; iZone++)
				std::fwrite(vars_buf[iZone].data() + iVar*nDOFsZone[iZone], sizeof(float), 3*nDOFsZone[iZone], fh);

			// Skip remainder of instructions in this loop index.
			continue;
    }

		// Check for a scalar variable.
    if( writevar )
    {
      // Scalar variable. Write the ASCII lines with the information.
      std::sprintf(str_buf, "\nSCALARS %s float 1\n", varname.c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      std::sprintf(str_buf, "LOOKUP_TABLE default\n");
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the scalar data.
			for(int iZone=0; iZone<nZone; iZone++)
				std::fwrite(vars_buf[iZone].data() + iVar*nDOFsZone[iZone], sizeof(float), nDOFsZone[iZone], fh);
    
			// Skip remainder of instructions in this loop index.
			continue;
		}
  }


	// Increment the file counter index.
	mFileNumber++;


  // Close the file again.
  std::fclose(fh);

	// Report: all is complete!
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CLegacyBinaryVTK::DetermineVisualizationData
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 as3vector2d<float>                    &vars_buf
)
 /*
	* Function that computes the required data for visualization in binary format.
	*/
{
	// Abbreviation involving gamma.
	const as3double gm1 = GAMMA_MINUS_ONE;

	// Extract total number of zones.
	const unsigned short nZone = config_container->GetnZone();

	for(size_t iZone=0; iZone<nZone; iZone++)
	{
		// Extract current grid zone.
		auto* zone   = geometry_container->GetZoneGeometry(iZone);
		// Extract current solver.
		auto* solver = solver_container[iZone].get();

		// Extract the total number of elements in this zone.
		const int nElem = static_cast<int>( zone->GetnElem() );

		// Extract the current zone information.
		const int nNode    = static_cast<int>( zone->GetnNodeGrid2D() );
		// Total DOFs and sub-elements in this zone.
		const int nDOFsTot = static_cast<int>( nNode*nElem );
		// Extract number of elements in x-direction in this zone.
		const int nxElem   = static_cast<int>( zone->GetnxElem() );


		// Loop over every element and process the data for visualization.
		for(size_t ijElem=0; ijElem<nElem; ijElem++)
		{
			// Extract solution in this element.
			auto& sol = solver->GetPhysicalElement(ijElem)->mSol2D;

			// Allocate data for the primitive variables on the solution DOFs for this element.
			// These contain: 1/rho, u, v, p. 
			// Note, the inverse of rho is computed for convenience.
			CMatrixAS3<as3double> primvar( 4, nNode );

			// Loop over the element DOFs and compute the primitive variables needed.
			for(size_t l=0; l<nNode; l++)
			{
				const as3double ovrho =   1.0/sol(0,l);                       
				const as3double u     = ovrho*sol(1,l);                     
				const as3double v     = ovrho*sol(2,l);                     
				const as3double p     = gm1*( sol(3,l) 
						                  -   0.5*sol(0,l)*(u*u + v*v) );
				
				// Store the values.
				primvar(0,l) = ovrho;
				primvar(1,l) = u;
				primvar(2,l) = v;
				primvar(3,l) = p;
			}

			// Counter for the number of variables.
			int iVar = 0;

			// Determine which variables to compute.
			for( auto& var: config_container->GetWriteVisVar() )
			{

				// Check which variable to compute.
				switch(var)
				{

					case(EWriteVariable::DENSITY): 
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
				
						// Compute and store the density in each DOF.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( sol(0,l) ); 
						}

						// Update variable counter for a scalar.
						iVar++;

						break;
					}

					case(EWriteVariable::MOMENTUM):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode*3;
	
						// Compute and store the momentum variables.
						for(size_t l=0; l<nNode; l++)
						{
							const int l3 = 3*l;
							buf[l3  ] = static_cast<float>( sol(1,l) );
							buf[l3+1] = static_cast<float>( sol(2,l) );
							buf[l3+2] = static_cast<float>(    0.0   );
						}
	
						// Update variable counter for a 3D vector.
						iVar += 3;
						
						break;
					}

					case(EWriteVariable::TOTAL_ENERGY):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
						// Compute and store the energy variable.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( sol(3,l) );
						}

						// Update variable counter for a scalar.
						iVar++;
						
						break;
					}

					case(EWriteVariable::PRESSURE):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
						
						// Compute and store the pressure variable.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( primvar(3,l) );
						}

						// Update variable counter for a scalar.
						iVar++;

						break;
					}

					case(EWriteVariable::VELOCITY):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode*3;
	
						// Compute and store the velocity variables.
						for(size_t l=0; l<nNode; l++)
						{
							const int l3 = 3*l;
							const as3double ovrho = primvar(0,l);
							buf[l3  ] = static_cast<float>( ovrho*sol(1,l) );
							buf[l3+1] = static_cast<float>( ovrho*sol(2,l) );
							buf[l3+2] = static_cast<float>(       0.0      );
						}

						// Update variable counter for a 3D vector.
						iVar += 3;
	
						break;
					}

					case(EWriteVariable::TEMPERATURE):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
						// Abbreviation for 1/R.
						const as3double ovR = 1.0/GAS_CONSTANT;
	
						// Compute and store the temperature variable.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( ovR*primvar(0,l)*primvar(3,l) );
						}

						// Update variable counter for a scalar.
						iVar++;

						break;
					}

					case(EWriteVariable::MACH):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
						// Compute and store the Mach number variable.
						for(size_t l=0; l<nNode; l++)
						{
							const as3double a2    = GAMMA*primvar(0,l)*primvar(3,l);
							const as3double umag2 = primvar(1,l)*primvar(1,l) 
								                    + primvar(2,l)*primvar(2,l); 
							buf[l] = static_cast<float>( std::sqrt( umag2/a2 ) ); 
						}

						// Update variable counter for a scalar.
						iVar++;

						break;
					}

					case(EWriteVariable::ENTROPY):
					{
						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;
	
						// Compute and store the Mach number variable.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( std::log( primvar(3,l)*std::pow( primvar(0,l), GAMMA ) ) );
						}

						// Update variable counter for a scalar.
						iVar++;

						break;
					}

					case(EWriteVariable::VORTICITY):
					{
						// Allocate memory for the velocity gradient.
						CMatrixAS3<as3double> dUDx( 2, nNode );
						CMatrixAS3<as3double> dUDy( 2, nNode );

						// Extract number of solution points in 1D.
						const size_t nSol1D = solver->GetStandardElement()->GetnSol1D();

						// Create an identity matrix, needed for the Lagrange interpolation over the DOFs.
						CMatrixAS3<as3double> identity = NLinearAlgebra::CreateIdentityMatrix<as3double>( nSol1D ); 

						// Extract the transposed differentiation matrix in 1D over the solution points.
						CMatrixAS3<as3double> derLagrangeSol1DTrans = solver->GetStandardElement()->GetDerLagrangeSol1DTrans();

						// Compute the gradient of the velocity.
						solver->GetTensorProduct()->CustomVolume(nSol1D, 2, nSol1D,
								                                     identity.data(), derLagrangeSol1DTrans.data(),
																                     &primvar[nNode], nullptr,
																                     dUDx.data(), dUDy.data());

						// Extract the metrics at the volume solution points.
						auto& metrics = solver->GetPhysicalElement(ijElem)->mMetricSol2D;

						// Convert the derivatives from parametric space to Cartesian space.
						for(size_t l=0; l<nNode; l++)
						{
							// Parametric derivatives.
							const as3double dudr = dUDx(0,l);
							const as3double dvdr = dUDx(1,l);
							const as3double duds = dUDy(0,l);
							const as3double dvds = dUDy(1,l);

							// Metrics of transformation.
							const as3double drdx = metrics(1,l);
							const as3double drdy = metrics(2,l);
							const as3double dsdx = metrics(3,l);
							const as3double dsdy = metrics(4,l);

							// Compute the Cartesian derivatives w.r.t. x.
							//dUDx(0,l) = dudr*drdx + duds*dsdx; // dudx (not needed here)
							dUDx(1,l) = dvdr*drdx + dvds*dsdx; // dvdx
							
							// Compute the Cartesian derivatives w.r.t. y.
							dUDy(0,l) = dudr*drdy + duds*dsdy; // dudy
							//dUDy(1,l) = dvdr*drdy + dvds*dsdy; // dvdy (not needed here)
						}

						// Set the pointer to the correct location.
						float *buf = vars_buf[iZone].data() + iVar*nDOFsTot + ijElem*nNode;

						// Compute and store the Vorticity variable.
						for(size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<float>( dUDx(1,l) - dUDy(0,l) ); 
						}

						// Update variable counter for a scalar.
						iVar += 1;
	
						break;
					}


					default: ERROR("Unknown variable detected.");
				}
			}
			
		}
	}
}


