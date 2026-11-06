#include "visualization/IFileVTK.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <type_traits>

#include "CConfig.hpp"
#include "compute_kernels/TensorProducts.hpp"
#include "data/CMatrixAS3.hpp"
#include "element/CPhysicalElement.hpp"
#include "element/CStandardElement.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "linear_algebra/NLinearAlgebra.hpp"
#include "solver/CMultizoneSolver.hpp"
#include "solver/ISinglezoneSolver.hpp"


//-----------------------------------------------------------------------------------
// IFileVTK member functions.
//-----------------------------------------------------------------------------------


IFileVTK::IFileVTK
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Constructor that selects the visualization fields and computes their buffer offsets.
	*/
{
	for(const auto variable : config_container->GetWriteVisVar())
	{
		CPointField field{variable, nullptr, mNPointComponents, 1};

		switch( variable )
		{
			case(EWriteVariable::DENSITY):      { field.mName = "Density";   break; }
			case(EWriteVariable::MOMENTUM):     { field.mName = "Momentum";  field.mNComponent = 3; break; }
			case(EWriteVariable::TOTAL_ENERGY): { field.mName = "Energy";    break; }
			case(EWriteVariable::PRESSURE):     { field.mName = "Pressure";  break; }
			case(EWriteVariable::VELOCITY):     { field.mName = "Velocity";  field.mNComponent = 3; break; }
			case(EWriteVariable::VORTICITY):    { field.mName = "Vorticity"; break; }
			case(EWriteVariable::MACH):         { field.mName = "Mach";      break; }
			case(EWriteVariable::TEMPERATURE):  { field.mName = "Temperature"; break; }
			case(EWriteVariable::ENTROPY):      { field.mName = "Entropy";   break; }

			default: ERROR("Unknown variable to write.");
		}

		mNPointComponents += field.mNComponent;
		mPointFields.push_back(field);
	}
}

//-----------------------------------------------------------------------------------

IFileVTK::~IFileVTK(void) = default;

//-----------------------------------------------------------------------------------

void IFileVTK::WriteSeriesIndex
(
 const CConfig      *config_container,
 const std::string &filename,
 as3double          time
)
 /*
	* Function that records physical times in a PVD collection or legacy VTK file series.
	*/
{
	if( !std::isfinite(time) ) ERROR("Cannot write a non-finite visualization time.");
	mSnapshots.push_back({std::filesystem::path(filename).filename().string(), time});

	const bool isXML = config_container->GetOutputVisFormat() == EVisualFormat::VTK_XML_BINARY;
	const std::string index = std::string(config_container->GetOutputVisFilename())
	                        + (isXML ? ".pvd" : ".vtk.series");
	const std::string temporary = index + ".tmp";

	// Replace the previous index only after the updated file has been written successfully.
	std::ofstream file(temporary);
	if( !file ) ERROR("Could not open visualization time-series index: " + temporary);
	file.imbue(std::locale::classic());
	file << std::setprecision(std::numeric_limits<as3double>::max_digits10);

	if( isXML )
	{
		file << "<?xml version=\"1.0\"?>\n"
		     << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
		     << "  <Collection>\n";
		for(const auto& snapshot : mSnapshots)
		{
			file << "    <DataSet timestep=\"" << snapshot.mTime << "\" group=\"\" part=\"0\" file=\"";
			for(char c : snapshot.mName)
			{
				switch( c )
				{
					case '&': file << "&amp;";  break;
					case '<': file << "&lt;";   break;
					case '"': file << "&quot;"; break;
					default:  file << c;       break;
				}
			}
			file << "\"/>\n";
		}
		file << "  </Collection>\n</VTKFile>\n";
	}
	else
	{
		file << "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n";
		for(std::size_t i=0; i<mSnapshots.size(); i++)
		{
			if( i != 0 ) file << ",\n";
			const auto& snapshot = mSnapshots[i];
			file << "    { \"name\": " << std::quoted(snapshot.mName)
			     << ", \"time\": " << snapshot.mTime << " }";
		}
		file << "\n  ]\n}\n";
	}

	file.close();
	if( !file ) ERROR("Could not write visualization time-series index: " + temporary);
	if( std::rename(temporary.c_str(), index.c_str()) != 0 )
		ERROR("Could not replace visualization time-series index: " + index);
}

//-----------------------------------------------------------------------------------

template<typename T>
void IFileVTK::DetermineVisualizationData
(
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container,
 as3vector2d<T>           &vars_buf
) const
 /*
	* Function that computes the requested conservative and derived visualization variables.
	*/
{
	static_assert((std::is_same_v<T,float> && sizeof(T) == 4) ||
	              (std::is_same_v<T,double> && sizeof(T) == 8));
	static_assert(std::numeric_limits<T>::is_iec559);

	if( mPointFields.empty() ) return;

	// Get the total number of elements in all zones.
	const std::size_t nElemTotal = multizone_geometry_container->GetnElemTotal();

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
	{
		// Each thread reuses its workspace as long as the element shape matches.
		CMatrixAS3<as3double> primvar;
		CMatrixAS3<as3double> dUDx;
		CMatrixAS3<as3double> dUDy;
		CMatrixAS3<as3double> identity;

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
		for(std::size_t i=0; i<nElemTotal; i++)
		{
			// Extract the element indices.
			const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

			// Deduce the current element's zone and index.
			const auto iZone  = elem_info.mIndexZone;
			const auto ijElem = elem_info.mIndexElem;

			// Extract current grid zone.
			auto* zone = multizone_geometry_container->GetSinglezoneGeometry(iZone);
			// Extract current solver.
			const auto* solver = multizone_solver_container->GetSinglezoneSolver(iZone);

			// Extract the total number of elements in this zone.
			const std::size_t nElem = zone->GetnElem();
			// Extract the current zone information.
			const std::size_t nNode = zone->GetnNodeGrid2D();
			// Total DOFs and sub-elements in this zone.
			const std::size_t nDOFs = nNode*nElem;

			// Extract solution in this element.
			auto& sol = solver->GetPhysicalElement(ijElem)->mSol2D;

			// Allocate data for the primitive variables on the solution DOFs for this element.
			// These contain: 1/rho, u, v, p.
			// Note, the inverse of rho is computed for convenience.
			if( primvar.col() != nNode ) primvar.resize(4, nNode);

			// Loop over the element DOFs and compute the primitive variables needed.
			for(std::size_t l=0; l<nNode; l++)
			{
				const as3double ovrho =  C_ONE/sol(0,l);
				const as3double u     =  ovrho*sol(1,l);
				const as3double v     =  ovrho*sol(2,l);
				const as3double p     =  C_GM1*( sol(3,l) - C_HALF*sol(0,l)*(u*u + v*v) );

				// Store the values.
				primvar(0,l) = ovrho;
				primvar(1,l) = u;
				primvar(2,l) = v;
				primvar(3,l) = p;
			}

			// Evaluate each configured field at its offset in the zone's buffer.
			for(const auto& field : mPointFields)
			{
				const std::size_t iVar = field.mComponentOffset;
				switch( field.mVariable )
				{

					case(EWriteVariable::DENSITY):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the density in each DOF.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( sol(0,l) );
						}

						break;
					}

					case(EWriteVariable::MOMENTUM):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode*3;

						// Compute and store the momentum variables.
						for(std::size_t l=0; l<nNode; l++)
						{
							const std::size_t l3 = 3*l;
							buf[l3  ] = static_cast<T>( sol(1,l) );
							buf[l3+1] = static_cast<T>( sol(2,l) );
							buf[l3+2] = static_cast<T>(0);
						}

						break;
					}

					case(EWriteVariable::TOTAL_ENERGY):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the energy variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( sol(3,l) );
						}

						break;
					}

					case(EWriteVariable::PRESSURE):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the pressure variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( primvar(3,l) );
						}

						break;
					}

					case(EWriteVariable::VELOCITY):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode*3;

						// Compute and store the velocity variables.
						for(std::size_t l=0; l<nNode; l++)
						{
							const std::size_t l3 = 3*l;
							const as3double ovrho = primvar(0,l);
							buf[l3  ] = static_cast<T>( ovrho*sol(1,l) );
							buf[l3+1] = static_cast<T>( ovrho*sol(2,l) );
							buf[l3+2] = static_cast<T>(0);
						}

						break;
					}

					case(EWriteVariable::TEMPERATURE):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Abbreviation for 1/R.
						const as3double ovR = C_ONE/C_RGAS;

						// Compute and store the temperature variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( ovR*primvar(0,l)*primvar(3,l) );
						}

						break;
					}

					case(EWriteVariable::MACH):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the Mach number variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							const as3double a2    = C_GMA*primvar(0,l)*primvar(3,l);
							const as3double umag2 = primvar(1,l)*primvar(1,l)
								                    + primvar(2,l)*primvar(2,l);
							buf[l] = static_cast<T>( std::sqrt( umag2/a2 ) );
						}

						break;
					}

					case(EWriteVariable::ENTROPY):
					{
						// Set the pointer to the correct location.
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the entropy variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( std::log( primvar(3,l)*std::pow( primvar(0,l), C_GMA ) ) );
						}

						break;
					}

					case(EWriteVariable::VORTICITY):
					{
						// Allocate memory for the velocity gradient.
						if( dUDx.col() != nNode )
						{
							dUDx.resize(2, nNode);
							dUDy.resize(2, nNode);
						}

						// Extract number of solution points in 1D.
						const std::size_t nSol1D = solver->GetStandardElement()->GetnSol1D();

						// Create an identity matrix, needed for the Lagrange interpolation over the DOFs.
						if( identity.row() != nSol1D )
							identity = NLinearAlgebra::CreateIdentityMatrix<as3double>(nSol1D);

						// Extract the transposed differentiation matrix in 1D over the solution points.
						const auto& derLagrangeSol1DTrans = solver->GetStandardElement()->GetDerLagrangeSol1DTrans();

						// Compute the gradient of the velocity.
						solver->GetTensorProduct()->RuntimeVolume(nSol1D, 2, nSol1D,
						                                          identity.data(), derLagrangeSol1DTrans.data(),
						                                          &primvar[nNode], nullptr,
						                                          dUDx.data(), dUDy.data());

						// Extract the metrics at the volume solution points.
						auto& metrics = solver->GetPhysicalElement(ijElem)->mMetricSol2D;

						// Convert the derivatives from parametric space to Cartesian space.
						for(std::size_t l=0; l<nNode; l++)
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
						T *buf = vars_buf[iZone].data() + iVar*nDOFs + ijElem*nNode;

						// Compute and store the Vorticity variable.
						for(std::size_t l=0; l<nNode; l++)
						{
							buf[l] = static_cast<T>( dUDx(1,l) - dUDy(0,l) );
						}

						break;
					}


					default: ERROR("Unknown variable detected.");
				}
			}


		}
	}
}

// Instantiate the two file precisions used by the binary writers.
template void IFileVTK::DetermineVisualizationData(const CMultizoneGeometry *,
                                                    const CMultizoneSolver   *,
                                                    as3vector2d<float>       &) const;

template void IFileVTK::DetermineVisualizationData(const CMultizoneGeometry *,
                                                    const CMultizoneSolver   *,
                                                    as3vector2d<double>      &) const;
