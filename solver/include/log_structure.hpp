#pragma once

/*!
 * @file
 * @brief Console output for configuration, grids, monitoring and parallel execution.
 */


#include <cstddef>
#include <ostream>
#include <string>

#include <iosfwd>
#include <iomanip>
#include "option_structure.hpp"
#include "config_structure.hpp"
#include "solver_structure.hpp"
#include "monitoring_structure.hpp"
#include "import_structure.hpp"
#include "marker_structure.hpp"


/*!
 * @brief A namespace used for storing logging functionalities.
 */
namespace NLogger
{

	/*!
	 * @brief Function that prints the entries of a matrix.
	 *
	 * @tparam T matrix scalar type.
	 *
	 * @param[in] m reference to the matrix.
	 * @param[in] g log message printed.
	 */
	template<typename T>
	void PrintMatrixAS3(CMatrixAS3<T> &m, std::string g)
	{
		using namespace std;
		cout << "\n" << g << "\n";
		for(std::size_t i=0; i<m.row(); i++)
		{
			for(std::size_t j=0; j<m.col(); j++)
			{
				cout << fixed << showpos << scientific << m(i,j) << " ";
			}
			cout << endl;
		}
	}

	/*!
	 * @brief Function that prints the monitored data and the header info.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 * @param[in] monitor_container data monitoring container.
	 * @param[in] iter data monitored.
	 * @param[in] time current physical time.
	 * @param[in] step current time step.
	 */
	void MonitorOutput(CConfig       *config_container,
			               CMonitorData  *monitor_container,
			               std::size_t  iter,
			               as3double      time,
										 as3double      step);

	/*!
	 * @brief Function that prints the type of solver in each zone.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 */
	void PrintInitSolver(CConfig *config_container);

	/*!
	 * @brief Function that displays the boundary condition information over all zones.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 * @param[in] multizone_geometry_container input geometry container.
	 * @param[in] multizone_solver_container solvers and residual storage for all zones.
	 */
	void DisplayBoundaryConditions(const CConfig            *config_container,
																 const CMultizoneGeometry *multizone_geometry_container,
																 const CMultizoneSolver   *multizone_solver_container);


	/*!
	 * @brief Function that displays the OpenMP information, if any.
	 *
	 * @param[in] openmp_container OpenMP shared parallelization container.
	 * @param[in] multizone_geometry_container geometry of all grid zones.
	 * @param[in] multizone_solver_container input multizone solver container.
	 */
	void DisplayOpenMPInfo(COpenMP                  *openmp_container,
			                   const CMultizoneGeometry *multizone_geometry_container,
												 const CMultizoneSolver   *multizone_solver_container);


  /*!
   * @brief Function that prints AS3 file metadata and identifier mappings.
   *
   * @param[in] metadata previously read AS3 header and zone metadata.
   * @param[in] out output stream.
   */
  void DisplayAS3BinaryMetadata(const NImportFile::NAS3BinaryFile::CAS3BinaryMetadata &metadata,
                                std::ostream                                          &out);

  /*!
   * @brief Function that prints internal zone-edge pairings.
   *
   * @param[in] families zone-edge marker families to display.
   */
  void DisplayInternalMarkers(const as3vector1d<CInternalFamilyMarker> &families);

  /*!
   * @brief Function that prints named physical boundary markers.
   *
   * @param[in] families zone-edge marker families to display.
   */
  void DisplayExternalMarkers(const as3vector1d<CExternalFamilyMarker> &families);

  /*!
   * @brief Function that prints periodic zone-edge pairings and translations.
   *
   * @param[in] families zone-edge marker families to display.
   */
  void DisplayPeriodicMarkers(const as3vector1d<CPeriodicFamilyMarker> &families);
}




