#pragma once

/*!
 * @file
 * @brief Legacy binary VTK writer.
 */


#include <string>

#include "common/AS3Options.hpp"
#include "visualization/IFileVTK.hpp"


/*!
 * @brief A class used for writing VTK files in binary and legacy format.
 */
class CLegacyBinaryVTK final : public IFileVTK
{
	public:
	

		/*!
		 * @brief Constructor of CLegacyBinaryVTK, which is responsible for VTK in binary and legacy
		 *   format.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		CLegacyBinaryVTK(const CConfig            *config_container,
				             const CMultizoneGeometry *multizone_geometry_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CLegacyBinaryVTK(void) override;

		/*!
		 * @brief Function that writes VTK data to a file using a binary and legacy format.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] multizone_geometry_container pointer to the geometry container.
		 * @param[in] openmp_container OpenMP shared parallelization container.
		 * @param[in] multizone_solver_container reference to the solver container.
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		void WriteFileVTK(const CConfig            *config_container,
		                  const CMultizoneGeometry *multizone_geometry_container,
		                  const COpenMP            *openmp_container,
		                  const CMultizoneSolver   *multizone_solver_container,
		                  as3double                 time) override;

	protected:

	private:
		bool mBigEndian; ///< Whether this machine uses big endian format.

		/*!
		 * @brief Function that writes a legacy binary snapshot in the selected precision.
		 *
		 * @tparam T floating-point type written to the file.
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] multizone_solver_container solvers and solution storage for all zones.
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		template<typename T>
		void WriteSnapshot(const CConfig            *config_container,
		                   const CMultizoneGeometry *multizone_geometry_container,
		                   const CMultizoneSolver   *multizone_solver_container,
		                   as3double                 time);
	
		// Disable default constructor.
		CLegacyBinaryVTK(void) = delete;
};
