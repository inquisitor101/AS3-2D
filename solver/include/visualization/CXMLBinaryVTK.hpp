#pragma once

/*!
 * @file
 * @brief VTK XML unstructured grid writer with raw appended binary arrays.
 */


#include <cstddef>
#include <cstdint>
#include <string>
#include <variant>

#include "visualization/IFileVTK.hpp"


/*!
 * @brief Writer for binary .vtu snapshots of the solution on a fixed grid.
 */
class CXMLBinaryVTK final : public IFileVTK
{
	public:

		/*!
		 * @brief Constructor that initializes the VTK XML writer.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		CXMLBinaryVTK(const CConfig            *config_container,
		              const CMultizoneGeometry *multizone_geometry_container);

		/*!
		 * @brief Function that writes the current solution and the cached grid to a .vtu file.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] openmp_container OpenMP team configuration.
		 * @param[in] multizone_solver_container solvers and solution storage for all zones.
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		void WriteFileVTK(const CConfig            *config_container,
		                  const CMultizoneGeometry *multizone_geometry_container,
		                  const COpenMP            *openmp_container,
		                  const CMultizoneSolver   *multizone_solver_container,
		                  as3double                 time) override;

	private:
		/*!
		 * @brief Cached coordinates and reusable point fields for one output precision.
		 *
		 * @tparam T floating-point type written to the file.
		 */
		template<typename T>
		struct CRealData
		{
			as3vector1d<T> mPoints;     ///< Cached XYZ coordinates.
			as3vector2d<T> mPointData;  ///< Reused zone-local buffers for solution fields.
		};

		std::variant<CRealData<float>, CRealData<double>> mRealData; ///< Buffers in the configured precision.
		as3vector1d<std::size_t>  mZonePointCounts;  ///< Number of output points in each zone.
		as3vector1d<std::int32_t> mConnectivity;     ///< Cached point indices for each quadrilateral.
		as3vector1d<std::int32_t> mOffsets;          ///< Cached end index of each cell in mConnectivity.
		as3vector1d<std::uint8_t> mCellTypes;        ///< Cached VTK cell type for each cell.
		std::string               mFileHeader;       ///< XML metadata and the raw appended-data marker.

		/*!
		 * @brief Function that builds the fixed mesh and allocates the reusable output buffers.
		 *
		 * @tparam T floating-point type written to the file.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[out] realData coordinates and point-field buffers in the selected precision.
		 */
		template<typename T>
		void InitializeMesh(const CMultizoneGeometry *multizone_geometry_container,
		                    CRealData<T>             &realData);

		/*!
		 * @brief Function that builds the XML header with byte offsets for all appended arrays.
		 *
		 * @tparam T floating-point type written to the file.
		 * @param[in] realData coordinates and point-field buffers in the selected precision.
		 */
		template<typename T>
		void InitializeFileHeader(const CRealData<T> &realData);

		/*!
		 * @brief Function that writes a VTK XML snapshot in the selected precision.
		 *
		 * @tparam T floating-point type written to the file.
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] multizone_solver_container solvers and solution storage for all zones.
		 * @param[in] time physical simulation time of the snapshot [s].
		 * @param[in,out] realData cached coordinates and reused point-field buffers.
		 */
		template<typename T>
		void WriteSnapshot(const CConfig            *config_container,
		                   const CMultizoneGeometry *multizone_geometry_container,
		                   const CMultizoneSolver   *multizone_solver_container,
		                   as3double                 time,
		                   CRealData<T>             &realData);
};
