#pragma once

/*!
 * @file
 * @brief VTK visualization writer interface.
 */


#include <cstddef>
#include <string>

#include "common/AS3Options.hpp"

class CConfig;
class CMultizoneGeometry;
class CMultizoneSolver;
class COpenMP;


/*!
 * @brief An interface class used for writing VTK files.
 */
class IFileVTK
{
	public:
	

		/*!
		 * @brief Constructor of IFileVTK, which serves as an interface for VTK-type formats.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		IFileVTK(const CConfig            *config_container,
				     const CMultizoneGeometry *multizone_geometry_container);
		

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IFileVTK(void);

		/*!
		 * @brief Pure virtual function that writes VTK data to a file. Must be implemented by a derived
		 *   class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] openmp_container OpenMP team configuration.
		 * @param[in] multizone_solver_container solvers and residual storage for all zones.
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		virtual void WriteFileVTK(const CConfig            *config_container,
		                          const CMultizoneGeometry *multizone_geometry_container,
		                          const COpenMP            *openmp_container,
		                          const CMultizoneSolver   *multizone_solver_container,
		                          as3double                 time) = 0;

	protected:
		/*!
		 * @brief Layout of one scalar or vector field in each zone's point-data buffer.
		 */
		struct CPointField
		{
			EWriteVariable mVariable;         ///< Solution variable to evaluate.
			const char    *mName;             ///< Name of the VTK array.
			std::size_t    mComponentOffset;  ///< Number of components in preceding fields.
			std::size_t    mNComponent;       ///< Number of interleaved components at each point.
		};

		as3vector1d<CPointField> mPointFields;            ///< Point fields in configured output order.
		std::size_t             mNPointComponents = 0;  ///< Total components per point across all fields.

		/*!
		 * @brief Function that returns the index of the next snapshot.
		 */
		std::size_t GetFileNumber(void) const { return mSnapshots.size(); }

		/*!
		 * @brief Function that adds a completed snapshot to the ParaView time-series index.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] filename path of the completed snapshot.
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		void WriteSeriesIndex(const CConfig      *config_container,
		                      const std::string &filename,
		                      as3double          time);

		/*!
		 * @brief Function that computes the requested visualization variables for both VTK formats.
		 *
		 * @tparam T floating-point type written to the file.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] multizone_solver_container solvers and solution storage for all zones.
		 * @param[out] vars_buf preallocated field buffers indexed by zone.
		 */
		template<typename T>
		void DetermineVisualizationData(const CMultizoneGeometry *multizone_geometry_container,
		                                const CMultizoneSolver   *multizone_solver_container,
		                                as3vector2d<T>           &vars_buf) const;

	private:
		/*!
		 * @brief One completed snapshot in the ParaView time-series index.
		 */
		struct CSnapshot
		{
			std::string mName;  ///< Snapshot filename relative to the index directory.
			as3double   mTime;  ///< Physical simulation time [s].
		};

		as3vector1d<CSnapshot> mSnapshots;  ///< Completed snapshots in output order.

		// Disable default constructor.
		IFileVTK(void) = delete;
		// Disable default copy constructor.
		IFileVTK(const IFileVTK&) = delete;
		// Disable default copy operator.
		IFileVTK& operator=(IFileVTK&) = delete;	
};
