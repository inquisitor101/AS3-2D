#pragma once

/*!
 * @file
 * @brief Legacy binary VTK visualization of the multizone solution.
 */


#include <cstddef>
#include <string>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "openmp_structure.hpp"


// Forward declaration to avoid compiler issues.
class CMultizoneSolver;
class ISolver;
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
		 * @brief Pure virtual function that writes VTK data a file. Must be implemented by a derived
		 *   class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] openmp_container OpenMP team configuration.
		 * @param[in] multizone_solver_container solvers and residual storage for all zones.
		 */
		virtual void WriteFileVTK(const CConfig            *config_container,
															const CMultizoneGeometry *multizone_geometry_container,
															const COpenMP            *openmp_container,
															const CMultizoneSolver   *multizone_solver_container) = 0;

	protected:
		as3vector1d<std::string> mVariableNames;  ///< Variable names for writing.
		std::size_t            mFileNumber;       ///< Written file index.

		bool                     mWriteDensity;      ///< Option for writing a density variable.
		bool                     mWriteMomentum;     ///< Option for writing a momentum variable.
		bool                     mWriteTotalEnergy;  ///< Option for writing a total energy variable.
		bool                     mWritePressure;     ///< Option for writing a pressure variable.
		bool                     mWriteVelocity;     ///< Option for writing a velocity variable.
		bool                     mWriteVorticity;    ///< Option for writing a vorticity variable.
		bool                     mWriteMach;         ///< Option for writing a Mach number variable.
		bool                     mWriteTemperature;  ///< Option for writing a temperature variable.
		bool                     mWriteEntropy;      ///< Option for writing a specific entropy variable.

	private:
		// Disable default constructor.
		IFileVTK(void) = delete;
		// Disable default copy constructor.
		IFileVTK(const IFileVTK&) = delete;
		// Disable default copy operator.
		IFileVTK& operator=(IFileVTK&) = delete;	
};

//-----------------------------------------------------------------------------------

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
		 */
		void WriteFileVTK(const CConfig            *config_container,
											const CMultizoneGeometry *multizone_geometry_container,
											const COpenMP            *openmp_container,
											const CMultizoneSolver   *multizone_solver_container) override;

	protected:

	private:
		bool mBigEndian; ///< Whether this machine uses big endian format.
	
		// Disable default constructor.
		CLegacyBinaryVTK(void) = delete;
	

	  /*!
	   * @brief Function that computes and stores the data required for writing a binary VTK file.
	   *
	   * @param[in] config_container pointer to the configuration container.
	   * @param[in] multizone_geometry_container pointer to geometry container.
	   * @param[in] openmp_container OpenMP shared parallelization container.
	   * @param[in] multizone_solver_container pointer to the multizone solver container.
	   * @param[out] vars_buf reference to the variables buffer data.
	   */
	  void DetermineVisualizationData(const CConfig            *config_container,
				                            const CMultizoneGeometry *multizone_geometry_container,
																		const COpenMP            *openmp_container,
				                            const CMultizoneSolver   *multizone_solver_container,
																		as3vector2d<float>       &vars_buf);
};

