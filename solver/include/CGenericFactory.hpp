#pragma once

/*!
 * @file
 * @brief Construction of solver components and supported tensor specializations.
 */


#include <cstddef>
#include <memory>

#include "visualization/IFileVTK.hpp"
#include "compute_kernels/TensorProducts.hpp"
#include "common/AS3Options.hpp"
#include "CConfig.hpp"
#include "temporal/ITemporal.hpp"
#include "boundary/IBoundary.hpp"
#include "interface/IInterface.hpp"
#include "CMonitorData.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "element/CStandardElement.hpp"
#include "element/CPhysicalElement.hpp"
#include "ic/IInitialCondition.hpp"

// Forward declaration to avoid compiler issues.
class CMultizoneGeometry;
class CSinglezoneGeometry;
class CElementGeometry;
class CBoundaryFaceBatch;
class CInterfaceFaceBatch;
class CMultizoneSolver;
class ISinglezoneSolver;
class IFileVTK;
class ITemporal;
class IInterface;
class IBoundary;
class IInitialCondition;

/*!
 * @brief A class of factories that creates different containers.
 */
class CGenericFactory
{
	public:

		/*!
		 * @brief Function that constructs the simulation monitoring data container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 *
		 * @return unique pointer to the specific data monitoring class.
		 */
		static std::unique_ptr<CMonitorData> CreateMonitoringContainer(const CConfig *config_container);

		/*!
		 * @brief Function that creates a specialized instance of a temporal container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 *
		 * @return unique pointer to the specific temporal class.
		 */
		static std::unique_ptr<ITemporal> CreateTemporalContainer(const CConfig *config_container);

		/*!
		 * @brief Function that creates a specialized instance of a vtk container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 *
		 * @return unique pointer to the specific vtk class.
		 */
		static std::unique_ptr<IFileVTK> CreateVTKContainer(const CConfig            *config_container,
				                                                const CMultizoneGeometry *multizone_geometry_container);

		/*!
		 * @brief Function that creates a specialized instance of an initial condition container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 *
		 * @return unique pointer to the specific initial condition class.
		 */
		static std::unique_ptr<IInitialCondition> CreateInitialConditionContainer(const CConfig *config_container);

		/*!
		 * @brief Function that creates a specialized instance of a Riemann solver container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] riemann type of Riemann solver.
		 *
		 * @return unique pointer to the specific Riemann solver class.
		 */
		static std::unique_ptr<IRiemannSolver> CreateRiemannSolverContainer(const CConfig      *config_container,
				                                                                ETypeRiemannSolver  riemann);

		/*!
		 * @brief Function that constructs a runtime reference element from the zone geometry order and
		 *   integration rule.
		 *
     * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 *
		 * @return unique pointer to specialized standard element class.
		 */
		static std::unique_ptr<CStandardElement> CreateStandardElement(const CConfig             *config_container,
                                                                   const CSinglezoneGeometry *singlezone_geometry_container);

		/*!
		 * @brief Function that creates a specialized instance of a physical element container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] standard_element standard element container.
		 * @param[in] tensor_container tensor product container.
		 * @param[in] element_geometry geometry of the physical element.
		 * @param[in] iZone zone ID of this container.
		 * @param[in] nVar number of solution variables to store.
		 *
		 * @return unique pointer to specialized physical element class.
		 */
		static std::unique_ptr<CPhysicalElement> CreatePhysicalElement(const CConfig    *config_container,
				                                                           CStandardElement *standard_element,
																																	 ITensorProduct   *tensor_container,
				                                                           CElementGeometry *element_geometry,
				                                                           std::size_t    iZone,
																																	 std::size_t    nVar);

		/*!
		 * @brief Function that creates a specialized instance of a templated tensor container.
		 *
		 * @param[in] standard_element pointer to the relevant standard element container.
		 * @param[in] nVar number of solution variables.
		 * @param[in] typeDGMethod DG formulation; DG-SEM requires collocated nodes.
		 *
		 * @return unique pointer to specialized template tensor-product class.
		 */
		static std::unique_ptr<ITensorProduct> CreateTensorContainer(CStandardElement *standard_element,
                                                                 std::size_t    nVar,
                                                                 ETypeDGMethod typeDGMethod = ETypeDGMethod::NODAL_DG);


		/*!
		 * @brief Function that creates a specialized instance of a boundary container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] zone_geometry geometry of the current zone.
		 * @param[in] zone_solver solver for the current zone.
		 * @param[in] boundary_batch physical boundary face batch.
		 *
		 * @return unique pointer to the specific boundary class.
		 */
		static std::unique_ptr<IBoundary> CreateBoundaryContainer(const CConfig              *config_container,
				                                                      const CSinglezoneGeometry  *zone_geometry,
                                                              const ISinglezoneSolver              *zone_solver,
																															const CBoundaryFaceBatch   &boundary_batch);


		/*!
		 * @brief Function that constructs the interface evaluator for a paired face batch.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_solver_container solvers and residual storage for all zones.
		 * @param[in] interface_batch paired-zone interface face batch.
		 *
		 * @return newly constructed object or matrix with ownership transferred to the caller.
		 */
		static std::unique_ptr<IInterface> CreateInterfaceContainer(const CConfig               *config_container,
																													      CMultizoneSolver            *multizone_solver_container,
                                                                const CInterfaceFaceBatch   &interface_batch);


		/*!
		 * @brief Function that constructs a zone solver for the configured equations.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry geometry of the current zone.
		 *
		 * @return newly constructed object or matrix with ownership transferred to the caller.
		 */
		static std::unique_ptr<ISinglezoneSolver> CreateSolverContainer(const CConfig             *config_container,
				                                                  const CSinglezoneGeometry *singlezone_geometry);

	private:
		// Disable default constructor.
		CGenericFactory(void)  = delete;		
		// Disable default destructor.
		~CGenericFactory(void) = delete;
		// Disable default copy constructor.
		CGenericFactory(const CGenericFactory&) = delete;
		// Disable default copy operator.
		CGenericFactory& operator=(CGenericFactory&) = delete;	
};

