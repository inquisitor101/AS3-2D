#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "temporal_structure.hpp"
#include "vtk_structure.hpp"
#include "standard_element_structure.hpp"
#include "physical_element_structure.hpp"
#include "solver_structure.hpp"
#include "spatial_structure.hpp"
#include "tensor_structure.hpp"
#include "initial_condition_structure.hpp"


// Forward declaration to avoid compiler issues.
class ISolver;
class IFileVTK;
class IInitialCondition;

/*!
 * @brief A class of factories that creates different containers. 
 */
class CGenericFactory
{
	public:

		/*!
		 * @brief Function that creates a specialized instance of a temporal container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 *
		 * @return unique pointer to the specific temporal class.
		 */
		static std::unique_ptr<ITemporal> CreateTemporalContainer(CConfig *config_container);

		/*!
		 * @brief Function that creates a specialized instance of a spatial container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 *
		 * @return unique pointer to the specific spatial class.
		 */
		static std::unique_ptr<ISpatial> CreateSpatialContainer(CConfig        *config_container,
				                                                    CGeometry      *geometry_container,
																														unsigned short  iZone);



		/*!
		 * @brief Function that creates a specialized instance of a vtk container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 *
		 * @return unique pointer to the specific vtk class.
		 */
		static std::unique_ptr<IFileVTK> CreateVTKContainer(CConfig   *config_container,
				                                                CGeometry *geometry_container);


		/*!
		 * @brief Function that creates a specialized instance of an initial condition container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 *
		 * @return unique pointer to the specific initial condition class.
		 */
		static std::unique_ptr<IInitialCondition> CreateInitialConditionContainer(CConfig *config_container);



		/*!
		 * @brief Function that creates a specialized instance of a standard element container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] iZone zone ID of this container.
		 *
		 * @return unique pointer to specialized standard element class.
		 */
		static std::unique_ptr<CStandardElement> CreateStandardElement(CConfig       *config_container,
				                                                           unsigned short iZone);


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
		static std::unique_ptr<CPhysicalElement> CreatePhysicalElement(CConfig          *config_container,
				                                                           CStandardElement *standard_element,
																																	 ITensorProduct   *tensor_container,
				                                                           CElementGeometry *element_geometry,
				                                                           unsigned short    iZone,
																																	 unsigned short    nVar);

		/*!
		 * @brief Function that creates a specialized instance of a templated tensor container.
		 *
		 * @param[in] standard_element pointer to the relevant standar element container.
		 *
		 * @return unique pointer to specialized template tensor-product class.
		 */
		static std::unique_ptr<ITensorProduct> CreateTensorContainer(CStandardElement *standard_element);

		/*!
		 * @brief Function that creates a specialized instance of a solver container.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iZone zone ID of this container.
		 *
		 * @return unique pointer to the specific solver class.
		 */
		static std::unique_ptr<ISolver> CreateSolverContainer(CConfig       *config_container,
				                                                  CGeometry     *geometry_container,
																													unsigned short iZone);

		/*!
		 * @brief Function that creates a vector of specialized instances of solver containers.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 *
		 * @return vector of specilized solver containers in each zone.
		 */
		static as3vector1d<std::unique_ptr<ISolver>> CreateMultizoneSolverContainer(CConfig   *config_container,
				                                                                        CGeometry *geometry_container);

		

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

