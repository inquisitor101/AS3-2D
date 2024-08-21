#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "tensor_structure.hpp"
#include "factory_structure.hpp"
#include "boundary_structure.hpp"
#include "standard_element_structure.hpp"
#include "physical_element_structure.hpp"


/*!
 * @brief An interface class used for the solver specification.
 */
class ISolver
{
	public:
		
		/*!
		 * @brief Constructor of ISolver, which serves as an interface for the solver.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iZone zone ID of this container.
		 */
		ISolver(CConfig       *config_container,
				    CGeometry     *geometry_container,
						unsigned short iZone);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ISolver(void);


		/*!
		 * @brief Pure virtual function that initializes the physical elements. Must be implemented in a derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		virtual void InitPhysicalElements(CConfig   *config_container,
				                              CGeometry *geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that initializes the boundary conditions. Must be implemented in a derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		virtual void InitBoundaryConditions(CConfig   *config_container,
				                                CGeometry *geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that computes the volume terms over a single element.
		 * 
		 * @param[in] iElem element index.
		 * @param[in] localtime current physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		virtual void ComputeVolumeResidual(size_t                     iElem,
				                               as3double                  localtime,
																			 as3vector1d<as3double>    &monitordata,
																			 CPoolMatrixAS3<as3double> &workarray) = 0;

		/*!
		 * @brief Pure virtual getter function which returns the number of working variables. Must be overridden.
		 */
		virtual unsigned short GetnVar(void) const = 0;

		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		unsigned short GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the tensor-product container of this zone.
		 *
		 * @return mTensorProductContainer.
		 */
		ITensorProduct *GetTensorProduct(void) const {return mTensorProductContainer.get();}

		/*!
		 * @brief Getter function which returns the standard element container of this zone.
		 *
		 * @return mStandardElementContainer.
		 */
		const CStandardElement *GetStandardElement(void) const {return mStandardElementContainer.get();}

		/*!
		 * @brief Getter function which returns a specific physical element container.
		 *
		 * @param[in] index index of the physical element.
		 *
		 * @return mPhysicalElementContainer[index].
		 */
		CPhysicalElement *GetPhysicalElement(size_t index) const {return mPhysicalElementContainer[index].get();}

		/*!
		 * @brief Getter function which returns the entire boundary container.
		 *
		 * @return mBoundaryContainer.
		 */
		as3vector1d<std::unique_ptr<IBoundary>> &GetBoundaryContainer(void) {return mBoundaryContainer;}

	protected:
		const unsigned short                           mZoneID;                    ///< Zone ID of this container.
	
		std::unique_ptr<ITensorProduct>                mTensorProductContainer;    ///< Tensor product container.
		std::unique_ptr<CStandardElement>              mStandardElementContainer;  ///< Standard element container.	
		as3vector1d<std::unique_ptr<CPhysicalElement>> mPhysicalElementContainer;  ///< Physical element container.
		as3vector1d<std::unique_ptr<IBoundary>>        mBoundaryContainer;         ///< Boundary condition container.

	private:
		// Disable default constructor.
		ISolver(void) = delete;
		// Disable default copy constructor.
		ISolver(const ISolver&) = delete;
		// Disable default copy operator.
		ISolver& operator=(ISolver&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a solver specification based on the (non-linear) Euler equations. 
 */
class CEESolver : public ISolver
{
	public:

		/*!
		 * @brief Constructor of CEESolver, which initializes an Euler equations solver.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iZone zone ID of this container.
		 */
		CEESolver(CConfig       *config_container,
				      CGeometry     *geometry_container,
							unsigned short iZone);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEESolver(void) override;

		/*!
		 * @brief Function that initializes the physical elements.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitPhysicalElements(CConfig   *config_container,
				                      CGeometry *geometry_container) override;

		/*!
		 * @brief Function that initializes the boundary conditions. 
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitBoundaryConditions(CConfig   *config_container,
		                            CGeometry *geometry_container) override;

		/*!
		 * @brief Function that computes the volume terms over a single element, based on the EE equations.
		 * 
		 * @param[in] iElem element index.
		 * @param[in] localtime current physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void ComputeVolumeResidual(size_t                     iElem,
		                           as3double                  localtime,
															 as3vector1d<as3double>    &monitordata,
															 CPoolMatrixAS3<as3double> &workarray) override;

		/*!
		 * @brief Getter function which returns the number of working variables.
		 *
		 * @return mNVar.
		 */
		unsigned short GetnVar(void) const override {return mNVar;}

	protected:

	private:
		unsigned short mNVar = 4; ///< Number of working variables
};


