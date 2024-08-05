#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "tensor_structure.hpp"
#include "spatial_structure.hpp"
#include "factory_structure.hpp"
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


		/*
		 * @brief Pure virtual function that initializes the physical elements. Must be implemented in a derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		virtual void InitPhysicalElements(CConfig   *config_container,
				                              CGeometry *geometry_container) = 0;
		
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



	protected:
		const unsigned short                           mZoneID;                    ///< Zone ID of this container.
		std::unique_ptr<ISpatial>                      mSpatialContainer;          ///< Spatial container.
		std::unique_ptr<ITensorProduct>                mTensorProductContainer;    ///< Tensor product container.
		std::unique_ptr<CStandardElement>              mStandardElementContainer;  ///< Standard element container.	
		as3vector1d<std::unique_ptr<CPhysicalElement>> mPhysicalElementContainer;  ///< Physical element container.

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

		/*
		 * @brief Function that initializes the physical elements.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitPhysicalElements(CConfig   *config_container,
				                      CGeometry *geometry_container) override;


	protected:

	private:
		unsigned short mNVar = 4; ///< Number of working variables
};


