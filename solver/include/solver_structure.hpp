#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "tensor_structure.hpp"
#include "factory_structure.hpp"
#include "boundary_structure.hpp"
#include "riemann_solver_structure.hpp"
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
		 * @brief Pure virtual function that returns the type of solver. Must be overridden.
		 */
		virtual ETypeSolver GetTypeSolver(void) const = 0;

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
		 * @brief Pure virtual function that computes the volume terms in the entire solver.
		 * 
		 * @param[in] localtime current physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 */
		virtual void ComputeVolumeResidual(as3double                  localtime,
																			 as3vector1d<as3double>    &monitordata,
																			 CPoolMatrixAS3<as3double> &workarray) = 0;

		/*!
		 * @brief Pure virtual function that computes the surface terms in the i-direction in the entire solver.
		 *
		 * @param[in] geometry_container input geometry container.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 */
		virtual void ComputeSurfaceResidualIDir(CGeometry                 *geometry_container,
																			      as3vector1d<as3double>    &monitordata,
																			      CPoolMatrixAS3<as3double> &workarray,
																						as3double                  localtime) = 0;

		/*!
		 * @brief Pure virtual function that computes the surface terms in the j-direction in the entire solver.
		 *
		 * @param[in] geometry_container input geometry container.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 */
		virtual void ComputeSurfaceResidualJDir(CGeometry                 *geometry_container,
																			      as3vector1d<as3double>    &monitordata,
																			      CPoolMatrixAS3<as3double> &workarray,
																						as3double                  localtime) = 0;

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
		 * @brief Getter function which returns the Riemann solver container of this zone.
		 *
		 * @return mRiemannSolverContainer.
		 */
		IRiemannSolver *GetRiemannSolver(void) const {return mRiemannSolverContainer.get();}

		/*!
		 * @brief Getter function which returns the standard element container of this zone.
		 *
		 * @return mStandardElementContainer.
		 */
		const CStandardElement *GetStandardElement(void) const {return mStandardElementContainer.get();}

		/*!
		 * @brief Getter function which returns the entire physical element container.
		 *
		 * @return mPhysicalElementContainer.
		 */
		as3vector1d<std::unique_ptr<CPhysicalElement>> &GetPhysicalElement(void) {return mPhysicalElementContainer;}

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
		std::unique_ptr<ITensorProduct>                mTensorProductContainer;    ///< Tensor product container.
		std::unique_ptr<IRiemannSolver>                mRiemannSolverContainer;    ///< Riemann solver container.
		std::unique_ptr<CStandardElement>              mStandardElementContainer;  ///< Standard element container.	
		as3vector1d<std::unique_ptr<CPhysicalElement>> mPhysicalElementContainer;  ///< Physical element container.
		
		as3vector1d<std::unique_ptr<IBoundary>>        mBoundaryIMINContainer;     ///< IMIN boundary container.
		as3vector1d<std::unique_ptr<IBoundary>>        mBoundaryIMAXContainer;     ///< IMAX boundary container.
		as3vector1d<std::unique_ptr<IBoundary>>        mBoundaryJMINContainer;     ///< JMIN boundary container.
		as3vector1d<std::unique_ptr<IBoundary>>        mBoundaryJMAXContainer;     ///< JMAX boundary container.

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
		 * @brief Function that returns the type of solver. 
		 */
		ETypeSolver GetTypeSolver(void) const override {return ETypeSolver::EE;}

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
		 * @brief Function that computes the volume terms in the entire solver, based on the EE.
		 * 
		 * @param[in] localtime current physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 */
		void ComputeVolumeResidual(as3double                  localtime,
															 as3vector1d<as3double>    &monitordata,
															 CPoolMatrixAS3<as3double> &workarray) override;

		/*!
		 * @brief Function that computes the surface terms in the i-direction in the entire solver, based on the EE.
		 *
		 * @param[in] geometry_container input geometry container.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 */
		void ComputeSurfaceResidualIDir(CGeometry                 *geometry_container,
															      as3vector1d<as3double>    &monitordata,
															      CPoolMatrixAS3<as3double> &workarray,
																		as3double                  localtime) override;

		/*!
		 * @brief Function that computes the surface terms in the j-direction in the entire solver, based on the EE.
		 *
		 * @param[in] geometry_container input geometry container.
		 * @param[out] monitordata vector of parameters to monitor.
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 */
		void ComputeSurfaceResidualJDir(CGeometry                 *geometry_container,
															      as3vector1d<as3double>    &monitordata,
															      CPoolMatrixAS3<as3double> &workarray,
																		as3double                  localtime) override;

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


