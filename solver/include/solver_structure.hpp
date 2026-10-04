#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "tensor_structure.hpp"
#include "factory_structure.hpp"
#include "boundary_structure.hpp"
#include "monitoring_structure.hpp"
#include "riemann_solver_structure.hpp"
#include "standard_element_structure.hpp"
#include "physical_element_structure.hpp"
#include "interface_structure.hpp"

// Forward declaration to avoid compiler problems.
class ISolver;
class IInterface;


class CMultizoneSolver
{
  public:
    
    CMultizoneSolver(const CConfig            *config_container,
                     const CMultizoneGeometry *multizone_geometry_container);

    void InitializeInterfaces(const CConfig            *config_container,
                              const CMultizoneGeometry *multizone_geometry_container);

    
    const ISolver* GetSinglezoneSolver(size_t iSolver) const 
    { 
      return mMultizoneSolverContainer[iSolver].get();
    }
    
    ISolver* GetSinglezoneSolver(size_t iSolver)
    {
      return mMultizoneSolverContainer[iSolver].get();
    }

    const auto& GetMultizoneInterface(void) const
    {
      return mMultizoneInterfaceContainer;
    }
    
    const auto* GetSinglezoneInterface(size_t iInterface) const 
    {
      return mMultizoneInterfaceContainer[iInterface].get();
    }

    size_t GetnZone(void) const { return mNZone; }

    size_t GetIndexInterfaceIFace(size_t iFamily) const
    {
#if DEBUG
      if(iFamily >= mIndexInterfaceIFace.size()) ERROR("I-Family index is out of range.");
#endif
      return mIndexInterfaceIFace[iFamily];
    }
    
    size_t GetIndexInterfaceJFace(size_t iFamily) const
    {
#if DEBUG
      if(iFamily >= mIndexInterfaceJFace.size()) ERROR("J-Family index is out of range.");
#endif
      return mIndexInterfaceJFace[iFamily];
    } 

  protected:

  private:
    size_t mNZone;
    as3vector1d<std::unique_ptr<ISolver>>    mMultizoneSolverContainer;
    as3vector1d<std::unique_ptr<IInterface>> mMultizoneInterfaceContainer;

    as3vector1d<size_t> mIndexInterfaceIFace; // index for interfaces whose i-face is is of type: IFace..
    as3vector1d<size_t> mIndexInterfaceJFace; // index for interfaces whose i-face is of type: JFace.
};





/*!
 * @brief An interface class used for the solver specification.
 */
class ISolver
{
	public:
		
		ISolver(const CConfig             *config_container,
				    const CSinglezoneGeometry *singlezone_geometry_container,
            unsigned short             nVar);
		
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
		 * @param[in] multizone_geometry_container input multizone geometry container.
		 */
		virtual void InitPhysicalElements(const CConfig            *config_container,
				                              const CMultizoneGeometry *multizone_geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that initializes the boundary conditions. Must be implemented in a derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		virtual void InitBoundaryConditions(const CConfig            *config_container,
				                                const CMultizoneGeometry *multizone_geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that computes the volume terms for a given element.
		 *
		 * @param[in] singlezone_container geometry of the current grid zone. 
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 */
		virtual void ComputeVolumeResidual(CSinglezoneGeometry       *singlezone_container,
																			 CPoolMatrixAS3<as3double> &workarray,
																			 as3double                  localtime,
																			 size_t                     iElem) = 0;

		virtual void ComputeSurfaceResidualIDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  size_t                     iElemL) = 0;
		
    virtual void ComputeSurfaceResidualJDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  size_t                     iElemB) = 0;


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

		CEESolver(const CConfig             *config_container,
				      const CSinglezoneGeometry *singlezone_geometry_container);
		
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
		void InitPhysicalElements(const CConfig            *config_container,
				                      const CMultizoneGeometry *multizone_geometry_container) override;

		/*!
		 * @brief Function that initializes the boundary conditions. 
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		void InitBoundaryConditions(const CConfig            *config_container,
		                            const CMultizoneGeometry *multizone_geometry_container) override;

		/*!
		 * @brief Function that computes the volume terms for a given element, based on the EE.
		 * 
		 * @param[in] singlezone_geometry_container geometry of the current grid zone. 
		 * @param[in] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 */
		void ComputeVolumeResidual(CSinglezoneGeometry       *singlezone_geometry_container,
				                       CPoolMatrixAS3<as3double> &workarray,
															 as3double                  localtime,
															 size_t                     iElem) override;


		void ComputeSurfaceResidualIDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  size_t                     iElemL) final;

		void ComputeSurfaceResidualJDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  size_t                     iElemB) final;

		/*!
		 * @brief Getter function which returns the number of working variables.
		 *
		 * @return mNVar.
		 */
		unsigned short GetnVar(void) const override {return mNVar;}

	protected:

	private:
		inline static constexpr unsigned short mNVar = 4; ///< Number of working variables
};



