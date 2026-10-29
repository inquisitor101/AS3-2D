#pragma once

/*!
 * @file
 * @brief Multizone solver ownership and zone-local Euler residual operations.
 */


#include <cstddef>
#include <memory>

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
class IBoundary;
class IInterface;


/*!
 * @brief Ownership of zone solvers, interfaces and physical boundary conditions.
 */
class CMultizoneSolver
{
  public:
    

    /*!
     * @brief Constructor that initializes CMultizoneSolver from the specified data.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] multizone_geometry_container geometry of all grid zones.
     */
    CMultizoneSolver(const CConfig            *config_container,
                     const CMultizoneGeometry *multizone_geometry_container);

    /*!
     * @brief Function that constructs interface evaluators and direction-local batch mappings.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] multizone_geometry_container geometry of all grid zones.
     */
    void InitializeInterfaces(const CConfig            *config_container,
                              const CMultizoneGeometry *multizone_geometry_container);

    /*!
     * @brief Function that constructs physical boundary evaluators and direction-local batch mappings.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] multizone_geometry_container geometry of all grid zones.
     */
    void InitializeBoundaries(const CConfig            *config_container,
                              const CMultizoneGeometry *multizone_geometry_container);

    

    /*!
     * @brief Getter function which returns the solver associated with the requested zone.
     *
     * @param[in] iSolver zero-based zone-solver index.
     *
     * @return solver associated with the requested zone.
     */
    const ISolver* GetSinglezoneSolver(std::size_t iSolver) const 
    { 
      return mMultizoneSolverContainer[iSolver].get();
    }
    

    /*!
     * @brief Getter function which returns the solver associated with the requested zone.
     *
     * @param[in] iSolver zero-based zone-solver index.
     *
     * @return solver associated with the requested zone.
     */
    ISolver* GetSinglezoneSolver(std::size_t iSolver)
    {
      return mMultizoneSolverContainer[iSolver].get();
    }

    /*!
     * @brief Getter function which returns the interface evaluators owned by the multizone solver.
     *
     * @return interface evaluators owned by the multizone solver.
     */
    const auto& GetMultizoneInterface(void) const
    {
      return mMultizoneInterfaceContainer;
    }
    

    /*!
     * @brief Getter function which returns one interface evaluator.
     *
     * @param[in] iInterface index in the interface container.
     *
     * @return one interface evaluator.
     */
    const auto* GetInterface(std::size_t iInterface) const 
    {
      return mMultizoneInterfaceContainer[iInterface].get();
    }
    

    /*!
     * @brief Getter function which returns one physical boundary evaluator.
     *
     * @param[in] iBoundary index in the boundary-condition container.
     *
     * @return one physical boundary evaluator.
     */
    const auto* GetBoundary(std::size_t iBoundary) const 
    {
      return mMultizoneBoundaryContainer[iBoundary].get();
    }

    /*!
     * @brief Getter function which returns the number of zones.
     *
     * @return number of zones.
     */
    std::size_t GetnZone(void) const { return mNZone; }

    /*!
     * @brief Getter function which returns the interface-container index associated with the I-direction
     *   interface batch.
     *
     * @param[in] iBatch batch index local to the selected face category.
     *
     * @return interface-container index associated with the I-direction interface batch.
     */
    std::size_t GetIndexInterfaceIFace(std::size_t iBatch) const
    {
#if DEBUG
      if(iBatch >= mIndexInterfaceIFace.size()) ERROR("I-face batch index is out of range.");
#endif
      return mIndexInterfaceIFace[iBatch];
    }
    

    /*!
     * @brief Getter function which returns the interface-container index associated with the J-direction
     *   interface batch.
     *
     * @param[in] iBatch batch index local to the selected face category.
     *
     * @return interface-container index associated with the J-direction interface batch.
     */
    std::size_t GetIndexInterfaceJFace(std::size_t iBatch) const
    {
#if DEBUG
      if(iBatch >= mIndexInterfaceJFace.size()) ERROR("J-face batch index is out of range.");
#endif
      return mIndexInterfaceJFace[iBatch];
    } 

    /*!
     * @brief Getter function which returns the boundary-container index associated with the I-direction
     *   boundary batch.
     *
     * @param[in] iBatch batch index local to the selected face category.
     *
     * @return boundary-container index associated with the I-direction boundary batch.
     */
    std::size_t GetIndexBoundaryIFace(std::size_t iBatch) const
    {
#if DEBUG
      if(iBatch >= mIndexBoundaryIFace.size()) ERROR("I-face batch index is out of range.");
#endif
      return mIndexBoundaryIFace[iBatch];
    }
    

    /*!
     * @brief Getter function which returns the boundary-container index associated with the J-direction
     *   boundary batch.
     *
     * @param[in] iBatch batch index local to the selected face category.
     *
     * @return boundary-container index associated with the J-direction boundary batch.
     */
    std::size_t GetIndexBoundaryJFace(std::size_t iBatch) const
    {
#if DEBUG
      if(iBatch >= mIndexBoundaryJFace.size()) ERROR("J-face batch index is out of range.");
#endif
      return mIndexBoundaryJFace[iBatch];
    } 

  protected:

  private:
    std::size_t mNZone;                                                     ///< Number of grid zones.
    as3vector1d<std::unique_ptr<ISolver>>    mMultizoneSolverContainer;     ///< Owned solvers indexed by zone.
    as3vector1d<std::unique_ptr<IInterface>> mMultizoneInterfaceContainer;  ///< Owned interface residual evaluators.
    as3vector1d<std::unique_ptr<IBoundary>>  mMultizoneBoundaryContainer;   ///< Owned physical boundary residual evaluators.

    as3vector1d<std::size_t> mIndexInterfaceIFace;  ///< index for interfaces whose i-face is of type: IFace.
    as3vector1d<std::size_t> mIndexInterfaceJFace;  ///< index for interfaces whose i-face is of type: JFace.

    as3vector1d<std::size_t> mIndexBoundaryIFace;  ///< index for boundaries whose face is of type: IFace.
    as3vector1d<std::size_t> mIndexBoundaryJFace;  ///< index for boundaries whose face is of type: JFace.


};



//-----------------------------------------------------------------------------------



/*!
 * @brief An interface class used for the solver specification.
 */
class ISolver
{
	public:
		

		/*!
		 * @brief Constructor that initializes ISolver from the specified data.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in] nVar number of solution variables.
		 */
		ISolver(const CConfig             *config_container,
				    const CSinglezoneGeometry *singlezone_geometry_container,
            std::size_t             nVar);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		virtual ~ISolver(void);

		/*!
		 * @brief Getter function which returns the governing-equation type.
		 *
		 * @return implemented governing-equation type.
		 */
		virtual ETypeSolver GetTypeSolver(void) const = 0;

		/*!
		 * @brief Pure virtual function that initializes the physical elements. Must be implemented in a
		 *   derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input multizone geometry container.
		 */
		virtual void InitPhysicalElements(const CConfig            *config_container,
				                              const CMultizoneGeometry *multizone_geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that computes the volume terms for a given element.
		 *
		 * @param[in] singlezone_container geometry of the current grid zone.
		 * @param[in,out] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 *
		 * @note Overwrites the element volume residual. Surface gathering and inverse mass application
		 *   are performed by the element residual operation.
		 */
		virtual void ComputeVolumeResidual(CSinglezoneGeometry       *singlezone_container,
																			 CPoolMatrixAS3<as3double> &workarray,
																			 as3double                  localtime,
																			 std::size_t                     iElem) = 0;

		/*!
		 * @brief Function that computes one same-zone I-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemL zone-local minus-side element adjacent to the I-face.
		 */
		virtual void ComputeSurfaceResidualIDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  std::size_t                     iElemL) = 0;
		

    /*!
     * @brief Function that computes one same-zone J-face flux and writes the two adjacent face residual
     *   buffers.
     *
     * @param[in] singlezone_container geometry of the current zone.
     * @param[in,out] workarray workspace pool owned by the calling thread.
     * @param[in] localtime physical time at which the residual is evaluated.
     * @param[in] iElemB zone-local minus-side element adjacent to the J-face.
     */
    virtual void ComputeSurfaceResidualJDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  std::size_t                     iElemB) = 0;


		/*!
		 * @brief Getter function which returns the number of solution variables.
		 *
		 * @return number of solution variables.
		 */
		virtual std::size_t GetnVar(void) const = 0;

		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		std::size_t GetZoneID(void) const {return mZoneID;}

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
		CPhysicalElement *GetPhysicalElement(std::size_t index) const {return mPhysicalElementContainer[index].get();}


	protected:
		const std::size_t                              mZoneID;                    ///< Zone ID of this container.
		std::unique_ptr<ITensorProduct>                mTensorProductContainer;    ///< Tensor product container.
		std::unique_ptr<IRiemannSolver>                mRiemannSolverContainer;    ///< Riemann solver container.
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
class CEESolver final: public ISolver
{
	public:

		/*!
		 * @brief Constructor that initializes CEESolver from the specified data.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 */
		CEESolver(const CConfig             *config_container,
				      const CSinglezoneGeometry *singlezone_geometry_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CEESolver(void) override;

		/*!
		 * @brief Getter function which returns the governing-equation type.
		 *
		 * @return implemented governing-equation type.
		 */
		ETypeSolver GetTypeSolver(void) const override {return ETypeSolver::EE;}

		/*!
		 * @brief Function that initializes the physical elements.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 */
		void InitPhysicalElements(const CConfig            *config_container,
				                      const CMultizoneGeometry *multizone_geometry_container) override;

		/*!
		 * @brief Function that computes the volume terms for a given element, based on the EE.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current grid zone.
		 * @param[in,out] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 *
		 * @note Overwrites the element volume residual. Surface gathering and inverse mass application
		 *   are performed by the element residual operation.
		 */
		void ComputeVolumeResidual(CSinglezoneGeometry       *singlezone_geometry_container,
				                       CPoolMatrixAS3<as3double> &workarray,
															 as3double                  localtime,
															 std::size_t                     iElem) override;


		/*!
		 * @brief Function that computes one same-zone I-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemL zone-local minus-side element adjacent to the I-face.
		 */
		void ComputeSurfaceResidualIDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  std::size_t                     iElemL) override final;

		/*!
		 * @brief Function that computes one same-zone J-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemB zone-local minus-side element adjacent to the J-face.
		 */
		void ComputeSurfaceResidualJDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  std::size_t                     iElemB) override final;

		/*!
		 * @brief Getter function which returns the number of solution variables.
		 *
		 * @return mNVar.
		 */
		std::size_t GetnVar(void) const override {return mNVar;}

	protected:

	private:
		inline static constexpr std::size_t mNVar = 4; ///< Number of working variables
};



