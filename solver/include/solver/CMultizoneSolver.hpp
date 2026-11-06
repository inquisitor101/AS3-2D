#pragma once

/*!
 * @file
 * @brief Multizone ownership of zone solvers, interfaces and physical boundaries.
 */


#include <cstddef>
#include <memory>

#include "common/AS3Options.hpp"

class CConfig;
class CMultizoneGeometry;
class ISinglezoneSolver;
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
     * @brief Destructor that releases the owned solvers, interfaces and boundaries.
     */
    ~CMultizoneSolver(void);

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
    const ISinglezoneSolver* GetSinglezoneSolver(std::size_t iSolver) const 
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
    ISinglezoneSolver* GetSinglezoneSolver(std::size_t iSolver)
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
    as3vector1d<std::unique_ptr<ISinglezoneSolver>> mMultizoneSolverContainer;     ///< Owned solvers indexed by zone.
    as3vector1d<std::unique_ptr<IInterface>> mMultizoneInterfaceContainer;  ///< Owned interface residual evaluators.
    as3vector1d<std::unique_ptr<IBoundary>>  mMultizoneBoundaryContainer;   ///< Owned physical boundary residual evaluators.

    as3vector1d<std::size_t> mIndexInterfaceIFace;  ///< index for interfaces whose i-face is of type: IFace.
    as3vector1d<std::size_t> mIndexInterfaceJFace;  ///< index for interfaces whose i-face is of type: JFace.

    as3vector1d<std::size_t> mIndexBoundaryIFace;  ///< index for boundaries whose face is of type: IFace.
    as3vector1d<std::size_t> mIndexBoundaryJFace;  ///< index for boundaries whose face is of type: JFace.


};
