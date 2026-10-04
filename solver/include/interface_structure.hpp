#pragma once

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "riemann_solver_structure.hpp"
#include <functional> 

// Forward declaration to avoid compiler issues.
class CMultizoneSolver;
class ISolver;


/*!
 * @brief An interface class used for the zone interface specification.
 */
class IInterface
{
	public:

    IInterface(const CConfig               *config_container,
               const CMultizoneSolver      *multizone_solver_container,
               const CInterfaceFacesFamily &interface_family);

	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IInterface(void);

    virtual void ComputeInterfaceResidual(CMultizoneSolver            *multizone_solver_container,
                                          const CInterfaceFacesFamily &family_face,
                                          CFlattenedFaceIndex          face_info,
                                          CPoolMatrixAS3<as3double>   &workarray,
                                          as3double                    localtime) const = 0;

		/*!
		 * @brief Getter function which returns the name of the owner marker.
		 *
		 * @return mIName.
		 */
		const std::string &GetIName(void) const {return mIName;}

		/*!
		 * @brief Getter function which returns the name of the matching marker.
		 *
		 * @return mJName.
		 */
		const std::string &GetJName(void) const {return mJName;}

		/*!
		 * @brief Getter function which returns the ID of the owner zone.
		 *
		 * @return mIZone.
		 */
		unsigned short GetIZone(void) const {return mIZone;}

		/*!
		 * @brief Getter function which returns the ID of the matching zone.
		 *
		 * @return mJZone.
		 */
		unsigned short GetJZone(void) const {return mJZone;}

		/*!
		 * @brief Getter function which returns the number of elements on this interface.
		 *
		 * @return mNFace.
		 */
		size_t GetnFace(void) const {return mNFace;}

	protected:	
		std::string    mIName;    ///< Name of the owner interface marker.
		std::string    mJName;    ///< Name of the matching interface marker.
		unsigned short mIZone;    ///< Zone ID of the owner interface marker.
		unsigned short mJZone;    ///< Zone ID of the matching interface marker.
		EFaceLocation  mIFace;    ///< Face ID of the owner interface marker.
		EFaceLocation  mJFace;    ///< Face ID of the matching interface marker.
		
		size_t         mNFace;    
		unsigned short mNInt1D;   ///< Number of integration points on this marker.	

    unsigned short mNPolyI;
    unsigned short mNPolyJ;

		CMatrixAS3<as3double>           mWInt1D;                  ///< Integration weights on the reference element in 1D.
		std::unique_ptr<ITensorProduct> mITensorProductContainer; ///< Tensor-product container of the iZone.
		std::unique_ptr<ITensorProduct> mJTensorProductContainer; ///< Tensor-product container of the jZone.
		std::unique_ptr<IRiemannSolver> mRiemannSolverContainer;  ///< Riemann solver container.
	
		/*!
		 * @brief Function which returns a function pointer for the interpolation on the (owner) i-face.
		 */
		inline auto GetFuncPointerInterpFaceI(void);

		/*!
		 * @brief Function which returns a function pointer for the interpolation on the (matching) j-face.
		 */
		inline auto GetFuncPointerInterpFaceJ(void);

		/*!
		 * @brief Function which returns a function pointer for the residual computation on the (owner) i-face.
		 */
		inline auto GetFuncPointerResidualFaceI(void);

		/*!
		 * @brief Function which returns a function pointer for the residual computation on the (matching) j-face.
		 */
		inline auto GetFuncPointerResidualFaceJ(void);

	private:

		// Disable default constructor.
		IInterface(void) = delete;
		// Disable default copy constructor.
		IInterface(const IInterface&) = delete;
		// Disable default copy operator.
		IInterface& operator=(IInterface&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for an interface specification based on the (non-linear) Euler equations. 
 */
class CEEInterface : public IInterface
{
  private:
    
    using AComputeResidualFace = std::function<void(const as3double*,
                                                          as3double*,
                                                          as3double*,
                                                          as3double*)>;
    
    using AInterpolateSurface = std::function<void(const as3double*,
                                                         as3double*,
                                                         as3double*,
                                                         as3double*)>;
	
  public:

    CEEInterface(const CConfig               *config_container,
                 const CMultizoneSolver      *multizone_solver_container,
                 const CInterfaceFacesFamily &interface_family);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEEInterface(void) override;

    void ComputeInterfaceResidual(CMultizoneSolver            *multizone_solver_container,
                                  const CInterfaceFacesFamily &family_face,
                                  CFlattenedFaceIndex          face_info,
                                  CPoolMatrixAS3<as3double>   &workarray,
                                  as3double                    localtime) const final;
	protected:

	private:
    unsigned short mNVar;

    AComputeResidualFace mlComputeResidualFaceI;
    AComputeResidualFace mlComputeResidualFaceJ;

    AInterpolateSurface  mlInterpolateSurfaceI;
    AInterpolateSurface  mlInterpolateSurfaceJ;

    void InitializeComputeKernels(void);
};



