#pragma once

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "riemann_solver_structure.hpp"
#include <functional> 

// Forward declaration to avoid compiler issues.
class ISolver;


/*!
 * @brief An interface class used for the zone interface specification.
 */
class IInterface
{
	public:
		
		/*!
		 * @brief Constructor of IInterface, which serves as an interface for the zone interface boundaries.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] param_container input interface parameter container.
		 * @param[in] imarker_container marker container of owner.
		 * @param[in] jmarker_container marker container of matching.
		 * @param[in] solver_container input vector of solver containers.
		 */
		IInterface(CConfig                               *config_container,
				       CGeometry                             *geometry_container,
							 CInterfaceParamMarker                 *param_container,
							 const CMarker                         *imarker_container,
							 const CMarker                         *jmarker_container,
						   as3vector1d<std::unique_ptr<ISolver>> &solver_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IInterface(void);


		/*!
		 * @brief Pure virtual function that computes the residual on the zone interface marker. Must be overridden.
		 *
		 * @param[in] solver_container input vector of solver containers.
		 * @param[in] workarray memory for the working array.
		 */
		virtual void ComputeInterfaceResidual(as3vector1d<std::unique_ptr<ISolver>> &solver_container,
				                                  CPoolMatrixAS3<as3double>             &workarray) = 0;

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

	protected:
		std::string    mIName; ///< Name of the owner interface marker.
		std::string    mJName; ///< Name of the matching interface marker.
		unsigned short mIZone; ///< Zone ID of the owner interface marker.
		unsigned short mJZone; ///< Zone ID of the matching interface marker.
		EFaceElement   mIFace; ///< Face ID of the owner interface marker.
		EFaceElement   mJFace; ///< Face ID of the matching interface marker.
		
		unsigned int   mNElem;    ///< Number of elements on this marker.
		unsigned short mNInt1D;   ///< Number of integration points on this marker.
		unsigned short mNVar = 4; ///< Number of working variables.

		as3vector1d<std::pair<unsigned int, unsigned int>> mIndexElement; ///< Indices of the pair of elements on this interface.

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


		/*!
		 * @brief Function that processes the marker pairs, such that their interface boundaries match.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] imarker_container marker container of owner.
		 * @param[in] jmarker_container marker container of matching.
		 * @param[in] param_container input interface parameter container.
		 */
		void ProcessMatchingMarkers(CConfig               *config_container,
				                        CGeometry             *geometry_container,
																const CMarker         *owner_marker,
																const CMarker         *match_marker,
																CInterfaceParamMarker *param_container);

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
	public:

		/*!
		 * @brief Constructor of CEEInterface, which initializes an Euler-equations interface boundary.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] param_container input interface parameter container.
		 * @param[in] imarker_container marker container of owner.
		 * @param[in] jmarker_container marker container of matching.
		 * @param[in] solver_container input vector of solver containers.
		 */
		CEEInterface(CConfig                               *config_container,
				         CGeometry                             *geometry_container,
								 CInterfaceParamMarker                 *param_container,
								 const CMarker                         *imarker_container,
								 const CMarker                         *jmarker_container,
						     as3vector1d<std::unique_ptr<ISolver>> &solver_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEEInterface(void) override;

		/*!
		 * @brief Function that computes the residual on the zone interface marker.
		 *
		 * @param[in] solver_container input vector of solver containers.
		 * @param[in] workarray memory for the working array.
		 */
		void ComputeInterfaceResidual(as3vector1d<std::unique_ptr<ISolver>> &solver_container,
				                          CPoolMatrixAS3<as3double>             &workarray) override;


	protected:

	private:

};

// Definitions of the inlined functions.
#include "interface_structure.inl"




