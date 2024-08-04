#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


// Forward declaration to avoid compiler issues.
class ISolver;

/*!
 * @brief An interface class used for writing VTK files. 
 */
class IFileVTK
{
	public:
	
		/*!
		 * @brief Constructor of IFileVTK, which serves as an interface for VTK-type formats.
		 */
		IFileVTK(CConfig   *config_container,
				     CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IFileVTK(void);

		/*!
		 * @brief Pure virtual function that writes VTK data a file. 
		 * Must be implemented by a derived class.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container reference to the solver container.
		 */
		virtual void WriteFileVTK(CConfig                               *config_container,
															CGeometry                             *geometry_container,
															as3vector1d<std::unique_ptr<ISolver>> &solver_container) = 0;

	protected:
		as3vector1d<std::string> mVariableNames;     ///< Variable names for writing.
		unsigned long            mFileNumber;        ///< Written file index.

		bool                     mWriteDensity;      ///< Option for writing a density variable.
		bool                     mWriteMomentum;     ///< Option for writing a momentum variable.
		bool                     mWriteTotalEnergy;  ///< Option for writing a total energy variable.
		bool                     mWritePressure;     ///< Option for writing a pressure variable.
		bool                     mWriteVelocity;     ///< Option for writing a velocity variable.
		bool                     mWriteVorticity;    ///< Option for writing a vorticity variable.
		bool                     mWriteMach;         ///< Option for writing a Mach number variable.
		bool                     mWriteTemperature;  ///< Option for writing a temperature variable.
		bool                     mWriteEntropy;      ///< Option for writing a specific entropy variable.

	private:
		// Disable default constructor.
		IFileVTK(void) = delete;
	
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class used for writing VTK files in binary and legacy format. 
 */
class CLegacyBinaryVTK final : public IFileVTK
{
	public:
	
		/*!
		 * @brief Constructor of CLegacyBinaryVTK, which is reponsible for VTK in binary and legacy format.
		 */
		CLegacyBinaryVTK(CConfig   *config_container,
				             CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CLegacyBinaryVTK(void) final;

		/*!
		 * @brief Function that writes VTK data to a file using a binary and legacy format.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 * @param[in] solver_container reference to the solver container.
		 */
		void WriteFileVTK(CConfig                               *config_container,
											CGeometry                             *geometry_container,
											as3vector1d<std::unique_ptr<ISolver>> &solver_container) final;

	protected:

	private:
		bool mBigEndian; ///< Whether this machine uses big endian format.
	
		// Disable default constructor.
		CLegacyBinaryVTK(void) = delete;
	
		/*!
		 * @brief Function that computes and stores the data required for writing a binary VTK file.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to geometry container.
		 * @param[in] solver_container reference to the solver container.
		 * @param[in] vars_buf reference to the variables buffer data.
		 */
	  void DetermineVisualizationData(CConfig                               *config_container,
				                            CGeometry                             *geometry_container,
				                            as3vector1d<std::unique_ptr<ISolver>> &solver_container,
																		as3vector2d<float>                    &vars_buf);
};

