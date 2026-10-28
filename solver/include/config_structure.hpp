#pragma once 

/*!
 * @file
 * @brief Simulation configuration, parameter lookup and option validation.
 */


#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "option_structure.hpp"
#include "input_structure.hpp"
#include "param_structure.hpp"


/*!
 * @brief A class used for storing the user-specified configuration options.
 */
class CConfig
{
	public:

		/*!
		 * @brief Constructor of CConfig, which is responsible for the configuration options.
		 *
		 * @param[in] filename input configuration file.
		 */
		CConfig(const char *filename);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CConfig(void);


		/*!
		 * @brief Getter function which returns the number of zones.
		 *
		 * @return mNZone
		 */ 
		std::size_t GetnZone(void)                                  const {return mNZone;}

		/*!
		 * @brief Getter function which returns the configured CFL stability number.
		 *
		 * @return mCFL
		 */ 
		as3double GetCFL(void)                                      const {return mCFL;}
	

		/*!
		 * @brief Getter function which returns the value of mOutputVisFilename.
		 *
		 * @return mOutputVisFilename
		 */	
		const char *GetOutputVisFilename(void)                      const {return mOutputVisFilename.c_str();}

		/*!
		 * @brief Getter function which returns the value of mOutputSolFilename.
		 *
		 * @return mOutputSolFilename
		 */	
		const char *GetOutputSolFilename(void)                      const {return mOutputSolFilename.c_str();}

		/*!
		 * @brief Getter function which returns the value of mMeshFormat.
		 *
		 * @return mMeshFormat
		 */	
		EMeshFormat GetMeshFormat(void)                             const {return mMeshFormat;}

		/*!
		 * @brief Getter function which returns the value of mInputGridFormat.
		 *
		 * @return mInputGridFormat
		 */	
		EFormatFile GetInputGridFormat(void)                        const {return mInputGridFormat;}

		/*!
		 * @brief Getter function which returns the value of mOutputVisFormat.
		 *
		 * @return mOutputVisFormat
		 */	
		EVisualFormat GetOutputVisFormat(void)                      const {return mOutputVisFormat;}

		/*!
		 * @brief Getter function which returns the value of mTemporalScheme.
		 *
		 * @return mTemporalScheme
		 */	
		ETemporalScheme GetTemporalScheme(void)                     const {return mTemporalScheme;}

		/*!
		 * @brief Getter function which returns the numerical flux method.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeRiemannSolver[iZone]
		 */   
		ETypeRiemannSolver GetTypeRiemannSolver(std::size_t iZone)  const {return mTypeRiemannSolver[iZone];}

		/*!
		 * @brief Getter function which returns the governing-equation type.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeSolver[iZone]
		 */   
		ETypeSolver GetTypeSolver(std::size_t iZone)                const {return mTypeSolver[iZone];}

		/*!
		 * @brief Getter function which returns the value of mTypeBufferLayer, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeBufferLayer[iZone]
		 */   
		ETypeBufferLayer GetTypeBufferLayer(std::size_t iZone)      const {return mTypeBufferLayer[iZone];}

		/*!
		 * @brief Getter function which returns the value of mTypeIC.
		 *
		 * @return mTypeIC
		 */   
		ETypeIC GetTypeIC(void)                                     const {return mTypeIC;}

		/*!
		 * @brief Getter function which returns the visualization variables to write.
		 *
		 * @return mWriteVisVar
		 */
		as3vector1d<EWriteVariable> GetWriteVisVar(void)            const {return mWriteVisVar;}

		/*!
		 * @brief Getter function which returns the simulation start time.
		 *
		 * @return mStartTime
		 */
		as3double GetStartTime(void)                                const {return mStartTime;}

		/*!
		 * @brief Getter function which returns the time step specified.
		 *
		 * @return mTimeStep
		 */
		as3double GetTimeStep(void)                                 const {return mTimeStep;}

		/*!
		 * @brief Getter function which returns the maximum number of (temporal) iterations.
		 *
		 * @return mMaxIterTime
		 */
		std::size_t GetMaxIterTime(void)                            const {return mMaxIterTime;}

		/*!
		 * @brief Getter function which returns the writing frequency of the visualization files.
		 *
		 * @return mWriteVisFreq
		 */
		std::size_t GetWriteVisFreq(void)                           const {return mWriteVisFreq;}

		/*!
		 * @brief Getter function which returns the value of mGridDirectory.
		 *
		 * @return mGridDirectory
		 */
		std::string GetGridDirectory(void)                          const {return mGridDirectory;}

		/*!
		 * @brief Getter function which returns the vector of marker tags.
		 *
		 * @return mMarkerTag
		 */
		const as3vector1d<std::pair<std::string, ETypeBC>> &GetMarkerTag(void) const {return mMarkerTag;}

		/*!
		 * @brief Getter function which returns a specific value in mDataIC, based on the index.
		 *
		 * @param[in] index index of the item in the data.
		 *
		 * @return mDataIC[index]
		 */   
		as3double GetDataIC(std::size_t index) const 
		{
			// A check is affordable here, as this function is called only in the preprocessing step.
			if( index >= mDataIC.size() ) ERROR("Index exceeds data in IC parameters.");

			// Return the corresponding value.
			return mDataIC[index];
		}

		/*!
		 * @brief Getter function which returns the vector of periodic parameter markers.
		 *
		 * @return mPeriodicParamMarker
		 */
		const as3vector1d<std::unique_ptr<CPeriodicParamMarker>> &GetPeriodicParamMarker(void) const {return mPeriodicParamMarker;}

		/*!
		 * @brief Getter function which returns the vector of boundary parameter markers.
		 *
		 * @return mBoundaryParamMarker
		 */
		const as3vector1d<std::unique_ptr<IBoundaryParamMarker>> &GetBoundaryParamMarker(void)  const {return mBoundaryParamMarker;}

		/*!
		 * @brief Getter function which returns a pointer to the specified boundary marker, if found.
		 *
		 * @param[in] imarker name of the marker to look for.
		 *
		 * @return pointer for the relevant mBoundaryParamMarker element
		 */
		IBoundaryParamMarker *GetBoundaryParamMarker(const std::string &imarker) const
		{
			// Search for the marker by name.
			for( auto& m: mBoundaryParamMarker )
			{
				if( m->mName == imarker ) return m.get();
			}
			// The marker could not be found, issue an error and exit.
			ERROR("Could not find the boundary marker: " + imarker);

			// Return a nullptr, to avoid a compiler issue.
			return nullptr;
		}

		/*!
		 * @brief Getter function which returns the bonus factor for quadrature exactness.
		 *
		 * @return non-negative nb in the nodal-DG target degree (2 + nb)*p; defaults to zero.
		 */
		std::size_t GetBonusQuadratureFactor(void) const {return mBonusQuadratureFactor;}

		/*!
		 * @brief Getter function which returns the configured quadrature rule.
		 *
		 * @return mTypeQuadrature; defaults to Gauss-Legendre.
		 */
		ETypeQuadrature GetTypeQuadrature(void) const {return mTypeQuadrature;}

		/*!
		 * @brief Getter function which returns the configured DG formulation.
		 *
		 * @return mTypeDGMethod; defaults to standard nodal DG.
		 */
		ETypeDGMethod GetTypeDGMethod(void) const {return mTypeDGMethod;}

	protected:

	private:
		EMeshFormat                     mMeshFormat;         ///< Mesh format specified.
		EFormatFile                     mInputGridFormat;    ///< Format of input grid file.
		EVisualFormat                   mOutputVisFormat;    ///< Format of output visualization file.
		ETemporalScheme                 mTemporalScheme;     ///< Type of temporal scheme.
		ETypeIC                         mTypeIC;             ///< Type of initial condition.
		as3double                       mStartTime;          ///< Simulation starting time.
		as3double                       mTimeStep;           ///< Time step.
		as3double                       mCFL;                ///< CFL stability number.
		std::size_t                     mMaxIterTime;        ///< Maximum temporal iterations during the simulation.
		std::size_t                     mWriteVisFreq;       ///< Visualization file writing frequency.
		std::string                     mOutputSolFilename;  ///< Output solution filename.
		std::string                     mOutputVisFilename;  ///< Output visualization filename.
		std::size_t                     mNZone;              ///< Total number of zones.
		as3vector1d<as3double>          mDataIC;             ///< Compressed (floating-type) data for the IC.
    as3vector1d<EWriteVariable>     mWriteVisVar;        ///< Vector of visualization variables to write.
		as3vector1d<ETypeSolver>        mTypeSolver;         ///< Type of PDEs in each zone.
    as3vector1d<ETypeBufferLayer>   mTypeBufferLayer;    ///< Type of buffer layer in each zone.
		as3vector1d<ETypeRiemannSolver> mTypeRiemannSolver;  ///< Type of Riemann solver per each zone.

		std::size_t     mBonusQuadratureFactor; ///< Non-negative nb in the nodal-DG target degree (2 + nb)*p.
		ETypeQuadrature mTypeQuadrature;        ///< Quadrature rule on the reference interval.
		ETypeDGMethod   mTypeDGMethod;          ///< DG formulation.

    std::string mGridDirectory;                          ///< Directory containing the configured AS3 grid file.
		as3vector1d<std::pair<std::string, ETypeBC>>        mMarkerTag;            ///< Vector of marker names and boundary conditions.
		as3vector1d<std::unique_ptr<CPeriodicParamMarker>>  mPeriodicParamMarker;  ///< Container of periodic marker parameters.
		as3vector1d<std::unique_ptr<IBoundaryParamMarker>>  mBoundaryParamMarker;  ///< Container of boundary marker parameters.

    /*!
     * @brief Function that reads the zone configuration options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadZoneConfigurationOptions(const char *filename);

    /*!
     * @brief Function that reads the output options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadOutputOptions(const char *filename);

    /*!
     * @brief Function that reads the temporal options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadTemporalOptions(const char *filename);

    /*!
     * @brief Function that reads the initial condition options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadInitialConditionOptions(const char *filename);

    /*!
     * @brief Function that reads the boundary condition options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadBoundaryConditionOptions(const char *filename);

    /*!
     * @brief Function that reads the spatial discretization options.
     *
     * @param[in] filename input configuration file.
     *
     * @return bool whether the operation succeeded.
     */
    bool ReadSpatialDiscretizationOptions(const char *filename);


		/*!
		 * @brief Function that extracts the parameters relevant to a Gaussian pressure initial
		 *   condition.
		 *
		 * @param[in] filename input configuration file.
		 */
		void IC_GaussianPressure(const char *filename);

		/*!
		 * @brief Function that extracts the parameters relevant to an isentropic vortex initial
		 *   condition.
		 *
		 * @param[in] filename input configuration file.
		 */
		void IC_IsentropicVortex(const char *filename);

		/*!
		 * @brief Function that checks for and extracts the necessary interface/periodic BC information.
		 *
		 * @param[in] filename input configuration file.
		 */
		void ExtractInfoInterfaceBC(const char *filename);

    /*!
     * @brief Function that reads parameters for the supported physical boundary conditions.
     *
     * @param[in] filename input filename.
     */
    void ExtractInfoBoundaryCondition(const char *filename);

		/*!
		 * @brief Function that maps a single string to an enum.
		 *
		 * @tparam Tenum enumeration type associated with the keyword map.
		 *
		 * @param[in] mapper type of map specified.
		 * @param[in] key string name for the key.
		 * @param[in] msg (optional) error message for the keyword.
		 *
		 * @return enum values, according to their associated map.
		 */
		template<class Tenum>
		Tenum GenericScalarMap(const std::map<std::string, Tenum> &mapper,
				                   std::string                        &key,
										       std::string                         msg = "Keyword");

		/*!
		 * @brief Function that maps a vector of strings to an enum.
		 *
		 * @tparam Tenum enumeration type associated with the keyword map.
		 *
		 * @param[in] mapper type of map specified.
		 * @param[in] key vector of string names for the keys.
		 * @param[in] msg (optional) error message for the keywords.
		 *
		 * @return vector of enum values, according to their associated map.
		 */
		template<class Tenum>
		as3vector1d<Tenum> GenericVectorMap(const std::map<std::string, Tenum> &mapper,
				                                as3vector1d<std::string>           &key,
										                    std::string                         msg = "Keyword");

    /*!
     * @brief Function that pads remaining entries in a vector.
     *
     * @tparam TValueType type used for parsed or padded values.
     *
     * @param[in,out] data reference to data that is padded.
     * @param[in] keyword input option string name.
     * @param[in] nExpected expected size of input/output vector.
     * @param[in] nCondition1 input condition on whether to pad.
     * @param[in] nCondition2 input condition on whether to pad.
     * @param[in] nCondition3 input condition on whether to pad.
     */
    template<class TValueType>
    void PadEntriesVectorData(std::vector<TValueType> &data,
                              std::string        			 keyword,
                              std::size_t           nExpected,
                              std::size_t           nCondition1 = 0,
                              std::size_t           nCondition2 = 0,
                              std::size_t           nCondition3 = 0);
    

    /*!
     * @brief Function that pads remaining entries in a vector based on default reference data.
     *
     * @tparam TValueType type used for parsed or padded values.
     *
     * @param[in,out] data reference to data that is padded.
     * @param[in] reference reference to default value in padding.
     * @param[in] keyword input option string name.
     */
    template<class TValueType>
    void PadEntriesVectorDefaultData(std::vector<TValueType> &data,
                                     std::vector<TValueType> &reference,
                                     std::string        			keyword);


		// Disable default constructor.
		CConfig(void) = delete;
		// Disable default copy constructor.
		CConfig(const CConfig&) = delete;
		// Disable default copy operator.
		CConfig& operator=(CConfig&) = delete;	
};

//-----------------------------------------------------------------------------------

// Definitions of the templated functions.
#include "config_structure.inl"


