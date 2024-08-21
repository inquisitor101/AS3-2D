#pragma once 

#include "option_structure.hpp"
#include "input_structure.hpp"
#include <algorithm>


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
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CConfig(void);

		/*!
		 * @brief Function that determines the boundary condition for a given marker.
		 * 
		 * @param[in] marker marker name.
		 *
		 * @return boundary condition
		 */
		ETypeBC DetermineMarkerBC(std::string &marker);

		/*!
		 * @brief Function that finds the name of the matching periodic boundary marker.
		 *
		 * @param[in] imarker name of the current marker.
		 * 
		 * @return name of the matching marker
		 */
		std::string FindMatchingPeriodicMarker(const std::string &imarker);

		/*!
		 * @brief Getter function which returns the value of nZone.
		 *
		 * @return nZone
		 */ 
		unsigned short GetnZone(void)                               const {return mNZone;}
	    
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
		 * @brief Getter function which returns the value of mTypeDOF, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeDOF[iZone]
		 */   
		ETypeDOF GetTypeDOF(unsigned short iZone)                   const {return mTypeDOF[iZone];}

		/*!
		 * @brief Getter function which returns the value of mTypeSolver, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeSolver[iZone]
		 */   
		ETypeSolver GetTypeSolver(unsigned short iZone)             const {return mTypeSolver[iZone];}

		/*!
		 * @brief Getter function which returns the value of mTypeBufferLayer, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mTypeBufferLayer[iZone]
		 */   
		ETypeBufferLayer GetTypeBufferLayer(unsigned short iZone)   const {return mTypeBufferLayer[iZone];}

		/*!
		 * @brief Getter function which returns the value of mNPolyGrid, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mNPolyGrid[iZone]
		 */   
		unsigned short GetnPolyGrid(unsigned short iZone)           const {return mNPolyGrid[iZone];}

		/*!
		 * @brief Getter function which returns the value of mNPolySol, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mNPolySol[iZone]
		 */   
		unsigned short GetnPolySol(unsigned short iZone)            const {return mNPolySol[iZone];}

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
		 * @brief Getter function which returns the simulation end time.
		 *
		 * @return mFinalTime
		 */
		as3double GetFinalTime(void)                                const {return mFinalTime;}

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
		size_t GetMaxIterTime(void)                                 const {return mMaxIterTime;}

		/*!
		 * @brief Getter function which returns the value of mZoneGridFilename, per input zone.
		 *
		 * @param[in] iZone zone ID.
		 *
		 * @return mZoneGridFilename[iZone]
		 */
		std::string GetZoneGridFilename(unsigned short iZone)       const {return mZoneGridFilename[iZone];}

		/*!
		 * @brief Getter function which returns the periodic marker names.
		 *
		 * @return mMarkerNamePeriodic
		 */
		const as3vector2d<std::string> &GetMarkerNamePeriodic(void) const {return mMarkerNamePeriodic;}

		/*!
		 * @brief Getter function which returns a specific value in mDataIC, based on the index.
		 *
		 * @param[in] index index of the item in the data.
		 *
		 * @return mDataIC[index]
		 */   
		as3double GetDataIC(size_t index) const 
		{
			// A check is affordable here, as this function is called only in the preprocessing step.
			if( index >= mDataIC.size() ) ERROR("Index exceeds data in IC parameters.");

			// Return the corresponding value.
			return mDataIC[index];
		}

	protected:

	private:
		EMeshFormat     mMeshFormat;                    ///< Mesh format specified.
		EFormatFile     mInputGridFormat;               ///< Format of input grid file.
		EVisualFormat   mOutputVisFormat;               ///< Format of output visualization file.
		ETemporalScheme mTemporalScheme;                ///< Type of temporal scheme.
		ETypeIC         mTypeIC;                        ///< Type of initial condition.

		as3double       mStartTime;                     ///< Simulation starting time.
		as3double       mFinalTime;                     ///< Simulation ending time.
		as3double       mTimeStep;                      ///< Time step.
		size_t          mMaxIterTime;                   ///< Maximum temporal iterations during the simulation.

		as3vector1d<as3double>      mDataIC;            ///< Compressed (floating-type) data for the IC.
    as3vector1d<EWriteVariable> mWriteVisVar;       ///< Vector of visualization variables to write.

		std::string     mOutputSolFilename;             ///< Output solution filename.
		std::string     mOutputVisFilename;             ///< Output visualization filename.
		unsigned short  mNZone;		    		              ///< Total number of zones.
		
		as3vector1d<unsigned short>   mZoneIndex;       ///< Indices for each zone.
		as3vector1d<ETypeSolver>      mTypeSolver;      ///< Type of PDEs in each zone.
    as3vector1d<ETypeBufferLayer> mTypeBufferLayer; ///< Type of buffer layer in each zone.
		as3vector1d<ETypeDOF>         mTypeDOF;         ///< Type of DOFs in each zone.
    as3vector1d<unsigned short>   mNPolyGrid;       ///< Polynomial order per zone based on the grid.
		as3vector1d<unsigned short>   mNPolySol;        ///< Polynomial order per zone based on solution.

		as3vector1d<std::string>                      mZoneGridFilename;    ///< Grid filename per each zone.
		as3vector2d<std::string>                      mMarkerNamePeriodic;  ///< Container of marker names for periodic BCs.
		as3vector1d<std::pair<std::string, ETypeBC>>  mBoundaryMarkers;     ///< Container of pairs of boundary markers: name and type.

    /*!
		 * @brief Function that reads the zone connectivity options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadZoneConnectivityOptions(const char *filename);

    /*!
		 * @brief Function that reads the zone configuration options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadZoneConfigurationOptions(const char *filename);

    /*!
		 * @brief Function that reads the output options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadOutputOptions(const char *filename);

    /*!
		 * @brief Function that reads the temporal options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadTemporalOptions(const char *filename);

    /*!
		 * @brief Function that reads the initial condition options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadInitialConditionOptions(const char *filename);

    /*!
		 * @brief Function that reads the boundary condition options.
		 *
		 * @param[in] filename input configuration file.
		 *
		 * @return bool whether the operation succeedeed.
		 */
    bool ReadBoundaryConditionOptions(const char *filename);

		/*!
		 * @brief Function that extracts the grid zone filenames.
		 *
		 * @param[in] filename input configuration file.
		 */
    void ExtractZoneGridFiles(const char *filename);

		/*!
		 * @brief Function that extracts the parameters relevant to a Gaussian pressure initial condition.
		 *
		 * @param[in] filename input configuration file.
		 */
		void IC_GaussianPressure(const char *filename);

		/*!
		 * @brief Function that extracts the parameters relevant to an isentropic vortex initial condition.
		 *
		 * @param[in] filename input configuration file.
		 */
		void IC_IsentropicVortex(const char *filename);

		/*!
		 * @brief Function that checks for and extracts the necessary periodic BC information.
		 *
		 * @param[in] filename input configuration file.
		 */
		void ExtractInfoPeriodicBC(const char *filename);

		/*!
		 * @brief Function that maps a single string to an enum.
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
     * @brief Function that does a consistency check for zone configuration.
     */
		void ConsistencyCheckZoneConfiguration(void);

    /*!
     * @brief Function that pads remaining entries in a vector.
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
                              unsigned short           nExpected,
                              unsigned short           nCondition1 = 0,
                              unsigned short           nCondition2 = 0,
                              unsigned short           nCondition3 = 0);
    
    
    /*!
     * @brief Function that pads remaining entries in a vector based on default reference data.
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


