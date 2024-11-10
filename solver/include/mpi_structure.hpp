#pragma once 

#include <mpi.h>
#include "config_structure.hpp"

/*!
 * @brief A class for the MPI (distributed memory) parallelization strategy.
 */
class CMPI
{
	public:

		/*!
		 * @brief Constructor of CMPI, which initializes the MPI class.
		 */
		CMPI(CConfig *config_container);

		/*!
		 * @brief Destructor, which frees any allocated data + ends the processes.
		 */
		~CMPI(void);

		// Initialize the MPI implementation.
		static inline void Init(int *argc, char ***argv)
		{
			MPI_Init(argc, argv);
			MPI_Comm_rank(MPI_COMM_WORLD, &mRankI);
			MPI_Comm_size(MPI_COMM_WORLD, &mNRank);
		}

		static inline int GetRankI(void) { return mRankI; }
		static inline int GetnRank(void) { return mNRank; }

	private:
		static inline int mRankI = 0;  ///< Current rank ID of this processor. 
		static inline int mNRank = 1;  ///< Overall total number of ranks specified.


		// Disable default constructor.
		CMPI(void)  = delete;		
		// Disable default copy constructor.
		CMPI(const CMPI&) = delete;
		// Disable default copy operator.
		CMPI& operator=(CMPI&) = delete;	
};


