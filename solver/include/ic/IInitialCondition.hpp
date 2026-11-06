#pragma once

/*!
 * @file
 * @brief Interface for initializing a zone solution.
 */



class CConfig;
class CSinglezoneGeometry;
class ISinglezoneSolver;


/*!
 * @brief An interface for initializing the conservative solution in one zone.
 */
class IInitialCondition
{
	public:

		/*!
		 * @brief Constructor of IInitialCondition, which serves as an interface for the initial
		 *   condition.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		IInitialCondition(const CConfig *config_container);

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IInitialCondition(void);


		/*!
		 * @brief A Pure virtual function that initializes the solution over a single zone. Must be
		 *   overridden.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry zone geometry container.
		 * @param[in,out] solver_container solver container.
		 */
		virtual void InitializeSolution(CConfig             *config_container,
				                            CSinglezoneGeometry *singlezone_geometry,
				                            ISinglezoneSolver  *solver_container) = 0;

	protected:

	private:

		// Disable default constructor.
		IInitialCondition(void) = delete;
		// Disable default copy constructor.
		IInitialCondition(const IInitialCondition&) = delete;
		// Disable default copy operator.
		IInitialCondition& operator=(IInitialCondition&) = delete;	
};
