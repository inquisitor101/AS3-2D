#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"


/*!
 * @brief An interface class used for the spatial discretization specification.
 */
class ISpatial
{
	public:
		
		/*!
		 * @brief Constructor of ISpatial, which serves as an interface for the spatial discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		ISpatial(CConfig   *config_container,
				     CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ISpatial(void);



	protected:

	private:
		// Disable default constructor.
		ISpatial(void) = delete;
		// Disable default copy constructor.
		ISpatial(const ISpatial&) = delete;
		// Disable default copy operator.
		ISpatial& operator=(ISpatial&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a spatial discretization purely based on the Euler equations. 
 */
class CEESpatial : public ISpatial
{
	public:

		/*!
		 * @brief Constructor of CEESpatial, which initializes a pure EE spatial class. 
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 */
		CEESpatial(CConfig   *config_container,
				       CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CEESpatial(void) override;

	protected:

	private:

};

