C++ Coding Style in AS3-2D
===============

The general rule we follow is:

1. Global variables are written in capital, e.g. `GAS_CONSTANT`.
2. Vector containers are aliased for readability, e.g. `std::vector<std::vector<T>>` is written as `as3vector2d<T>`. 
3. The working precision of the solution is typdef'd into `as3double`, by default it uses `double-precision`.
4. Classes (or structures) can be two types: interfaces or regular classes. Interfaces start with `I` and regular classes with `C`, followed by their name. For instance, the solver class has an interface class `ISolver` and an implementation of type `CEESolver`.
5. Namespaces always start with `N`, followed by their name. For instance, a linear algebra namespace is: `NLinearAlgebra`.
6. Global functionalities should be defined inside namespaces or classes. If they require no member variables, it's probably better to define them inside namespaces. Otherwise, define them as static member functions in classes, e.g. see `CGenericFactory`.
7. Functions should always start with a capital letter, using the *pascal case* style. For instance, `ComputeResidual`.
8. Member variables in a class/struct are always prefixed with a small-case `m`. For instance, `mLagrangeInt1D`. 
9. Addition of user-defined features in the code require adjusting the `.cfg` file. Thus, ensure your addition is well-commented and placed in its relevant section in the `.cfg` file.
10. For core functionalities, please include approriate tests. These are located in `tests/` and utilize *Google Test* (gtest).
11. Header files use the extension `.hpp` while template implementations use `.inl`. Both are placed in the directory `include/`.
12. Source files use the extension `.cpp` and are situated in the directory `src/`.
13. Header files are documented in line with doxygen-supported comments (see below).
14. Functions in source files are documented between their input arguments and the implementation (i.e. between `)` and `{`). See below.
15. Pay close attention to the vertical alignment of the input arguments to functions, as it will complicate readability. 
16. In the root directory of **AS3-2D**, there exists a `.editorconfig` file, which sets the relevant editor-related conventions. For instance, tab indentation uses 2 characters.


### Example header file.

``solver_structure.hpp``

```C++
/*!
 * @brief A class for a solver specification based on the (non-linear) Euler equations. 
 */
class CEESolver : public ISolver
{
  public:
  
    /*!
     * @brief Constructor of CEESolver, which initializes an Euler equations solver.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] geometry_container input geometry container.
     * @param[in] iZone zone ID of this container.
     */
    CEESolver(CConfig       *config_container,
              CGeometry     *geometry_container,
              unsigned short iZone);
    
    /*!
     * @brief Destructor, which frees any allocated memory.
     */
    ~CEESolver(void) override;
    
    /*!
     * @brief Function that initializes the physical elements.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] geometry_container input geometry container.
     */
    void InitPhysicalElements(CConfig   *config_container,
                              CGeometry *geometry_container) override;
    
    /*!
     * @brief Function that computes the volume terms over a single element, based on the EE equations.
     * 
     * @param[in] iElem element index.
     * @param[in] localtime current physical time.
     * @param[out] monitordata vector of parameters to monitor.
     * @param[out] workarray reference to the work array.
     */
    void ComputeVolumeResidual(size_t                     iElem,
                               as3double                  localtime,
                               as3vector1d<as3double>    &monitordata,
                               CPoolMatrixAS3<as3double> &workarray) override;
    
    /*!
     * @brief Getter function which returns the number of working variables.
     *
     * @return mNVar.
     */
    unsigned short GetnVar(void) const override {return mNVar;}
  
  protected:
  
  private:
  	unsigned short mNVar = 4; ///< Number of working variables
};
```


### Example source file.

``solver_structure.cpp``

```C++
#include "solver_structure.hpp"


//-----------------------------------------------------------------------------------
// ISolver member functions.
//-----------------------------------------------------------------------------------


ISolver::ISolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
  :
    mZoneID(iZone)
 /*
  * Constructor for the interface solver class.
  */
{

}

//-----------------------------------------------------------------------------------

ISolver::~ISolver
(
 void
)
 /*
  * Destructor, which cleans up after the interface solver class. 
  */
{

}


//-----------------------------------------------------------------------------------
// CEESolver member functions.
//-----------------------------------------------------------------------------------


CEESolver::CEESolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
  :
    ISolver(config_container, geometry_container, iZone)
 /*
  * Constructor for the (non-linear) Euler equations class.
  */
{

}

//-----------------------------------------------------------------------------------

CEESolver::~CEESolver
(
 void
)
 /*
  * Destructor, which cleans up after the Euler equations class.
  */
{

}

//-----------------------------------------------------------------------------------

void CEESolver::InitPhysicalElements
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
  * Function that initializes the physical elements. 
  */
{

}

//-----------------------------------------------------------------------------------

void CEESolver::ComputeVolumeResidual
(
 size_t                     iElem,
 as3double                  localtime,
 as3vector1d<as3double>    &monitordata,
 CPoolMatrixAS3<as3double> &workarray
)
 /*
  * Function that computes the volume terms in a EE-type PDE. 
  */
{

}
```


