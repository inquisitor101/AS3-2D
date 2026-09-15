#include "face_element.hpp"


//-----------------------------------------------------------------------------------
// IElementFace member functions.
//-----------------------------------------------------------------------------------

IElementFace::IElementFace
(
 size_t iElemM,
 size_t iElemP,
 size_t nVar,
 size_t nSol2D
)
 /*
  *
  */
{

}

//-----------------------------------------------------------------------------------

IElementFace::~IElementFace
(
 void
)
 /*
  *
  */
{

}


//-----------------------------------------------------------------------------------
// CInternalElementFace member functions.
//-----------------------------------------------------------------------------------

CInternalElementFace::CInternalElementFace
(
 size_t iElemM,
 size_t iElemP,
 size_t nVar,
 size_t nSol2D
)
  :
    IElementFace(iElemM, iElemP, nVar, nSol2D)
 /*
  *
  */
{
  mElementIndexM = iElemM;
  mElementIndexP = iElemP;

  mResMinus.resize(nVar, nSol2D);
}

//-----------------------------------------------------------------------------------

CInternalElementFace::~CInternalElementFace
(
 void
)
 /*
  *
  */
{

}


//-----------------------------------------------------------------------------------
// CInterfaceElementFace member functions.
//-----------------------------------------------------------------------------------

CInterfaceElementFace::CInterfaceElementFace
(
 size_t index_interface
)
  :
    IElementFace(index_interface)
 /*
  *
  */
{
  mIndexInterface = index_interface;
}

//-----------------------------------------------------------------------------------

CInterfaceElementFace::~CInterfaceElementFace
(
 void
)
 /*
  *
  */
{

}

//-----------------------------------------------------------------------------------




