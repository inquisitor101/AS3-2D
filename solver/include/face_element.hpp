#pragma once

#include "option_structure.hpp"


class IElementFace
{
  public:

    IElementFace(size_t index_interface) {} // for testing now
    IElementFace(size_t iElemM, 
                 size_t iElemP,
                 size_t nVar,
                 size_t nSol2D);

    virtual ~IElementFace(void);

    virtual ETypeElementFace GetTypeElementFace(void) const = 0;


    CMatrixAS3<as3double>& GetResMinus(void) { return mResMinus; } 
    

  protected:
    CMatrixAS3<as3double> mResMinus; // for the minus (left/bottom) elements.
  
  private:

};

//-----------------------------------------------------------------------------------

class CInternalElementFace : public IElementFace
{
  public:

    CInternalElementFace(size_t iElemM,
                         size_t iElemP,
                         size_t nVar,
                         size_t nSol2D);
    
    ~CInternalElementFace(void) override;

    ETypeElementFace GetTypeElementFace(void) const override { return ETypeElementFace::INTERNAL; }

    size_t GetElementIndexM(void) const {return mElementIndexM;}
    size_t GetElementIndexP(void) const {return mElementIndexP;}
  protected:

  private:
    size_t mElementIndexM; // minus (left or bottom)
    size_t mElementIndexP; // plus  (right or top)
};

//-----------------------------------------------------------------------------------

class CBoundaryElementFace : public IElementFace
{
  public:

    ETypeElementFace GetTypeElementFace(void) const override { return ETypeElementFace::BOUNDARY; }

  protected:

  private:

};

//-----------------------------------------------------------------------------------

class CInterfaceElementFace : public IElementFace
{
  public:

    CInterfaceElementFace(size_t index_interface);

    ~CInterfaceElementFace(void) override;

    ETypeElementFace GetTypeElementFace(void) const override { return ETypeElementFace::INTERFACE; }

  protected:

  private:
    size_t mIndexInterface;
};




