#pragma once

#include "option_structure.hpp"


struct CFlattenedFaceIndex
{
  size_t            mIndexFamily;
  size_t            mIndexFace;
  ETypeFaceGeometry mFaceType;
};


struct CFlattenedElementIndex
{
	unsigned short mIndexZone;
	size_t         mIndexElem;
};


