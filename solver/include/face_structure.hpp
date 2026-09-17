#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometry;


struct CInternalFaceGeometry
{
	size_t mIndexElementM;
	size_t mIndexElementP;
	
	unsigned short mIndexZone;
};

struct CBoundaryFaceGeometry
{
	size_t mIndexElement;
	size_t mIndexZone;

	EFaceLocation mFaceLocation;
};

struct CInterfaceFaceGeometry
{
	size_t mIndexElementM;
	size_t mIndexElementP;
	
	unsigned short mIndexZoneM;
	unsigned short mIndexZoneP;

	EFaceLocation mFaceLocationM;
	EFaceLocation mFaceLocationP;
};



class CMultizoneFaceGeometry
{
	public:


		void InitializeInternalFaces(const CGeometry *geometry_container);

		void InitializeInterfaceFaces(const CConfig   *config_container,
				                          const CGeometry *geometry_structure);


		ETypeFaceGeometry GetFaceTypeFromIndex(size_t i) const
		{
			if( i < GetnInternalFaces() )                       return ETypeFaceGeometry::INTERNAL;
			if( i < GetnInternalFaces() + GetnBoundaryFaces() ) return ETypeFaceGeometry::BOUNDARY;
			if( i < GetnFacesTotal() )                          return ETypeFaceGeometry::INTERFACE;

			ERROR("Invalid face index specified.");
		}

		CInternalFaceGeometry& GetInternalFace(size_t i)
		{
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInternalFaces[i]; 
		}
		const CInternalFaceGeometry& GetInternalFace(size_t i) const 
		{
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInternalFaces[i]; 
		}

		CBoundaryFaceGeometry& GetBoundaryFace(size_t i)
		{ 
#if DEBUG
			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mBoundaryFaces[i]; 
		}
		const CBoundaryFaceGeometry& GetBoundaryFace(size_t i) const 
		{ 
#if DEBUG
			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mBoundaryFaces[i]; 
		}

		CInterfaceFaceGeometry& GetInterfaceFace(size_t i)
		{
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInterfaceFaces[i]; 
		}
		const CInterfaceFaceGeometry& GetInterfaceFace(size_t i) const 
		{
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInterfaceFaces[i]; 
		}

		size_t GetnInternalFaces(void)  const { return mInternalFaces.size();  }
		size_t GetnBoundaryFaces(void)  const { return mBoundaryFaces.size();  }
		size_t GetnInterfaceFaces(void) const { return mInterfaceFaces.size(); }

		size_t GetnFacesTotal(void) const
		{
			return GetnInternalFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

	private:
		as3vector1d<CInternalFaceGeometry>  mInternalFaces;
		as3vector1d<CBoundaryFaceGeometry>  mBoundaryFaces;
		as3vector1d<CInterfaceFaceGeometry> mInterfaceFaces;
};





