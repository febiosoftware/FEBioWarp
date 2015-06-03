#pragma once
#include "FEWarpConstraint.h"
#include <FECore/FESurface.h>
#include <FECore/FEClosestPointProjection.h>

//-----------------------------------------------------------------------------
// This class implements a warping constraint on a surface definition
class FEWarpSurface : public FESurface
{
public:
	FEWarpSurface(FEMesh* pm);
	~FEWarpSurface(){}

public:
	vec3d Project(vec3d x);
	double value(vec3d x);

	void Update();

public:
	bool	m_binit;	// has NQ been initialized
	double	m_beta;
};

//-----------------------------------------------------------------------------
// Warp constraint between two surfaces
class FEWarpSurfaceConstraint : public FEWarpConstraint
{
public:
	FEWarpSurfaceConstraint(FEModel* pfem);
	~FEWarpSurfaceConstraint() {}

	// Initialization
	bool Init();

	// update
	void Update();

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt);

public:
	FESurface* GetSurface(const char* sz);

	FEWarpSurface* GetTemplate() { return m_ptmp; }
	FEWarpSurface* GetTarget  () { return m_ptrg; }

	void ShallowCopy(DumpStream& dmp, bool bsave) {};

protected:
	FEWarpSurface*	m_ptmp;
	FEWarpSurface*	m_ptrg;

public:
	DECLARE_PARAMETER_LIST();
};
