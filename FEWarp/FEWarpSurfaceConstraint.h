#include "FEWarpConstraint.h"
#include <FECore/FESurface.h>
#include <FECore/FEClosestPointProjection.h>

//-----------------------------------------------------------------------------
// This class implements a warping constraint on a surface definition
class FEWarpSurface
{
public:
	FEWarpSurface();
	~FEWarpSurface(){}

	// Initialization
	void Init();

	// project a point onto the surface
	vec3d Project(vec3d r);

protected:
	FEMesh	m_msh;	// the actual surface mesh

	FESurface*	m_psurf;	// surface definition of mesh

	FEClosestPointProjection*	m_cpp;	// closest point projector
};

//-----------------------------------------------------------------------------
// Warp constraint between two surfaces
class FEWarpSurfaceConstraint : public FEWarpConstraint
{
public:
	FEWarpSurfaceConstraint(FEModel* pfem);
	~FEWarpSurfaceConstraint() {}

	// Initialization
	void Init();

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt);

protected:
	FEWarpSurface*	m_ptmp;
	FEWarpSurface*	m_ptrg;

public:
	DECLARE_PARAMETER_LIST();
};
