#include "stdafx.h"
#include "FEWarpSurfaceConstraint.h"
#include <FECore/FEElasticMaterial.h>

//=============================================================================
FEWarpSurface::FEWarpSurface()
{
}

//-----------------------------------------------------------------------------
void FEWarpSurface::Init()
{
	// create a new surface

	// allocate and initialize closest-point structure

}

//-----------------------------------------------------------------------------
// Find the closest-point projection of this point
vec3d FEWarpSurface::Project(vec3d r)
{
	vec3d q;
	vec2d s;
	m_cpp->Project(r, q, s);
	return q;
}

//=============================================================================
BEGIN_PARAMETER_LIST(FEWarpSurfaceConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_k      , FE_PARAM_DOUBLE  , "penalty" );
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL    , "laugon"  );
	ADD_PARAMETER(m_altol  , FE_PARAM_DOUBLE  , "altol"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEWarpSurfaceConstraint::FEWarpSurfaceConstraint(FEModel* pfem) : FEWarpConstraint(pfem)
{
	m_ptmp = 0;
	m_ptrg = 0;
}

//-----------------------------------------------------------------------------
void FEWarpSurfaceConstraint::Init()
{
	// read the two models
}

//-----------------------------------------------------------------------------
vec3d FEWarpSurfaceConstraint::wrpForce(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the initial and current position
	vec3d r0 = pt.r0;
	vec3d rt = pt.rt;

	vec3d qt, qs;

	vec3d qT = m_ptmp->Project(r0);
	vec3d qS = m_ptrg->Project(rt);

	double lT = (r0 - qT).norm();
	double lS = (rt - qS).norm();

	double T = 1.0/(1.0 + lT);
	double S = 1.0/(1.0 + lS);

	// get the target gradient
	// NOTE: this gradient must be perpendicular to the target surface at the closest point
	vec3d nS = rt - qS; nS.unit();

	return nS*(T - S);
}

//-----------------------------------------------------------------------------
mat3ds FEWarpSurfaceConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	return mat3ds();
}
