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
double FEWarpSurface::Project(vec3d r)
{
	return 0;
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
}

//-----------------------------------------------------------------------------
vec3d FEWarpSurfaceConstraint::wrpForce(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the initial and current position
	vec3d r0 = pt.r0;
	vec3d rt = pt.rt;

	vec3d qt, qs;

	double T = m_ptmp->Project(r0);
	double S = m_ptrg->Project(rt);

	// get the target gradient
	// NOTE: this gradient must be perpendicular to the target surface at the closest point
	vec3d nS;

	return nS*(T - S);
}

//-----------------------------------------------------------------------------
mat3ds FEWarpSurfaceConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	return mat3ds();
}
