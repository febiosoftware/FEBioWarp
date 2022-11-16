#include "stdafx.h"
#include "FEWarpSurfaceConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FEModel.h>

//=============================================================================
FEWarpSurface::FEWarpSurface(FEModel* fem) : FESurface(fem)
{
	m_binit = true;
	m_beta = 50.0;
}

//-----------------------------------------------------------------------------
// Find the closest-point projection of this point
vec3d FEWarpSurface::Project(vec3d x)
{
	vec3d q;
	vec2d s;
//	ClosestReferencePointProjection(x, q, s, m_binit, 0.001);
	m_binit = false;
	return q;
}

//-----------------------------------------------------------------------------
double FEWarpSurface::value(vec3d x)
{
	vec3d q = Project(x);
	double l = (q - x).norm();
	return 1.0/(1.0 + m_beta*l);
}

//-----------------------------------------------------------------------------
void FEWarpSurface::Update()
{
	m_binit = true;
}

//=============================================================================
BEGIN_FECORE_CLASS(FEWarpSurfaceConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEWarpSurfaceConstraint::FEWarpSurfaceConstraint(FEModel* pfem) : FEWarpConstraint(pfem)
{
	m_ptmp = new FEWarpSurface(pfem);
	m_ptrg = new FEWarpSurface(pfem);
}

//-----------------------------------------------------------------------------
FESurface* FEWarpSurfaceConstraint::GetSurface(const char* sz)
{
	if (strcmp(sz,"template") == 0) return m_ptmp;
	if (strcmp(sz, "target" ) == 0) return m_ptrg;
	return 0;
}

//-----------------------------------------------------------------------------
bool FEWarpSurfaceConstraint::Init()
{
	FEWarpConstraint::Init();

	m_ptmp->Init();
	m_ptrg->Init();

	return true;
}

//-----------------------------------------------------------------------------
vec3d FEWarpSurfaceConstraint::wrpForce(FEMaterialPoint& mp)
{
	// get the initial and current position
	vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

	vec3d qt, qs;

	vec3d qT = m_ptmp->Project(r0);
	vec3d qS = m_ptrg->Project(rt);

	double lT = (r0 - qT).norm();
	double lS = (rt - qS).norm();

	double T = 1.0/(1.0 + m_ptmp->m_beta*lT);
	double S = 1.0/(1.0 + m_ptmp->m_beta*lS);

	// get the target gradient
	// NOTE: this gradient must be perpendicular to the target surface at the closest point
	vec3d nS = rt - qS; nS.unit();

	return nS*(T - S)*m_k;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpSurfaceConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the initial and current position
	vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

	vec3d qt, qs;

	vec3d qT = m_ptmp->Project(r0);
	vec3d qS = m_ptrg->Project(rt);

	double lT = (r0 - qT).norm();
	double lS = (rt - qS).norm();

	double T = 1.0/(1.0 + m_ptmp->m_beta*lT);
	double S = 1.0/(1.0 + m_ptrg->m_beta*lS);

	// get the target gradient
	// NOTE: this gradient must be perpendicular to the target surface at the closest point
	vec3d nS = rt - qS; nS.unit();

//	return H*((T - S)*m_k) - dyad(dS)*m_k;
	return -dyad(nS)*m_k;
}

//-----------------------------------------------------------------------------
void FEWarpSurfaceConstraint::Update()
{
	m_ptmp->Update();
	m_ptrg->Update();
}
