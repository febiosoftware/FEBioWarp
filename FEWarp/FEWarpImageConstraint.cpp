#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEWarpImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_tmp    , FE_PARAM_IMAGE_3D, "template");
	ADD_PARAMETER(m_trg    , FE_PARAM_IMAGE_3D, "target"  );
	ADD_PARAMETER(m_k      , FE_PARAM_DOUBLE  , "penalty" );
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL    , "laugon"  );
	ADD_PARAMETER(m_altol  , FE_PARAM_DOUBLE  , "altol"   );
	ADD_PARAMETERV(m_r0    , FE_PARAM_DOUBLEV, 3, "range_min");
	ADD_PARAMETERV(m_r1    , FE_PARAM_DOUBLEV, 3, "range_max");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{

}

//-----------------------------------------------------------------------------
bool FEWarpImageConstraint::Init()
{
	FEWarpConstraint::Init();

	int nx = m_tmp.width ();
	int ny = m_tmp.height();
	int nz = m_tmp.depth ();

	vec3d r0(m_r0[0], m_r0[1], m_r0[2]);
	vec3d r1(m_r1[0], m_r1[1], m_r1[2]);

	m_tmap.SetRange(r0, r1);
	m_smap.SetRange(r0, r1);

	return true;
}

//-----------------------------------------------------------------------------
vec3d FEWarpImageConstraint::wrpForce(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// evaluate template
	double T = m_tmap.value(pt.m_r0);

	// evaluate target
	double S = m_smap.value(pt.m_rt);

	// evaluate target gradient
	vec3d G = m_smap.gradient(pt.m_rt);

	// evaluate force
	vec3d Fw = G*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpImageConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// template value
	double T = m_tmap.value(pt.m_r0);

	// target value
	double S = m_smap.value(pt.m_rt);

	// calculate target gradient
	vec3d dS = m_smap.gradient(pt.m_rt);

	// calculate target hessian
	mat3ds H = m_smap.hessian(pt.m_rt);

	// warping stiffness
	return H*((T - S)*m_k) - dyad(dS)*m_k;
}
