#include "stdafx.h"
#include "FEWarpImageConstraint2.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEWarpImageConstraint2, FEWarpConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
	ADD_PARAMETER(m_blur   , "blur"    );
	ADD_PARAMETER(m_r0    , 3, "range_min");
	ADD_PARAMETER(m_r1    , 3, "range_max");

	ADD_PARAMETER(m_tmp, "template");
	ADD_PARAMETER(m_trg, "target");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEWarpImageConstraint2::FEWarpImageConstraint2(FEModel* pfem) : FEWarpConstraint(pfem)
{
	m_tmp = 0.0;
	m_trg = 0.0;

	m_blur = 0.0;
	m_blur_cur = 0.0;
}

//-----------------------------------------------------------------------------
vec3d FEWarpImageConstraint2::wrpForce(FEMaterialPoint& mp)
{
	// evaluate template
	double T = m_tmp(mp);

	// evaluate target
	double S = m_trg(mp);

	// evaluate target gradient
	vec3d G;// = m_trg.gradient(pt.m_rt);

	// evaluate force
	vec3d Fw = G*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpImageConstraint2::wrpStiffness(FEMaterialPoint& mp)
{
	// template value
	double T = m_tmp(mp);

	// target value
	double S = m_trg(mp);

	// calculate target gradient
	vec3d dS;// = m_trg.gradient(pt.m_rt);

	// calculate target hessian
	mat3ds H;// = m_trg.hessian(pt.m_rt);

	// warping stiffness
	return H*((T - S)*m_k) - dyad(dS)*m_k;
}
