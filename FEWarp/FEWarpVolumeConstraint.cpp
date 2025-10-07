#include "stdafx.h"
#include "FEWarpVolumeConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FEWarpVolumeConstraint::FEWarpVolumeConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{
	m_tr0[0] = m_tr0[1] = m_tr0[2] = 0.0;
	m_tr1[0] = m_tr1[1] = m_tr1[2] = 1.0;
	m_sr0[0] = m_sr0[1] = m_sr0[2] = 0.0;
	m_sr1[0] = m_sr1[1] = m_sr1[2] = 1.0;
}

//-----------------------------------------------------------------------------
FEWarpVolumeConstraint::~FEWarpVolumeConstraint(void) {}

//-----------------------------------------------------------------------------
vec3d FEWarpVolumeConstraint::wrpForce(FEMaterialPoint& mp)
{
	vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

	// evaluate template
	double T = m_tmap.valid(r0) ? m_tmap.value(mp.m_r0) : 0.0;

	// evaluate target
	double S = m_smap.valid(rt) ? m_smap.value(mp.m_rt) : 0.0;

	// evaluate target gradient
	vec3d G = m_smap.valid(rt) ? m_smap.gradient(mp.m_rt) : vec3d(0.0);

	// evaluate force
	vec3d Fw = G*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpVolumeConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

	// template value
	double T = m_tmap.valid(r0) ? m_tmap.value(mp.m_r0) : 0.0;

	// target value
	double S = m_smap.valid(rt) ? m_smap.value(mp.m_rt) : 0.0;

	// calculate target gradient
	vec3d dS = m_smap.valid(rt) ? m_smap.gradient(mp.m_rt) : vec3d(0.0);

	// calculate target hessian
	mat3ds H = m_smap.valid(rt) ? m_smap.hessian(mp.m_rt) : mat3ds(0.0);

	// warping stiffness
	return H*((T - S)*m_k) - dyad(dS)*m_k;
}