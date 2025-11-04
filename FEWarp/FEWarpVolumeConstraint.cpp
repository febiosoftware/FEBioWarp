#include "stdafx.h"
#include "FEWarpVolumeConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FEWarpVolumeConstraint::FEWarpVolumeConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg), m_tmap0(m_tmp0), m_smap0(m_trg0)
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

	if (m_tmap.valid(r0) && m_smap.valid(rt))
	{
		// evaluate template
		double T = m_tmap.value(mp.m_r0);

		// evaluate target
		double S = m_smap.value(mp.m_rt);

		// evaluate target gradient
		vec3d G = m_smap.gradient(mp.m_rt);

		// evaluate force
		vec3d Fw = G * ((S - T) * m_k);

		return Fw;
	}
	else
		return vec3d(0.0);
}

//-----------------------------------------------------------------------------
mat3ds FEWarpVolumeConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

	if (m_tmap.valid(r0) && m_smap.valid(rt))
	{
		// template value
		double T = m_tmap.value(mp.m_r0);

		// target value
		double S = m_smap.value(mp.m_rt);

		// calculate target gradient
		vec3d dS = m_smap.gradient(mp.m_rt);

		// calculate target hessian
		mat3ds H = m_smap.hessian(mp.m_rt);

		// warping stiffness
		return H * ((T - S) * m_k) - dyad(dS) * m_k;
	}
	else
		return mat3ds(0.0);
}