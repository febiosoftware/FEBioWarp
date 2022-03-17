#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "image_tools.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEWarpImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
	ADD_PARAMETER(m_blur   , "blur"    );
	ADD_PARAMETER(m_r0    , 3, "range_min");
	ADD_PARAMETER(m_r1    , 3, "range_max");

	ADD_PROPERTY(m_tmpReader, "template")->SetDefaultType("raw");
	ADD_PROPERTY(m_trgReader, "target"  )->SetDefaultType("raw");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{
	m_tmpReader = nullptr;
	m_trgReader = nullptr;

	m_blur = 0.0;
	m_blur_cur = 0.0;
}

//-----------------------------------------------------------------------------
bool FEWarpImageConstraint::Init()
{
	if ((m_tmpReader == nullptr) || (m_trgReader == nullptr)) return false;

	if (FEWarpConstraint::Init() == false) return false;

	if (m_tmpReader->GetImage3D(m_tmp0) == false) return false;
	if (m_trgReader->GetImage3D(m_trg0) == false) return false;

	int nx = m_tmp0.width ();
	int ny = m_tmp0.height();
	int nz = m_tmp0.depth ();

	vec3d r0(m_r0[0], m_r0[1], m_r0[2]);
	vec3d r1(m_r1[0], m_r1[1], m_r1[2]);

	m_tmap.SetRange(r0, r1);
	m_smap.SetRange(r0, r1);

	m_tmp = m_tmp0;
	m_trg = m_trg0;

	m_blur_cur = m_blur;
	if (m_blur > 0)
	{
		if (m_tmp0.depth() == 1) blur_image_2d(m_tmp, m_tmp0, (float) m_blur); else blur_image(m_tmp, m_tmp0, (float) m_blur);
		if (m_trg0.depth() == 1) blur_image_2d(m_trg, m_trg0, (float) m_blur); else blur_image(m_trg, m_trg0, (float) m_blur);
	}

	return true;
}

//-----------------------------------------------------------------------------
// This is called at the beginning of each time step
// NOTE: This is currently called at the end of each iteration as well.
//       In future versions of FEBio this will no longer be the case.
void FEWarpImageConstraint::Update()
{
	if (m_blur !=  m_blur_cur)
	{
		m_blur_cur = m_blur;
		if (m_tmp0.depth() == 1) blur_image_2d(m_tmp, m_tmp0, (float) m_blur); else blur_image(m_tmp, m_tmp0, (float) m_blur);
		if (m_trg0.depth() == 1) blur_image_2d(m_trg, m_trg0, (float) m_blur); else blur_image(m_trg, m_trg0, (float) m_blur);
	}
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
