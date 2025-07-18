#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FEImgLib/ImageFilter.h>
#include <FEImgLib/Image.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{
	//m_blur = 0.0;
	m_blur_cur = 0.0;

	m_r0[0] = m_r0[1] = m_r0[2] = 0.0;
	m_r1[0] = m_r1[1] = m_r1[2] = 1.0;
}

//-----------------------------------------------------------------------------
FEWarpImageConstraint::~FEWarpImageConstraint() {}

//-----------------------------------------------------------------------------
bool FEWarpImageConstraint::Init()
{
	if (FEWarpConstraint::Init() == false) return false;

	int nx = m_tmp0.width ();
	int ny = m_tmp0.height();
	int nz = m_tmp0.depth ();

	vec3d r0(m_r0[0], m_r0[1], m_r0[2]);
	vec3d r1(m_r1[0], m_r1[1], m_r1[2]);

	m_tmap.SetRange(r0, r1);
	m_smap.SetRange(r0, r1);

	m_tmp = m_tmp0;
	m_trg = m_trg0;

	// SL: Added load curve update. Check with Steve for how to do just the single parameter so we don't overwrite anything.
	// update load curve
	FEModel* fem = GetFEModel();
	fem->EvaluateLoadParameters();
	
	// Initialize filters
	if (m_filt.size() < 1) { feLog("No filters provided.\n"); }
	for (auto iter = m_filt.begin(); iter != m_filt.end(); ++iter)
	{
		if ((*iter)->Init() == false) return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
// This is called at the beginning of each time step
// NOTE: This is currently called at the end of each iteration as well.
//       In future versions of FEBio this will no longer be the case.
void FEWarpImageConstraint::Update()
{
	for (auto iter = m_filt.begin(); iter != m_filt.end(); ++iter)
	{
		(*iter)->Update(m_trg, m_trg0);
	}
}

//-----------------------------------------------------------------------------
vec3d FEWarpImageConstraint::wrpForce(FEMaterialPoint& mp)
{
	// evaluate template
	double T = m_tmap.value(mp.m_r0);

	// evaluate target
	double S = m_smap.value(mp.m_rt);

	// evaluate target gradient
	vec3d G = m_smap.gradient(mp.m_rt);

	// evaluate force
	vec3d Fw = G*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpImageConstraint::wrpStiffness(FEMaterialPoint& mp)
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
	return H*((T - S)*m_k) - dyad(dS)*m_k;
}

//=====================================================================
BEGIN_FECORE_CLASS(FEWarpSingleImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
	ADD_PARAMETER(m_r0    , 3, "range_min");
	ADD_PARAMETER(m_r1    , 3, "range_max");

	ADD_PROPERTY(m_tmpReader, "template")->SetDefaultType("raw");
	ADD_PROPERTY(m_trgReader, "target"  )->SetDefaultType("raw");
	ADD_PROPERTY(m_filt, "filter");
END_FECORE_CLASS();

FEWarpSingleImageConstraint::FEWarpSingleImageConstraint(FEModel* fem) : FEWarpImageConstraint(fem)
{
	m_tmpReader = nullptr;
	m_trgReader = nullptr;
}

//-----------------------------------------------------------------------------
bool FEWarpSingleImageConstraint::Init()
{
	if ((m_tmpReader == nullptr) || (m_trgReader == nullptr)) return false;

	if (m_tmpReader->GetImage3D(m_tmp0) == false) return false;
	if (m_trgReader->GetImage3D(m_trg0) == false) return false;

	return FEWarpImageConstraint::Init();
}
