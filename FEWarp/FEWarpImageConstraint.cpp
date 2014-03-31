#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEWarpImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_tmp0   , FE_PARAM_IMAGE_3D, "template");
	ADD_PARAMETER(m_trg0   , FE_PARAM_IMAGE_3D, "target"  );
	ADD_PARAMETER(m_k      , FE_PARAM_DOUBLE  , "penalty" );
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL    , "laugon"  );
	ADD_PARAMETER(m_altol  , FE_PARAM_DOUBLE  , "altol"   );
	ADD_PARAMETER(m_blur   , FE_PARAM_DOUBLE  , "blur"    );
	ADD_PARAMETERV(m_r0    , FE_PARAM_DOUBLEV, 3, "range_min");
	ADD_PARAMETERV(m_r1    , FE_PARAM_DOUBLEV, 3, "range_max");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void blur_image_2d(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int) d;
	float w = d - (float) n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg.Create(nx, ny, nz);
	Image tmp(src);
	float f[4];
	for (int l=0; l<n; ++l)
	{
		for (int k=0; k<nz; ++k)
			for (int j=0; j<ny; ++j)
				for (int i=0; i<nx; ++i)
				{
					if (i>0   ) f[0] = tmp.value(i-1,j,k); else f[0] = tmp.value(i, j, k);
					if (i<nx-1) f[1] = tmp.value(i+1,j,k); else f[1] = tmp.value(i, j, k);
					if (j>0   ) f[2] = tmp.value(i,j-1,k); else f[2] = tmp.value(i, j, k);
					if (j<ny-1) f[3] = tmp.value(i,j+1,k); else f[3] = tmp.value(i, j, k);
					trg.value(i,j,k) = 0.25f*(f[0] + f[1] + f[2] + f[3]);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
		for (int k=0; k<nz; ++k)
			for (int j=0; j<ny; ++j)
				for (int i=0; i<nx; ++i)
				{
					if (i>0   ) f[0] = tmp.value(i-1,j,k); else f[0] = tmp.value(i, j, k);
					if (i<nx-1) f[1] = tmp.value(i+1,j,k); else f[1] = tmp.value(i, j, k);
					if (j>0   ) f[2] = tmp.value(i,j-1,k); else f[2] = tmp.value(i, j, k);
					if (j<ny-1) f[3] = tmp.value(i,j+1,k); else f[3] = tmp.value(i, j, k);
					float f1 = 0.25f*(f[0] + f[1] + f[2] + f[3]);
					float f2 = trg.value(i,j,k);
					trg.value(i,j,k) = f1*w + f2*(1.f - w);
				}
	}
}

//-----------------------------------------------------------------------------
void blur_image(Image& trg, Image& src, float d)
{
	if (d <= 0) { trg = src; return; }

	int n = (int) d;
	float w = d - (float) n;

	int nx = src.width();
	int ny = src.height();
	int nz = src.depth();

	trg.Create(nx, ny, nz);
	Image tmp(src);
	float f[6];
	for (int l=0; l<n; ++l)
	{
		for (int k=0; k<nz; ++k)
			for (int j=0; j<ny; ++j)
				for (int i=0; i<nx; ++i)
				{
					if (i>0   ) f[0] = tmp.value(i-1,j,k); else f[0] = tmp.value(i, j, k);
					if (i<nx-1) f[1] = tmp.value(i+1,j,k); else f[1] = tmp.value(i, j, k);
					if (j>0   ) f[2] = tmp.value(i,j-1,k); else f[2] = tmp.value(i, j, k);
					if (j<ny-1) f[3] = tmp.value(i,j+1,k); else f[3] = tmp.value(i, j, k);
					if (k>0   ) f[4] = tmp.value(i,j,k-1); else f[4] = tmp.value(i, j, k);
					if (k<nz-1) f[5] = tmp.value(i,j,k+1); else f[5] = tmp.value(i, j, k);
					trg.value(i,j,k) = 0.1666667f*(f[0] + f[1] + f[2] + f[3] + f[4] + f[5]);
				}
		tmp = trg;
	}

	if (w > 0.0)
	{
		for (int k=0; k<nz; ++k)
			for (int j=0; j<ny; ++j)
				for (int i=0; i<nx; ++i)
				{
					if (i>0   ) f[0] = tmp.value(i-1,j,k); else f[0] = tmp.value(i, j, k);
					if (i<nx-1) f[1] = tmp.value(i+1,j,k); else f[1] = tmp.value(i, j, k);
					if (j>0   ) f[2] = tmp.value(i,j-1,k); else f[2] = tmp.value(i, j, k);
					if (j<ny-1) f[3] = tmp.value(i,j+1,k); else f[3] = tmp.value(i, j, k);
					if (k>0   ) f[4] = tmp.value(i,j,k-1); else f[4] = tmp.value(i, j, k);
					if (k<nz-1) f[5] = tmp.value(i,j,k+1); else f[5] = tmp.value(i, j, k);
					float f1 = 0.1666667f*(f[0] + f[1] + f[2] + f[3] + f[4] + f[5]);
					float f2 = trg.value(i,j,k);
					trg.value(i,j,k) = f1*w + f2*(1.f - w);
				}
	}
}

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{
	m_blur = 0.0;
	m_blur_cur = 0.0;
}

//-----------------------------------------------------------------------------
bool FEWarpImageConstraint::Init()
{
	FEWarpConstraint::Init();

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
