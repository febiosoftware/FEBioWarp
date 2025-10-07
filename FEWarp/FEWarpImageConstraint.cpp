#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpVolumeConstraint(pfem)
{
	m_blur = 0.0;
	m_blur_cur = 0.0;
	m_blur_method = 0; // default average blur

	m_tr0[0] = m_tr0[1] = m_tr0[2] = 0.0;
	m_tr1[0] = m_tr1[1] = m_tr1[2] = 1.0;
	m_sr0[0] = m_sr0[1] = m_sr0[2] = 0.0;
	m_sr1[0] = m_sr1[1] = m_sr1[2] = 1.0;
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

	vec3d tr0(m_tr0[0], m_tr0[1], m_tr0[2]);
	vec3d tr1(m_tr1[0], m_tr1[1], m_tr1[2]);
	// SL: Leave this unchanged for now by copying the template range.
	m_sr0[0] = m_tr0[0]; m_sr0[1] = m_tr0[1]; m_sr0[2] = m_tr0[2];
	m_sr1[0] = m_tr1[0]; m_sr1[1] = m_tr1[1]; m_sr1[2] = m_tr1[2];
	vec3d sr0(m_sr0[0], m_sr0[1], m_sr0[2]);
	vec3d sr1(m_sr1[0], m_sr1[1], m_sr1[2]);

	m_tmap.SetRange(tr0, tr1);
	m_smap.SetRange(sr0, sr1);

	m_tmp = m_tmp0;
	m_trg = m_trg0;

	BlurMethod blurMethod = (BlurMethod)m_blur_method;

	m_blur_cur = m_blur;
	if (m_blur > 0)
	{
		feLog("Blurring images, blur factor %lg\n", m_blur);

		if (m_tmp0.depth() == 1) blur_image_2d(m_tmp, m_tmp0, (float)m_blur, blurMethod); else blur_image_3d(m_tmp, m_tmp0, (float)m_blur, blurMethod);
		if (m_trg0.depth() == 1) blur_image_2d(m_trg, m_trg0, (float)m_blur, blurMethod); else blur_image_3d(m_trg, m_trg0, (float)m_blur, blurMethod);
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

		feLog("Blurring images, blur factor %lg\n", m_blur);

		BlurMethod blurMethod = (BlurMethod)m_blur_method;

		if (m_tmp0.depth() == 1) blur_image_2d(m_tmp, m_tmp0, (float)m_blur, blurMethod); else blur_image_3d(m_tmp, m_tmp0, (float)m_blur, blurMethod);
		if (m_trg0.depth() == 1) blur_image_2d(m_trg, m_trg0, (float)m_blur, blurMethod); else blur_image_3d(m_trg, m_trg0, (float)m_blur, blurMethod);
	}
}

//=====================================================================
BEGIN_FECORE_CLASS(FEWarpSingleImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
	ADD_PARAMETER(m_blur   , "blur"    );
	ADD_PARAMETER(m_blur_method   , "blur_method"    )->setEnums("AVERAGE\0FFT\0");
	//SL: Leave parameter name unchanged for now.
	ADD_PARAMETER(m_tr0    , 3, "range_min");
	ADD_PARAMETER(m_tr1    , 3, "range_max");
	//SL: Don't expose for now to leave this unchanged.
	//ADD_PARAMETER(m_sr0, 3, "target_range_min");
	//ADD_PARAMETER(m_sr1, 3, "target_range_max");

	ADD_PROPERTY(m_tmpReader, "template")->SetDefaultType("raw");
	ADD_PROPERTY(m_trgReader, "target"  )->SetDefaultType("raw");
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
