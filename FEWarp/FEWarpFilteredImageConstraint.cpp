#include "stdafx.h"
#include "FEWarpFilteredImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FEImgLib/ImageFilter.h>
#include <FEImgLib/Image.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <chrono>
#include <iostream>
#include <map>

//-----------------------------------------------------------------------------
FEWarpFilteredImageConstraint::FEWarpFilteredImageConstraint(FEModel* pfem) : FEWarpVolumeConstraint(pfem)
{
	m_pt = 0.0;

	m_tr0[0] = m_tr0[1] = m_tr0[2] = 0.0;
	m_tr1[0] = m_tr1[1] = m_tr1[2] = 1.0;
	m_sr0[0] = m_sr0[1] = m_sr0[2] = 0.0;
	m_sr1[0] = m_sr1[1] = m_sr1[2] = 1.0;

}

//-----------------------------------------------------------------------------
FEWarpFilteredImageConstraint::~FEWarpFilteredImageConstraint() {}

//-----------------------------------------------------------------------------
bool FEWarpFilteredImageConstraint::Init()
{
	if (FEWarpVolumeConstraint::Init() == false) return false;

	int nx = m_tmp0.width ();
	int ny = m_tmp0.height();
	int nz = m_tmp0.depth ();

	vec3d tr0(m_tr0[0], m_tr0[1], m_tr0[2]);
	vec3d tr1(m_tr1[0], m_tr1[1], m_tr1[2]);
	vec3d sr0(m_sr0[0], m_sr0[1], m_sr0[2]);
	vec3d sr1(m_sr1[0], m_sr1[1], m_sr1[2]);

	m_tmap.SetRange(tr0, tr1);
	m_smap.SetRange(sr0, sr1);

	m_tmp = m_tmp0;
	m_trg = m_trg0;

	m_tmap0.SetRange(tr0, tr1);
	m_smap0.SetRange(sr0, sr1);

	m_pt = 0.0;
	m_pr = 0.0;

	// SL: Added load curve update. Check with Steve for how to do just the single parameter so we don't overwrite anything.
	// update load curve
	FEModel* fem = GetFEModel();
	fem->EvaluateLoadParameters();
	
	// Initialize filters
	if (m_filt->Init() == false) return false;

	// blur images
	Update();

	return true;
}

//-----------------------------------------------------------------------------
// This is called at the beginning of each time step
// NOTE: This is currently called at the end of each iteration as well.
//       In future versions of FEBio this will no longer be the case.
//		 SL: Added check against current time and previous blur to skip blur if 
//		 the time or blur haven't changed
void FEWarpFilteredImageConstraint::Update()
{
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<double, std::ratio<1> >;
	std::chrono::time_point<Clock> m_beg;
	double elapsed;
	double ct = this->CurrentTime();
#ifndef NDEBUG
	std::cout << "current time is " << ct << std::endl;
#endif
	bool test1 = ct <= m_pt;
	bool test2 = ct == 0.0;
	bool test3 = m_filt->GetBlur() != m_pr;
#ifndef NDEBUG
	std::cout << "bools are " << test1 << ", " << test2 << ", " << test3 << std::endl;
#endif
	if ((ct <= m_pt) || (ct == 0.0)) { return; }
	m_pt = ct;
	
	if (m_filt->GetBlur() != m_pr)
	{
		//update i_pblur
		m_pr = m_filt->GetBlur();
		// template image
		m_beg = Clock::now();
		m_filt->Update(m_tmp, m_tmp0);
		elapsed = std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
		std::cout << "template image: took " << elapsed << " seconds" << std::endl;
		// target image
		m_beg = Clock::now();
		m_filt->Update(m_trg, m_trg0);
		elapsed = std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
		std::cout << "target image: took " << elapsed << " seconds" << std::endl;
	}
}

//=====================================================================
BEGIN_FECORE_CLASS(FEWarpSingleFilteredImageConstraint, FEWarpVolumeConstraint);
	ADD_PARAMETER(m_k      , "penalty" );
	ADD_PARAMETER(m_blaugon, "laugon"  );
	ADD_PARAMETER(m_altol  , "altol"   );
	ADD_PARAMETER(m_tr0    , 3, "template_range_min");
	ADD_PARAMETER(m_tr1    , 3, "template_range_max");
	ADD_PARAMETER(m_sr0, 3, "target_range_min");
	ADD_PARAMETER(m_sr1, 3, "target_range_max");

	ADD_PROPERTY(m_tmpReader, "template")->SetDefaultType("raw");
	ADD_PROPERTY(m_trgReader, "target"  )->SetDefaultType("raw");
	ADD_PROPERTY(m_filt, "filter");
END_FECORE_CLASS();

FEWarpSingleFilteredImageConstraint::FEWarpSingleFilteredImageConstraint(FEModel* fem) : FEWarpFilteredImageConstraint(fem)
{
	m_tmpReader = nullptr;
	m_trgReader = nullptr;
}

//-----------------------------------------------------------------------------
bool FEWarpSingleFilteredImageConstraint::Init()
{
	if ((m_tmpReader == nullptr) || (m_trgReader == nullptr)) return false;

	if (m_tmpReader->GetImage3D(m_tmp0) == false) return false;
	if (m_trgReader->GetImage3D(m_trg0) == false) return false;

	return FEWarpFilteredImageConstraint::Init();
}
