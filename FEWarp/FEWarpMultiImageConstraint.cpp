#include "stdafx.h"
#include "FEWarpMultiImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEWarpMultiImageConstraint, FEWarpImageConstraint);
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
FEWarpMultiImageConstraint::FEWarpMultiImageConstraint(FEModel* pfem) : FEWarpImageConstraint(pfem)
{
	m_tmpReader = nullptr;
	m_ntrg = -1;
}

//-----------------------------------------------------------------------------
bool FEWarpMultiImageConstraint::Init()
{
	if ((m_tmpReader == nullptr) || (m_trgReader.empty())) return false;

	// read the template image
	std::string tmpName = m_tmpReader->GetFileName();
	feLog("Reading template image %s\n", tmpName.c_str());
	if (m_tmpReader->GetImage3D(m_tmp0) == false) return false;

	// read the target image
	m_ntrg = 0;
	std::string trgName = m_trgReader[m_ntrg]->GetFileName();
	feLog("Reading target image %s\n", trgName.c_str());
	if (m_trgReader[m_ntrg]->GetImage3D(m_trg0) == false) return false;

	return FEWarpImageConstraint::Init();
}

//-----------------------------------------------------------------------------
// This is called at the beginning of each time step
// NOTE: This is currently called at the end of each iteration as well.
//       In future versions of FEBio this will no longer be the case.
void FEWarpMultiImageConstraint::Update()
{
	// see if we need to switch targets
	double t = GetTimeInfo().currentTime;
	double dt = 1.0 / (double) m_trgReader.size();

	int n = (int) floor(t / dt);
	if (n > m_trgReader.size()) n = m_trgReader.size();
	
	if (n != m_ntrg)
	{
		// switch targets
		m_ntrg = n;

		std::string trgName = m_trgReader[m_ntrg]->GetFileName();
		feLog("Reading target image %s\n", trgName.c_str());
		if (m_trgReader[m_ntrg]->GetImage3D(m_trg0) == false) throw std::exception("Can't read next target image");
		m_trg = m_trg0;

		// we need to reset the blur
		m_blur_cur = 0;
	}

	FEWarpImageConstraint::Update();
}
