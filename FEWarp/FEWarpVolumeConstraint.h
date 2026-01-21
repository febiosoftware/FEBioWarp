#pragma once
#include "FEWarpConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>

//-----------------------------------------------------------------------------
class FEWarpVolumeConstraint : public FEWarpConstraint
{
public:
	FEWarpVolumeConstraint(FEModel* pfem);
	~FEWarpVolumeConstraint(void);

public: // inherited members
	void ShallowCopy(DumpStream& dmp, bool bsave) {};

public:
	ImageMap& GetTemplateMap() { return m_tmap; }
	ImageMap& GetTargetMap  () { return m_smap; }
	ImageMap& GetRawTemplateMap() { return m_tmap0; }
	ImageMap& GetRawTargetMap() { return m_smap0; }

public:
	//! Calculate the force at a material point
	virtual vec3d wrpForce(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	virtual mat3ds wrpStiffness(FEMaterialPoint& pt);

protected:
	double	m_tr0[3];	//!< minimum range
	double	m_tr1[3];	//!< maximum range
	double	m_sr0[3];	//!< minimum range
	double	m_sr1[3];	//!< maximum range
	

	Image	m_tmp;	//!< template image
	Image	m_trg;	//!< target image

	ImageMap m_tmap;
	ImageMap m_smap;

	Image m_tmp0;	//!< original template image data (pre-blurred)
	Image m_trg0;	//!< original target image data (pre-blurred)

	ImageMap m_tmap0;
	ImageMap m_smap0;
};
