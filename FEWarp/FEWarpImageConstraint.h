#pragma once
#include "FEWarpConstraint.h"

//-----------------------------------------------------------------------------
class FEWarpImageConstraint : public FEWarpConstraint
{
public:
	FEWarpImageConstraint(FEModel* pfem);
	~FEWarpImageConstraint() {}

	// initialization
	bool Init();

	// update
	void Update();

public: // inherited members
	void ShallowCopy(DumpStream& dmp, bool bsave) {};

public:
	ImageMap& GetTemplateMap() { return m_tmap; }
	ImageMap& GetTargetMap  () { return m_smap; }

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt);

protected:
	double	m_r0[3];	//!< minimum range
	double	m_r1[3];	//!< maximum range
	double	m_blur;		//!< blur factor
	
	double		m_blur_cur;	//!< current blur factor

	Image	m_tmp;	//!< template image
	Image	m_trg;	//!< target image

	ImageMap m_tmap;
	ImageMap m_smap;
	
	Image m_tmp0;	//!< original template image data (pre-blurred)
	Image m_trg0;	//!< original target image data (pre-blurred)

public:
	DECLARE_FECORE_CLASS();
};
