#pragma once
#include "FEWarpConstraint.h"

//-----------------------------------------------------------------------------
class FEWarpImageConstraint : public FEWarpConstraint
{
public:
	FEWarpImageConstraint(FEModel* pfem);
	~FEWarpImageConstraint() {}

	// initialization
	void Init();

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

	Image	m_tmp;	//!< template image
	Image	m_trg;	//!< target image

	ImageMap m_tmap;
	ImageMap m_smap;

public:
	DECLARE_PARAMETER_LIST();
};
