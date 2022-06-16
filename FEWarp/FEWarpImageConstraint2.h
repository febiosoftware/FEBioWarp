#pragma once
#include "FEWarpConstraint.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
class FEWarpImageConstraint2 : public FEWarpConstraint
{
public:
	FEWarpImageConstraint2(FEModel* pfem);
	~FEWarpImageConstraint2() {}

public: // inherited members
	void ShallowCopy(DumpStream& dmp, bool bsave) {};

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

	FEParamDouble	m_tmp;
	FEParamDouble	m_trg;

public:
	DECLARE_FECORE_CLASS();
};
