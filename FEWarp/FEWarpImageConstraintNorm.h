#pragma once
#include "FEWarpImageConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>

//-----------------------------------------------------------------------------
class FEWarpImageConstraintNorm : public FEWarpSingleImageConstraint
{
public:
	FEWarpImageConstraintNorm(FEModel* pfem);
	~FEWarpImageConstraintNorm(void);

public: // inherited members
	void ShallowCopy(DumpStream& dmp, bool bsave) {};

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt) override;

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt) override;

private:
    bool Init();

protected:
    vec3d   m_windowSize; //!< neighborhood size for normalization

    bool m_init;

    double m_tBar;   //!< template image mean
    double m_tSigma; //!< template image std. dev.
    double m_sBar;   //!< source image mean
    double m_sSigma; //!< source image std. dev.

    DECLARE_FECORE_CLASS();
};
