#pragma once
#include "FEWarpImageConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>
// #include <FECore/FEOctreeSearch.h>

//-----------------------------------------------------------------------------
class FEWarpImageConstraintNorm : public FEWarpSingleImageConstraint
{
public:
	FEWarpImageConstraintNorm(FEModel* pfem);
	~FEWarpImageConstraintNorm(void);

    void Update() override;

public: // inherited members
	void ShallowCopy(DumpStream& dmp, bool bsave) {};

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt) override;

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt) override;

    double getTBar() const { return m_tBar; }
    double getTSigma() const { return m_tSigma; }
    double getSBar() const { return m_sBar; }
    double getSSigma() const { return m_sSigma; }

private:
    bool Init();

    void ComputeStatistics();

protected:
    vec3d   m_windowSize; //!< neighborhood size for normalization

    double m_tBar;   //!< template image mean
    double m_tSigma; //!< template image std. dev.
    double m_sBar;   //!< source image mean
    double m_sSigma; //!< source image std. dev.

    // FEOctreeSearch m_octree;

    DECLARE_FECORE_CLASS();
};
