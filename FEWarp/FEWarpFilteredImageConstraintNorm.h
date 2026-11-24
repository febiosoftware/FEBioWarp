#pragma once

#include "FEWarpFilteredImageConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>
#include <FEImgLib/ImageFilter.h> 
#include <map>

class FEWarpFilteredImageConstraintNorm : public FEWarpSingleFilteredImageConstraint
{
public:
    FEWarpFilteredImageConstraintNorm(FEModel* fem);
    ~FEWarpFilteredImageConstraintNorm() = default;

    void Update() override;

    vec3d wrpForce(FEMaterialPoint& mp) override;
    mat3ds wrpStiffness(FEMaterialPoint& mp) override;

    double GetTemplateMean() const { return m_tBar; }
    double GetTemplateSigma() const { return m_tSigma; }
    double GetTargetMean() const { return m_sBar; }
    double GetTargetSigma() const { return m_sSigma; }

private:
    bool Init() override;
    void ComputeStatistics();

private:
    double m_tBar = 0.0;
    double m_tSigma = 0.0;
    double m_sBar = 0.0;
    double m_sSigma = 0.0;

    DECLARE_FECORE_CLASS();
};