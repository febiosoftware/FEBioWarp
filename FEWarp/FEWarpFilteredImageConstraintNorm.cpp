#include "FEWarpFilteredImageConstraint.h"
#include "FEWarpFilteredImageConstraintNorm.h"
#include <FECore/log.h>
#include <cmath>

BEGIN_FECORE_CLASS(FEWarpFilteredImageConstraintNorm, FEWarpSingleFilteredImageConstraint);
END_FECORE_CLASS();

FEWarpFilteredImageConstraintNorm::FEWarpFilteredImageConstraintNorm(FEModel* fem) : FEWarpSingleFilteredImageConstraint(fem)
{
}

bool FEWarpFilteredImageConstraintNorm::Init()
{
    if (FEWarpSingleFilteredImageConstraint::Init() == false)
    {
        return false;
    }

    ComputeStatistics();
    return true;
}

void FEWarpFilteredImageConstraintNorm::Update()
{
    const double previousBlur = m_pr;

    FEWarpFilteredImageConstraint::Update();

    if (m_pr != previousBlur)
    {
        ComputeStatistics();
    }
}

void FEWarpFilteredImageConstraintNorm::ComputeStatistics()
{
    Image& tImg = GetTemplateMap().GetImage();
    Image& sImg = GetTargetMap().GetImage();

    auto calculateStatistics = [](Image& image, double& mean, double& sigma)
        {
            mean = 0.0;
            sigma = 0.0;

            const int nx = image.width();
            const int ny = image.height();
            const int nz = image.depth();
            const int n = nx * ny * nz;

            for (int k = 0; k < nz; ++k)
            {
                for (int j = 0; j < ny; ++j)
                {
                    for (int i = 0; i < nx; ++i)
                    {
                        mean += image.value(i, j, k);
                    }
                }
            }

            mean /= static_cast<double>(n);

            for (int k = 0; k < nz; ++k)
            {
                for (int j = 0; j < ny; ++j)
                {
                    for (int i = 0; i < nx; ++i)
                    {
                        const double d = image.value(i, j, k) - mean;
                        sigma += d * d;
                    }
                }
            }

            sigma = std::sqrt(sigma / static_cast<double>(n));
        };

    calculateStatistics(tImg, m_tBar, m_tSigma);
    calculateStatistics(sImg, m_sBar, m_sSigma);

    feLog(
        "Filtered global image statistics: "
        "T=(%lg, %lg), S=(%lg, %lg)\n",
        m_tBar, m_tSigma, m_sBar, m_sSigma);
}

vec3d FEWarpFilteredImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
{
    ImageMap& tmap = GetTemplateMap();
    ImageMap& smap = GetTargetMap();

    const vec3d r0 = mp.m_r0;
    const vec3d rt = mp.m_rt;

    if (!tmap.valid(r0) || !smap.valid(rt))
    {
        return vec3d(0.0);
    }

    constexpr double eps = 1e-12;

    const double tHat =
        (tmap.value(r0) - m_tBar) / (m_tSigma + eps);

    const double sHat =
        (smap.value(rt) - m_sBar) / (m_sSigma + eps);

    return smap.gradient(rt)
        * (-m_k / (m_sSigma + eps))
        * (tHat - sHat);
}

mat3ds FEWarpFilteredImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
{
    // get the image map
    ImageMap& tmap = GetTemplateMap();
    ImageMap& smap = GetTargetMap();

    vec3d r0 = mp.m_r0;
    vec3d rt = mp.m_rt;

    if (!tmap.valid(r0) || !smap.valid(rt)) return mat3ds(0);

    double eps = 1e-12;

    vec3d gradS = smap.gradient(rt);
    mat3ds H = smap.hessian(rt);

    double T = tmap.value(r0);
    double S = smap.value(rt);

    double tHat = (tmap.value(r0) - m_tBar) / (m_tSigma + eps);
    double sHat = (smap.value(rt) - m_sBar) / (m_sSigma + eps);
    mat3ds K = (dyad(gradS) * (m_k / (m_sSigma * m_sSigma)) - H * (m_k / (m_sSigma + eps)) * (tHat - sHat)) * (-1);

    return K;
}