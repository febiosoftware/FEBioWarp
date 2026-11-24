#include "stdafx.h"
#include "FEWarpImageConstraintNorm.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>
#include <FECore/ad.h>

// #include <FECore/FEOctreeSearch.h>
#include <FECore/FEModel.h>
#include <iostream>

using std::cout;

BEGIN_FECORE_CLASS(FEWarpImageConstraintNorm, FEWarpSingleImageConstraint);
	ADD_PARAMETER(m_windowSize, "window_size" );
END_FECORE_CLASS();

//SL: Used for autodiff features. Commented out for now.
//class TargetImage
//{
//public:
//	TargetImage(ImageMap& smap) : m_smap(smap) {}
//	
//	ad::number operator () (const ad::vec3d& r) const
//	{
//		vec3d rr = r.values();
//		assert(m_smap.valid(rr));
//
//		double S = m_smap.value(rr);
//		vec3d G = m_smap.gradient(rr);
//		vec3d dr = r.partials();
//		if (dr.x != 0) return ad::number(S, G.x);
//		if (dr.y != 0) return ad::number(S, G.y);
//		if (dr.z != 0) return ad::number(S, G.z);
//		assert(false);
//		return 0;
//	}
//
//	ad::vec3d grad(const ad::vec3d& r) const
//	{
//		vec3d rr = r.values();
//		assert(m_smap.valid(rr));
//		mat3ds H = m_smap.hessian(rr);
//		vec3d G = m_smap.gradient(rr);
//		vec3d dr = r.partials();
//		if (dr.x != 0)
//		{
//			ad::number Gx(G.x, H.xx());
//			ad::number Gy(G.y, H.xy());
//			ad::number Gz(G.z, H.xz());
//			return ad::vec3d(Gx, Gy, Gz);
//		}
//		else if (dr.y != 0)
//		{
//			ad::number Gx(G.x, H.xy());
//			ad::number Gy(G.y, H.yy());
//			ad::number Gz(G.z, H.yz());
//			return ad::vec3d(Gx, Gy, Gz);
//		}
//		else if (dr.z != 0)
//		{
//			ad::number Gx(G.x, H.xz());
//			ad::number Gy(G.y, H.yz());
//			ad::number Gz(G.z, H.zz());
//			return ad::vec3d(Gx, Gy, Gz);
//		}
//		assert(false);
//		return ad::vec3d(0,0,0);
//	}
//
//private:
//	ImageMap& m_smap;
//};

//-----------------------------------------------------------------------------
FEWarpImageConstraintNorm::FEWarpImageConstraintNorm(FEModel* pfem) 
    : FEWarpSingleImageConstraint(pfem), m_tBar(0.0), m_tSigma(0.0), m_sBar(0.0), m_sSigma(0.0)//, m_octree(&pfem->GetMesh())
{
	m_tr0[0] = m_tr0[1] = m_tr0[2] = 0.0;
	m_tr1[0] = m_tr1[1] = m_tr1[2] = 1.0;
	m_sr0[0] = m_sr0[1] = m_sr0[2] = 0.0;
	m_sr1[0] = m_sr1[1] = m_sr1[2] = 1.0;
}

//-----------------------------------------------------------------------------
FEWarpImageConstraintNorm::~FEWarpImageConstraintNorm(void) {}

void FEWarpImageConstraintNorm::Update()
{
    //bool updateStats = m_blur_cur != m_blur;

    FEWarpImageConstraint::Update();

    /*if(updateStats)
    {*/
    ComputeStatistics();
    //}
}

bool FEWarpImageConstraintNorm::Init()
{
    if (!FEWarpSingleImageConstraint::Init()) return false;

    ComputeStatistics();

    // m_octree.Init();

    return true;
}

void FEWarpImageConstraintNorm::ComputeStatistics()
{
    Image& tImg = GetTemplateMap().GetImage();
    Image& sImg = GetTargetMap().GetImage();

    m_tBar = 0.0;
    m_tSigma = 0.0;
    m_sBar = 0.0;
    m_sSigma = 0.0;

    // Precompute global statistics
    int nVoxels = tImg.width() * tImg.height() * tImg.depth();
    if (nVoxels <= 0)
    {
        std::cerr << "Error: Template image has no voxels." << std::endl;
        return;
	}
    for(int k = 0; k < tImg.depth(); ++k)
    {
        for(int j = 0; j < tImg.height(); ++j)
        {
            for(int i = 0; i < tImg.width(); ++i)
            {
                m_tBar += tImg.value(i, j, k);
                m_sBar += sImg.value(i, j, k);
            }
        }
    }
    m_tBar /= nVoxels;
    m_sBar /= nVoxels;

    for(int k = 0; k < tImg.depth(); ++k)
    {
        for(int j = 0; j < tImg.height(); ++j)
        {
            for(int i = 0; i < tImg.width(); ++i)
            {
                double tDiff = tImg.value(i, j, k) - m_tBar;
                double sDiff = sImg.value(i, j, k) - m_sBar;

                m_tSigma += tDiff * tDiff;
                m_sSigma += sDiff * sDiff;
            }
        }
    }
    
    m_tSigma = sqrt(m_tSigma / nVoxels);
    m_sSigma = sqrt(m_sSigma / nVoxels);

#ifndef NDEBUG
    std::cout << nVoxels << " voxels processed for global statistics." << std::endl;
    std::cout << "Template Image: mean = " << m_tBar << ", std. dev. = " << m_tSigma << std::endl;
    std::cout << "Source Image: mean = " << m_sBar << ", std. dev. = " << m_sSigma << std::endl;
#endif // !NDEBUG
}

vec3d FEWarpImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
{
	// get the image map
	ImageMap& tmap = GetTemplateMap();
	ImageMap& smap = GetTargetMap();

    vec3d r0 = mp.m_r0;
	vec3d rt = mp.m_rt;

    constexpr double eps = 1e-12;

    if(!tmap.valid(r0) || !smap.valid(rt)) return vec3d(0);
    

    double tHat = (tmap.value(mp.m_r0) - m_tBar)/(m_tSigma + eps);
    //double tHat = (tmap.value(mp.m_r0) - m_tBar) / (m_tSigma);
    double sHat = (smap.value(mp.m_rt) - m_sBar)/(m_sSigma + eps);
    //double sHat = (smap.value(mp.m_rt) - m_sBar) / (m_sSigma);

    //vec3d Fw = smap.gradient(mp.m_rt)*(-m_k/(m_sSigma + eps))*(tHat - sHat);
    vec3d Fw = smap.gradient(mp.m_rt) * (-m_k / (m_sSigma)) * (tHat - sHat);

    return Fw;
}

mat3ds FEWarpImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
{
	// get the image map
	ImageMap& tmap = GetTemplateMap();
	ImageMap& smap = GetTargetMap();

    vec3d r0 = mp.m_r0;
 	vec3d rt = mp.m_rt;

     if(!tmap.valid(r0) || !smap.valid(rt)) return mat3ds(0);

    constexpr double eps = 1e-12;

    vec3d gradS = smap.gradient(rt);
    mat3ds H = smap.hessian(rt);

    double T = tmap.value(r0);
    double S = smap.value(rt);

    double tHat = (tmap.value(r0) - m_tBar) / (m_tSigma + eps);
    double sHat = (smap.value(rt) - m_sBar) / (m_sSigma + eps);
    mat3ds K = (dyad(gradS) * (m_k / (m_sSigma * m_sSigma)) - H * (m_k / (m_sSigma + eps)) * (tHat - sHat)) * (-1);
    
    return K;
}
