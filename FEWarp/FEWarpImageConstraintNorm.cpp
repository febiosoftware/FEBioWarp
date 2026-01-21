#include "stdafx.h"
#include "FEWarpImageConstraintNorm.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FEImgLib/image_tools.h>
#include <FECore/log.h>
#include <FECore/ad.h>

#include <FECore/FEOctreeSearch.h>
#include <FECore/FEModel.h>
#include <iostream>

using std::cout;

BEGIN_FECORE_CLASS(FEWarpImageConstraintNorm, FEWarpSingleImageConstraint);
	ADD_PARAMETER(m_windowSize, "window_size" );
END_FECORE_CLASS();

class TargetImage
{
public:
	TargetImage(ImageMap& smap) : m_smap(smap) {}
	
	ad::number operator () (const ad::vec3d& r) const
	{
		vec3d rr = r.values();
		assert(m_smap.valid(rr));

		double S = m_smap.value(rr);
		vec3d G = m_smap.gradient(rr);
		vec3d dr = r.partials();
		if (dr.x != 0) return ad::number(S, G.x);
		if (dr.y != 0) return ad::number(S, G.y);
		if (dr.z != 0) return ad::number(S, G.z);
		assert(false);
		return 0;
	}

	ad::vec3d grad(const ad::vec3d& r) const
	{
		vec3d rr = r.values();
		assert(m_smap.valid(rr));
		mat3ds H = m_smap.hessian(rr);
		vec3d G = m_smap.gradient(rr);
		vec3d dr = r.partials();
		if (dr.x != 0)
		{
			ad::number Gx(G.x, H.xx());
			ad::number Gy(G.y, H.xy());
			ad::number Gz(G.z, H.xz());
			return ad::vec3d(Gx, Gy, Gz);
		}
		else if (dr.y != 0)
		{
			ad::number Gx(G.x, H.xy());
			ad::number Gy(G.y, H.yy());
			ad::number Gz(G.z, H.yz());
			return ad::vec3d(Gx, Gy, Gz);
		}
		else if (dr.z != 0)
		{
			ad::number Gx(G.x, H.xz());
			ad::number Gy(G.y, H.yz());
			ad::number Gz(G.z, H.zz());
			return ad::vec3d(Gx, Gy, Gz);
		}
		assert(false);
		return ad::vec3d(0,0,0);
	}

private:
	ImageMap& m_smap;
};

//-----------------------------------------------------------------------------
FEWarpImageConstraintNorm::FEWarpImageConstraintNorm(FEModel* pfem) 
    : FEWarpSingleImageConstraint(pfem), m_tBar(0.0), m_tSigma(0.0), m_sBar(0.0), m_sSigma(0.0), m_init(false)
{
	m_tr0[0] = m_tr0[1] = m_tr0[2] = 0.0;
	m_tr1[0] = m_tr1[1] = m_tr1[2] = 1.0;
	m_sr0[0] = m_sr0[1] = m_sr0[2] = 0.0;
	m_sr1[0] = m_sr1[1] = m_sr1[2] = 1.0;

    
}

bool FEWarpImageConstraintNorm::Init()
{
    if (!FEWarpSingleImageConstraint::Init()) return false;

    Image& tImg = GetTemplateMap().GetImage();
    Image& sImg = GetTargetMap().GetImage();

    // Precompute global statistics
    int nVoxels = tImg.width() * tImg.height() * tImg.depth();
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

    std::cout << nVoxels << " voxels processed for global statistics." << std::endl;

    std:cout << "Template Image: mean = " << m_tBar << ", std. dev. = " << m_tSigma << std::endl;
    std::cout << "Source Image: mean = " << m_sBar << ", std. dev. = " << m_sSigma << std::endl;

    return true;
}

//-----------------------------------------------------------------------------
FEWarpImageConstraintNorm::~FEWarpImageConstraintNorm(void) {}

//-----------------------------------------------------------------------------
// vec3d FEWarpImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
// {
// 	// find the warping constraint
// 	// FEModel& fem = *GetFEModel();
// 	// FEWarpImageConstraintNorm* pc = 0;
// 	// for (int i=0; i<fem.NonlinearConstraints(); ++i)
// 	// {
// 	// 	pc = dynamic_cast<FEWarpImageConstraintNorm*>(fem.NonlinearConstraint(i));
// 	// 	if (pc) break;
// 	// }
// 	// if (pc == 0) return false;

//     FEModel& fem = *GetFEModel();

// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     double varianceThreshold = 1e-3; // tune this
//     double eps = 1e-12;

//     FEMesh& m = fem.GetMesh();

//     // FEOctreeSearch octree(&m);
//     // octree.Init();

// 	// int N = m.Nodes();
//     // int bad = 0;
//     // int mismatches = 0;
//     // #pragma omp parallel for
// 	// for (int i=0; i<N; ++i) 
// 	// {
//     vec3d r0 = mp.m_r0;
//  	vec3d rt = mp.m_rt;

//     // Build neighborhoods
//     int radius = 2;
    // int size = radius * 2 + 1;
    // double xScale = m_windowSize.x / size;
    // double yScale = m_windowSize.y / size;
    // double zScale = m_windowSize.z / size;
    // vector<vec3d> tNeighborhood;
    // vector<vec3d> sNeighborhood;
    // for(int dz = -radius; dz <= radius; ++dz)
    // {
    //     for(int dy = -radius; dy <= radius; ++dy)
    //     {
    //         for(int dx = -radius; dx <= radius; ++dx)
    //         {
    //             vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
    //             vec3d sTemp = rt + vec3d(dx * xScale, dy * yScale, dz * zScale);

    //             if(tmap.valid(tTemp) && smap.valid(sTemp))
    //             {
    //                 tNeighborhood.push_back(tTemp);
    //                 sNeighborhood.push_back(sTemp);
    //             }
    //         }
    //     }
    // }

//     // vector<vec3d> sNeighborhood;
//     // for(vec3d& tPos : tNeighborhood)
//     // {
//     //     double r[3];
//     //     FEElement* el = octree.FindElement(tPos, r);

//     //     if(!el)
//     //     {
//     //         continue;
//     //     }

//     //     vec3d v[FEElement::MAX_NODES];
//     //     for (int j = 0; j < el->Nodes(); ++j) v[j] = fem.GetMesh().Node(el->m_node[j]).m_rt;

//     //     vec3d sTemp = ((FESolidElement*)el)->evaluate(v, r[0], r[1], r[2]);  

//     //     if(smap.valid(sTemp))
//     //     {
//     //         sNeighborhood.push_back(sTemp);
//     //     }
//     // }

//     // Calc means
//     double tBar= 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         tBar += tmap.value(pos);
//     }
//     tBar /= tNeighborhood.size();

//     double sBar= 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         sBar += smap.value(pos);
//     }
//     sBar /= sNeighborhood.size();

//     // Calc other stats
//     double tSigma = 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         double diff = tmap.value(pos) - tBar;
//         tSigma += diff * diff;
//     }
//     tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

//     double sSigma = 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         double diff = smap.value(pos) - sBar;
//         sSigma += diff * diff;
//     }
//     sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

//     if (tSigma < varianceThreshold || sSigma < varianceThreshold)
//     {
//         // No image-based force in featureless regions
//         // Rely purely on elastic regularization
//         return vec3d(0, 0, 0);
//     }

//     double C_ts = 0.0;
//     for(int index = 0; index < tNeighborhood.size(); ++index)
//     {
//         double val =(tmap.value(tNeighborhood[index]) - tBar) * (smap.value(sNeighborhood[index]) - sBar);
//         C_ts += val;
//     }

//     if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
//     {
//         // std::cout << "Neighborhood size mismatch!" << std::endl;
//         // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//         // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//         // mismatches++;
//         return vec3d(0.0);
//     }

//     if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
//     {
//         // std::cout << "Small neighborhood size!" << std::endl;
//         // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//         // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//         // bad++;
//         return vec3d(0.0);
//     }

//     // Calc force
//     vec3d fw(0.0, 0.0, 0.0);
//     for(int index = 0; index < tNeighborhood.size(); ++index)
//     {
//         double tDiff = tmap.value(tNeighborhood[index]) - tBar;
//         double sDiff = smap.value(sNeighborhood[index]) - sBar;
//         vec3d G = smap.gradient(sNeighborhood[index]);

//         fw += G /(2*tSigma*sSigma + eps) * (tDiff - C_ts*sDiff / (sSigma * sSigma + eps));
//     }

//     // if(fw.norm() > 1)
//     // {
//     //     std::cout << "Large NCC force detected: " << fw.norm() << std::endl;
//     //     std::cout << "Node index: " << i << std::endl;
//     //     std::cout << "tSigma: " << tSigma << ", sSigma: " << sSigma << std::endl;
//     //     std::cout << "tBar: " << tBar << ", sBar: " << sBar << std::endl;
//     //     std::cout << "C_ts: " << C_ts << std::endl;
//     //     std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//     //     std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//     // }

//     return fw * (-m_k);
// 	// }

//     // std::cout << "Number of bad nodes: " << bad << " out of " << N << std::endl;
//     // std::cout << "Number of neighborhood size mismatches: " << mismatches << " out of " << N << std::endl;
// 	// return true;
// }

// vec3d FEWarpImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
// {
//     FEModel& fem = *GetFEModel();

// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     double varianceThreshold = 1e-6; // tune this
//     double eps = 1e-12;

//     FEMesh& m = fem.GetMesh();

//     vec3d r0 = mp.m_r0;
//  	vec3d rt = mp.m_rt;

//     // Build neighborhoods
//     int radius = 2;
//     int size = radius * 2 + 1;
//     double xScale = m_windowSize.x / size;
//     double yScale = m_windowSize.y / size;
//     double zScale = m_windowSize.z / size;
//     vector<vec3d> tNeighborhood;
//     vector<vec3d> sNeighborhood;
//     for(int dz = -radius; dz <= radius; ++dz)
//     {
//         for(int dy = -radius; dy <= radius; ++dy)
//         {
//             for(int dx = -radius; dx <= radius; ++dx)
//             {

//                 // vec3d tTemp = r0 + vec3d(dx * tmap.dx()*xScale, dy * tmap.dy()*yScale, dz * tmap.dz()*zScale);
//                 // vec3d sTemp = rt + vec3d(dx * smap.dx()*xScale, dy * smap.dy()*yScale, dz * smap.dz()*zScale);

//                 vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
//                 vec3d sTemp = rt + vec3d(dx * xScale, dy * yScale, dz * zScale);

//                 if(tmap.valid(tTemp) && smap.valid(sTemp))
//                 {
//                     tNeighborhood.push_back(tTemp);
//                     sNeighborhood.push_back(sTemp);
//                 }
//             }
//         }
//     }

//     // vector<vec3d> sNeighborhood;
//     // for(vec3d& tPos : tNeighborhood)
//     // {
//     //     double r[3];
//     //     FEElement* el = octree.FindElement(tPos, r);

//     //     if(!el)
//     //     {
//     //         continue;
//     //     }

//     //     vec3d v[FEElement::MAX_NODES];
//     //     for (int j = 0; j < el->Nodes(); ++j) v[j] = fem.GetMesh().Node(el->m_node[j]).m_rt;

//     //     vec3d sTemp = ((FESolidElement*)el)->evaluate(v, r[0], r[1], r[2]);  

//     //     if(smap.valid(sTemp))
//     //     {
//     //         sNeighborhood.push_back(sTemp);
//     //     }
//     // }

//     // Calc means
//     double tBar= 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         tBar += tmap.value(pos);
//     }
//     tBar /= tNeighborhood.size();

//     double sBar= 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         sBar += smap.value(pos);
//     }
//     sBar /= sNeighborhood.size();

//     // Calc other stats
//     double tSigma = 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         double diff = tmap.value(pos) - tBar;
//         tSigma += diff * diff;
//     }
//     tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

//     double sSigma = 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         double diff = smap.value(pos) - sBar;
//         sSigma += diff * diff;
//     }
//     sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

//     // if (tSigma < varianceThreshold || sSigma < varianceThreshold)
//     // {
//     //     // No image-based force in featureless regions
//     //     // Rely purely on elastic regularization
//     //     return vec3d(0, 0, 0);
//     // }

//     // if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
//     // {
//     //     // std::cout << "Neighborhood size mismatch!" << std::endl;
//     //     // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//     //     // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//     //     // mismatches++;
//     //     return vec3d(0.0);
//     // }

//     // if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
//     // {
//     //     // std::cout << "Small neighborhood size!" << std::endl;
//     //     // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//     //     // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//     //     // bad++;
//     //     return vec3d(0.0);
//     // }

//     double tHat = (tmap.value(r0) - tBar)/(tSigma + eps);
//     double sHat = (smap.value(rt) - sBar)/(sSigma + eps);

//     return smap.gradient(rt)*(-m_k/(sSigma + eps))*(tHat - sHat);
// }

// vec3d FEWarpImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
// {
// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     vec3d r0 = mp.m_r0;
// 	vec3d rt = mp.m_rt;

//     double eps = 1e-12;

//     // if(!tmap.valid(r0) || !smap.valid(rt)) return vec3d(0);
    

//     double tHat = (tmap.value(mp.m_r0) - m_tBar)/(m_tSigma + eps);
//     double sHat = (smap.value(mp.m_rt) - m_sBar)/(m_sSigma + eps);

//     vec3d Fw = smap.gradient(mp.m_rt)*(-m_k/(m_sSigma + eps))*(tHat - sHat);

//     // auto-differentation version
// 	// ad::vec3d rad = rt;
// 	// TargetImage Sad(m_smap);
// 	// auto E = [=](ad::vec3d& rt) { 
// 	// 	ad::number TS = tHat - (Sad(rt) - m_sBar)/(m_sSigma + eps);
// 	// 	return 0.5*m_k*TS*TS; 
// 	// 	};
// 	// vec3d Fad = ad::Grad(E, rt);
// 	// (end ad version)

// 	// double err = (Fad - Fw).norm();
// 	// assert(err < 1e-6);

//     return Fw;
// }

vec3d FEWarpImageConstraintNorm::wrpForce(FEMaterialPoint& mp)
{
    FEModel& fem = *GetFEModel();

	// get the image map
	ImageMap& tmap = GetTemplateMap();
	ImageMap& smap = GetTargetMap();

    double varianceThreshold = 1e-3; // tune this
    double eps = 1e-12;

    FEMesh& m = fem.GetMesh();

    vec3d r0 = mp.m_r0;
 	vec3d rt = mp.m_rt;

    if(!tmap.valid(r0) || !smap.valid(rt)) return vec3d(0);

    // Build neighborhoods
    int radius = 2;
    int size = radius * 2 + 1;
    double xScale = m_windowSize.x / size;
    double yScale = m_windowSize.y / size;
    double zScale = m_windowSize.z / size;
    vector<vec3d> tNeighborhood;
    vector<vec3d> sNeighborhood;
    for(int dz = -radius; dz <= radius; ++dz)
    {
        for(int dy = -radius; dy <= radius; ++dy)
        {
            for(int dx = -radius; dx <= radius; ++dx)
            {
                vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
                vec3d sTemp = vec3d(dx * xScale, dy * yScale, dz * zScale);

                if(tmap.valid(tTemp) && smap.valid(sTemp + rt))
                {
                    tNeighborhood.push_back(tTemp);
                    sNeighborhood.push_back(sTemp);
                }
            }
        }
    }

    // vector<vec3d> sNeighborhood;
    // for(vec3d& tPos : tNeighborhood)
    // {
    //     double r[3];
    //     FEElement* el = octree.FindElement(tPos, r);

    //     if(!el)
    //     {
    //         continue;
    //     }

    //     vec3d v[FEElement::MAX_NODES];
    //     for (int j = 0; j < el->Nodes(); ++j) v[j] = fem.GetMesh().Node(el->m_node[j]).m_rt;

    //     vec3d sTemp = ((FESolidElement*)el)->evaluate(v, r[0], r[1], r[2]);  

    //     if(smap.valid(sTemp))
    //     {
    //         sNeighborhood.push_back(sTemp);
    //     }
    // }

    if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
    {
        // std::cout << "Neighborhood size mismatch!" << std::endl;
        // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
        // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
        // mismatches++;
        return vec3d(0.0);
    }

    if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
    {
        // std::cout << "Small neighborhood size!" << std::endl;
        // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
        // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
        // bad++;
        return vec3d(0.0);
    }

    // Calc means
    double tBar= 0.0;
    for(vec3d& pos : tNeighborhood)
    {
        tBar += tmap.value(pos);
    }
    tBar /= tNeighborhood.size();

    double sBar= 0.0;
    for(vec3d& pos : sNeighborhood)
    {
        sBar += smap.value(pos + rt);
    }
    sBar /= sNeighborhood.size();

    // Calc other stats
    double tSigma = 0.0;
    for(vec3d& pos : tNeighborhood)
    {
        double diff = tmap.value(pos) - tBar;
        tSigma += diff * diff;
    }
    tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

    double sSigma = 0.0;
    for(vec3d& pos : sNeighborhood)
    {
        double diff = smap.value(pos + rt) - sBar;
        sSigma += diff * diff;
    }
    sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

    // if (tSigma < varianceThreshold || sSigma < varianceThreshold)
    // {
    //     // No image-based force in featureless regions
    //     // Rely purely on elastic regularization
    //     return vec3d(0, 0, 0);
    // }

    vec3d gradBar(0.0);
    for(vec3d& pos : sNeighborhood)
    {
        gradBar += smap.gradient(pos + rt);
    }
    gradBar /= sNeighborhood.size();

    vec3d dSigma = vec3d(0.0);
    for(vec3d& pos : sNeighborhood)
    {
        dSigma += (smap.gradient(pos + rt) - gradBar) * (smap.value(pos + rt) - sBar);
    }
    dSigma /= (sNeighborhood.size() * sSigma + eps);

    double tHat = (tmap.value(r0) - tBar)/(tSigma + eps);
    double sHat = (smap.value(rt) - sBar)/(sSigma + eps);

    vec3d Fw = (smap.gradient(rt) - gradBar - dSigma * (smap.value(rt) - sBar) / (sSigma + eps)) * (-m_k * (tHat - sHat) / (sSigma + eps));

    // auto-differentation version
	// TargetImage Sad(smap);
	// auto E = [&](ad::vec3d& rt) 
    // { 
    //     ad::number sBar= 0.0;
    //     for(vec3d& pos : sNeighborhood)
    //     {
    //         sBar = sBar + Sad(pos + rt);
    //     }
    //     sBar = sBar / sNeighborhood.size();

    //     ad::number sSigma = 0.0;
    //     for(vec3d& pos : sNeighborhood)
    //     {
    //         ad::number diff = Sad(pos + rt) - sBar;
    //         sSigma = sSigma + diff * diff;
    //     }
    //     sSigma = sqrt(sSigma / sNeighborhood.size() + eps);


	// 	ad::number TS = tHat - (Sad(rt) - sBar)/(sSigma + eps);
	// 	return 0.5*m_k*TS*TS; 
	
    
    // };
	// vec3d Fad = ad::Grad(E, rt);

	// double err = (Fad - Fw).norm();
	// assert(err < 1e-6);

    return Fw;
}

// mat3ds FEWarpImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
// {
//     FEModel& fem = *GetFEModel();

// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     double varianceThreshold = 1e-3; // tune this
//     double eps = 1e-12;

//     FEMesh& m = fem.GetMesh();

//     vec3d r0 = mp.m_r0;
//  	vec3d rt = mp.m_rt;

//     // Build neighborhoods
//     int radius = 2;
    // int size = radius * 2 + 1;
    // double xScale = m_windowSize.x / size;
    // double yScale = m_windowSize.y / size;
    // double zScale = m_windowSize.z / size;
    // vector<vec3d> tNeighborhood;
    // vector<vec3d> sNeighborhood;
    // for(int dz = -radius; dz <= radius; ++dz)
    // {
    //     for(int dy = -radius; dy <= radius; ++dy)
    //     {
    //         for(int dx = -radius; dx <= radius; ++dx)
    //         {
    //             vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
    //             vec3d sTemp = rt + vec3d(dx * xScale, dy * yScale, dz * zScale);

    //             if(tmap.valid(tTemp) && smap.valid(sTemp))
    //             {
    //                 tNeighborhood.push_back(tTemp);
    //                 sNeighborhood.push_back(sTemp);
    //             }
    //         }
    //     }
    // }

//     if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
//     {
//         // std::cout << "Neighborhood size mismatch!" << std::endl;
//         // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//         // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//         // mismatches++;
//         return mat3ds(0.0);
//     }

//     if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
//     {
//         // std::cout << "Small neighborhood size!" << std::endl;
//         // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//         // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//         // bad++;
//         return mat3ds(0.0);
//     }

// 	 // Calc means
//     double tBar= 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         tBar += tmap.value(pos);
//     }
//     tBar /= tNeighborhood.size();

//     double sBar= 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         sBar += smap.value(pos);
//     }
//     sBar /= sNeighborhood.size();

//     // Calc other stats
//     double tSigma = 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         double diff = tmap.value(pos) - tBar;
//         tSigma += diff * diff;
//     }
//     tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

//     double sSigma = 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         double diff = smap.value(pos) - sBar;
//         sSigma += diff * diff;
//     }
//     sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

//     if (tSigma < varianceThreshold || sSigma < varianceThreshold)
//     {
//         // No image-based force in featureless regions
//         // Rely purely on elastic regularization
//         return mat3ds(0.0);
//     }
	
//     mat3ds K(0.0);
//     for (int index = 0; index < tNeighborhood.size(); ++index)
//     {
//         vec3d gradS = smap.gradient(sNeighborhood[index]);
//         double magS = gradS.norm();

//         K += dyad(gradS) * (m_k / (2 * tSigma * sSigma + eps));
//     }

//     return K;
// }

// mat3ds FEWarpImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
// {
//     FEModel& fem = *GetFEModel();

// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     double varianceThreshold = 1e-6; // tune this
//     double eps = 1e-12;

//     FEMesh& m = fem.GetMesh();

//     vec3d r0 = mp.m_r0;
//  	vec3d rt = mp.m_rt;

//     // Build neighborhoods
//     int radius = 2;
//     int size = radius * 2 + 1;
//     double xScale = m_windowSize.x / size;
//     double yScale = m_windowSize.y / size;
//     double zScale = m_windowSize.z / size;
//     vector<vec3d> tNeighborhood;
//     vector<vec3d> sNeighborhood;
//     for(int dz = -radius; dz <= radius; ++dz)
//     {
//         for(int dy = -radius; dy <= radius; ++dy)
//         {
//             for(int dx = -radius; dx <= radius; ++dx)
//             {
//                 vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
//                 vec3d sTemp = rt + vec3d(dx * xScale, dy * yScale, dz * zScale);

//                 if(tmap.valid(tTemp) && smap.valid(sTemp))
//                 {
//                     tNeighborhood.push_back(tTemp);
//                     sNeighborhood.push_back(sTemp);
//                 }
//             }
//         }
//     }

//     // vector<vec3d> sNeighborhood;
//     // for(vec3d& tPos : tNeighborhood)
//     // {
//     //     double r[3];
//     //     FEElement* el = octree.FindElement(tPos, r);

//     //     if(!el)
//     //     {
//     //         continue;
//     //     }

//     //     vec3d v[FEElement::MAX_NODES];
//     //     for (int j = 0; j < el->Nodes(); ++j) v[j] = fem.GetMesh().Node(el->m_node[j]).m_rt;

//     //     vec3d sTemp = ((FESolidElement*)el)->evaluate(v, r[0], r[1], r[2]);  

//     //     if(smap.valid(sTemp))
//     //     {
//     //         sNeighborhood.push_back(sTemp);
//     //     }
//     // }

//     // Calc means
//     double tBar= 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         tBar += tmap.value(pos);
//     }
//     tBar /= tNeighborhood.size();

//     double sBar= 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         sBar += smap.value(pos);
//     }
//     sBar /= sNeighborhood.size();

//     // Calc other stats
//     double tSigma = 0.0;
//     for(vec3d& pos : tNeighborhood)
//     {
//         double diff = tmap.value(pos) - tBar;
//         tSigma += diff * diff;
//     }
//     tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

//     double sSigma = 0.0;
//     for(vec3d& pos : sNeighborhood)
//     {
//         double diff = smap.value(pos) - sBar;
//         sSigma += diff * diff;
//     }
//     sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

//     // if (tSigma < varianceThreshold || sSigma < varianceThreshold)
//     // {
//     //     // No image-based force in featureless regions
//     //     // Rely purely on elastic regularization
//     //     return mat3ds(0.0);
//     // }

//     // if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
//     // {
//     //     // std::cout << "Neighborhood size mismatch!" << std::endl;
//     //     // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//     //     // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//     //     // mismatches++;
//     //     return mat3ds(0.0);
//     // }

//     // if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
//     // {
//     //     // std::cout << "Small neighborhood size!" << std::endl;
//     //     // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
//     //     // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
//     //     // bad++;
//     //     return mat3ds(0.0);
//     // }

//     vec3d gradS = smap.gradient(rt);
//     mat3ds H = smap.hessian(rt);

//     double tHat = (tmap.value(r0) - tBar)/(tSigma + eps);
//     double sHat = (smap.value(rt) - sBar)/(sSigma + eps);

//     return  (dyad(gradS)*(m_k/(sSigma * sSigma + eps)) - H*(m_k/(sSigma + eps))*(tHat - sHat))*(-1);
// }

// mat3ds FEWarpImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
// {
// 	// get the image map
// 	ImageMap& tmap = GetTemplateMap();
// 	ImageMap& smap = GetTargetMap();

//     vec3d r0 = mp.m_r0;
//  	vec3d rt = mp.m_rt;

//     // if(!tmap.valid(r0) || !smap.valid(rt)) return mat3ds(0);

//     double eps = 1e-12;

//     vec3d gradS = smap.gradient(rt);
//     mat3ds H = smap.hessian(rt);

//     double T = tmap.value(r0);
//     double S = smap.value(rt);

//     double tHat = (tmap.value(r0) - m_tBar)/(m_tSigma + eps);
//     double sHat = (smap.value(rt) - m_sBar)/(m_sSigma + eps);

//     mat3ds K = (dyad(gradS)*(m_k/(m_sSigma * m_sSigma + eps)) - H*(m_k/(m_sSigma + eps))*(tHat - sHat))*(-1);


//     // TargetImage Sad(smap);

// 	// auto F = [=](ad::vec3d& rt) {
// 	// 	ad::number TS = tHat - (Sad(rt)- m_sBar)/(m_sSigma + eps) ;
// 	// 	return Sad.grad(rt)*(m_k/(m_sSigma + eps))*(TS);
// 	// 	};

// 	// mat3d Kad = ad::Grad(F, rt);
// 	// mat3d Ke = (K);

// 	// double err = (Kad - Ke).norm();
//     // std::cout << err << std::endl;
//     // assert(err < 1e-6);

//     // return  Kad.sym();
//     return K;
// }

mat3ds FEWarpImageConstraintNorm::wrpStiffness(FEMaterialPoint& mp)
{
    FEModel& fem = *GetFEModel();

	// get the image map
	ImageMap& tmap = GetTemplateMap();
	ImageMap& smap = GetTargetMap();

    double varianceThreshold = 1e-3; // tune this
    double eps = 1e-12;

    FEMesh& m = fem.GetMesh();

    vec3d r0 = mp.m_r0;
 	vec3d rt = mp.m_rt;

    if(!tmap.valid(r0) || !smap.valid(rt)) return mat3ds(0);

    // Build neighborhoods
    int radius = 2;
    int size = radius * 2 + 1;
    double xScale = m_windowSize.x / size;
    double yScale = m_windowSize.y / size;
    double zScale = m_windowSize.z / size;
    vector<vec3d> tNeighborhood;
    vector<vec3d> sNeighborhood;
    for(int dz = -radius; dz <= radius; ++dz)
    {
        for(int dy = -radius; dy <= radius; ++dy)
        {
            for(int dx = -radius; dx <= radius; ++dx)
            {
                vec3d tTemp = r0 + vec3d(dx * xScale, dy * yScale, dz * zScale);
                vec3d sTemp = vec3d(dx * xScale, dy * yScale, dz * zScale);

                if(tmap.valid(tTemp) && smap.valid(sTemp + rt))
                {
                    tNeighborhood.push_back(tTemp);
                    sNeighborhood.push_back(sTemp);
                }
            }
        }
    }

    // vector<vec3d> sNeighborhood;
    // for(vec3d& tPos : tNeighborhood)
    // {
    //     double r[3];
    //     FEElement* el = octree.FindElement(tPos, r);

    //     if(!el)
    //     {
    //         continue;
    //     }

    //     vec3d v[FEElement::MAX_NODES];
    //     for (int j = 0; j < el->Nodes(); ++j) v[j] = fem.GetMesh().Node(el->m_node[j]).m_rt;

    //     vec3d sTemp = ((FESolidElement*)el)->evaluate(v, r[0], r[1], r[2]);  

    //     if(smap.valid(sTemp))
    //     {
    //         sNeighborhood.push_back(sTemp);
    //     }
    // }

    if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
    {
        // std::cout << "Neighborhood size mismatch!" << std::endl;
        // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
        // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
        // mismatches++;
        return mat3ds(0.0);
    }

    if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
    {
        // std::cout << "Small neighborhood size!" << std::endl;
        // std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
        // std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
        // bad++;
        return mat3ds(0.0);
    }

    // Calc means
    double tBar= 0.0;
    for(vec3d& pos : tNeighborhood)
    {
        tBar += tmap.value(pos);
    }
    tBar /= tNeighborhood.size();

    // Calc other stats
    double tSigma = 0.0;
    for(vec3d& pos : tNeighborhood)
    {
        double diff = tmap.value(pos) - tBar;
        tSigma += diff * diff;
    }
    tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

    double tHat = (tmap.value(r0) - tBar)/(tSigma + eps);

    // auto-differentation version
	TargetImage Sad(smap);
	auto E = [&](ad::vec3d& rt) 
    { 
        ad::number sBar = 0.0;
        for(vec3d& pos : sNeighborhood)
        {
            sBar = sBar + Sad(pos + rt);
        }
        sBar = sBar / sNeighborhood.size();

        ad::number sSigma = 0.0;
        for(vec3d& pos : sNeighborhood)
        {
            ad::number diff = Sad(pos + rt) - sBar;
            sSigma = sSigma + diff * diff;
        }
        sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

        ad::vec3d gradBar(0,0,0);
        for(vec3d& pos : sNeighborhood)
        {
            gradBar = gradBar + Sad.grad(pos + rt);
        }
        gradBar = gradBar * (1.0 / sNeighborhood.size());

        ad::vec3d dSigma(0,0,0);
        for(vec3d& pos : sNeighborhood)
        {
            dSigma = dSigma + (Sad.grad(pos + rt) - gradBar) * (Sad(pos + rt) - sBar);
        }
        dSigma = dSigma * (1.0 / (sNeighborhood.size() * sSigma + eps));

        ad::number sHat = (Sad(rt) - sBar)/(sSigma + eps);

       return (Sad.grad(rt) - gradBar - dSigma * (Sad(rt) - sBar) *(1.0 / (sSigma + eps))) * (-m_k * (tHat - sHat) / (sSigma + eps));
    };
	mat3d K = ad::Grad(E, rt);

    return K.sym()*(-1);

    // Calc means
    // double tBar= 0.0;
    // for(vec3d& pos : tNeighborhood)
    // {
    //     tBar += tmap.value(pos);
    // }
    // tBar /= tNeighborhood.size();

    // double sBar= 0.0;
    // for(vec3d& pos : sNeighborhood)
    // {
    //     sBar += smap.value(pos + rt);
    // }
    // sBar /= sNeighborhood.size();

    // // Calc other stats
    // double tSigma = 0.0;
    // for(vec3d& pos : tNeighborhood)
    // {
    //     double diff = tmap.value(pos) - tBar;
    //     tSigma += diff * diff;
    // }
    // tSigma = sqrt(tSigma / tNeighborhood.size() + eps);

    // double sSigma = 0.0;
    // for(vec3d& pos : sNeighborhood)
    // {
    //     double diff = smap.value(pos + rt) - sBar;
    //     sSigma += diff * diff;
    // }
    // sSigma = sqrt(sSigma / sNeighborhood.size() + eps);

    // // if (tSigma < varianceThreshold || sSigma < varianceThreshold)
    // // {
    // //     // No image-based force in featureless regions
    // //     // Rely purely on elastic regularization
    // //     return vec3d(0, 0, 0);
    // // }

    // vec3d gradBar(0.0);
    // for(vec3d& pos : sNeighborhood)
    // {
    //     gradBar += smap.gradient(pos + rt);
    // }
    // gradBar /= sNeighborhood.size();

    // mat3ds H = smap.hessian(rt);
    // mat3ds hBar = smap.hessian(rt);
    // for(vec3d& pos : sNeighborhood)
    // {
    //     hBar += smap.hessian(pos + rt);
    // }
    // hBar /= sNeighborhood.size();

    // vec3d dSigma = vec3d(0.0);
    // for(vec3d& pos : sNeighborhood)
    // {
    //     dSigma += (smap.gradient(pos + rt) - gradBar) * (smap.value(pos + rt) - sBar);
    // }
    // dSigma /= (sNeighborhood.size() * sSigma + eps);

    // double tHat = (tmap.value(r0) - tBar)/(tSigma + eps);
    // double sHat = (smap.value(rt) - sBar)/(sSigma + eps);

    // vec3d V = (smap.gradient(rt) - gradBar - dSigma * (smap.value(rt) - sBar) / (sSigma + eps));

    // return (dyad(V)*m_k + (H - hBar)*(m_k * (tHat - sHat) / (sSigma + eps)))*(-1);
}