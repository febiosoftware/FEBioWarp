#include "stdafx.h"
#include "FEWarpPlot.h"
#include "FEWarpConstraint.h"
#include "FEWarpVolumeConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FEDataStream.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEOctreeSearch.h>

#include <iostream>

#define SCALAR 1

bool FEPlotTemplate::Save(FEMesh &m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);

		FEWarpVolumeConstraint* pci = dynamic_cast<FEWarpVolumeConstraint*>(pc);
		if (pci && pci->IsActive()) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs && pcs->IsActive()) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}

bool FEPlotTemplate::SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint* pc, FEDataStream& s)
{
	// get the template image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		if (tmap.valid(r0))
			s << tmap.value(r0);
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotTemplate::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s)
{
	FEWarpSurface* ps = pc->GetTemplate();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		s << ps->value(m.Node(i).m_r0);
	}
	return true;
}

bool FEPlotTarget::Save(FEMesh &m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);

		FEWarpVolumeConstraint* pci = dynamic_cast<FEWarpVolumeConstraint*>(pc);
		if (pci && pci->IsActive()) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs && pcs->IsActive()) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}


bool FEPlotTarget::SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint* pc, FEDataStream& s)
{
	// get the target image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d rt = m.Node(i).m_rt;
		if (smap.valid(rt))
			s << smap.value(rt);
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotTarget::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s)
{
	FEWarpSurface* ps = pc->GetTarget();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
		s << ps->value(m.Node(i).m_rt);
	return true;
}

bool FEPlotEnergy::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);
			s << (0.5 * (T - S) * (T - S));
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotNCCEnergy::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

    double varianceThreshold = 1e-3; // tune this
    double eps = 1e-12;

	int N = m.Nodes();
    int bad = 0;
    #pragma omp parallel for
	for (int i=0; i<N; ++i) 
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

        // Build neighborhoods
        int radius = 2;
        vector<vec3d> tNeighborhood;
        vector<vec3d> sNeighborhood;
        for(int dz = -radius; dz <= radius; ++dz)
        {
            for(int dy = -radius; dy <= radius; ++dy)
            {
                for(int dx = -radius; dx <= radius; ++dx)
                {
                    vec3d tTemp = r0 + vec3d(dx * tmap.dx()*SCALAR, dy * tmap.dy()*SCALAR, dz * tmap.dz()*SCALAR);
                    vec3d sTemp = rt + vec3d(dx * smap.dx()*SCALAR, dy * smap.dy()*SCALAR, dz * smap.dz()*SCALAR);

                    if(tmap.valid(tTemp) && smap.valid(sTemp))
                    {
                        tNeighborhood.push_back(tTemp);
                        sNeighborhood.push_back(sTemp);
                    }
                }
            }
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
            sBar += smap.value(pos);
        }
        sBar /= sNeighborhood.size();

        // Calc other stats
        double tSigma = 0.0;
        for(vec3d& pos : tNeighborhood)
        {
            double diff = tmap.value(pos) - tBar;
            tSigma += diff * diff;
        }
        tSigma = sqrt(tSigma + eps);

        double sSigma = 0.0;
        for(vec3d& pos : sNeighborhood)
        {
            double diff = smap.value(pos) - sBar;
            sSigma += diff * diff;
        }
        sSigma = sqrt(sSigma + eps);

        if (tSigma < varianceThreshold || sSigma < varianceThreshold) 
        {
            s << 0;
            continue;
        }

        double C_ts = 0.0;
        for(int index = 0; index < tNeighborhood.size(); ++index)
        {
            C_ts += (tmap.value(tNeighborhood[index]) - tBar) * (smap.value(sNeighborhood[index]) - sBar);
        }

        s << 0.5 *(1 -  C_ts / (tSigma * sSigma + eps));
	}

    std::cout << "Number of bad nodes: " << bad << " out of " << N << std::endl;
	return true;
}

bool FEPlotForce::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	double lam = ppc->GetPenalty();
	int N = m.Nodes();
	for (int i=0; i<N; ++i) 
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);
			vec3d G = smap.gradient(rt);
			vec3d fw = - G * ((T - S)) * lam;

			s << fw;
		}
		else 
			s << vec3d(0.0);
	}
	return true;
}

bool FEPlotNCCForce::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

    double varianceThreshold = 1e-3; // tune this
    double eps = 1e-12;

    FEOctreeSearch octree(&fem.GetMesh());
    octree.Init();

	int N = m.Nodes();
    int bad = 0;
    int mismatches = 0;
    #pragma omp parallel for
	for (int i=0; i<N; ++i) 
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

        // Build neighborhoods
        int radius = 2;
        vector<vec3d> tNeighborhood;
        vector<vec3d> sNeighborhood;
        for(int dz = -radius; dz <= radius; ++dz)
        {
            for(int dy = -radius; dy <= radius; ++dy)
            {
                for(int dx = -radius; dx <= radius; ++dx)
                {
                    vec3d tTemp = r0 + vec3d(dx * tmap.dx()*SCALAR, dy * tmap.dy()*SCALAR, dz * tmap.dz()*SCALAR);
                    vec3d sTemp = rt + vec3d(dx * smap.dx()*SCALAR, dy * smap.dy()*SCALAR, dz * smap.dz()*SCALAR);

                    if(tmap.valid(tTemp) && smap.valid(sTemp))
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
            sBar += smap.value(pos);
        }
        sBar /= sNeighborhood.size();

        // Calc other stats
        double tSigma = 0.0;
        for(vec3d& pos : tNeighborhood)
        {
            double diff = tmap.value(pos) - tBar;
            tSigma += diff * diff;
        }
        tSigma = sqrt(tSigma + eps);

        double sSigma = 0.0;
        for(vec3d& pos : sNeighborhood)
        {
            double diff = smap.value(pos) - sBar;
            sSigma += diff * diff;
        }
        sSigma = sqrt(sSigma + eps);

        if (tSigma < varianceThreshold || sSigma < varianceThreshold)
        {
            // No image-based force in featureless regions
            // Rely purely on elastic regularization
            s << vec3d(0, 0, 0);
            continue;
        }

        double C_ts = 0.0;
        for(int index = 0; index < tNeighborhood.size(); ++index)
        {
            C_ts += (tmap.value(tNeighborhood[index]) - tBar) * (smap.value(sNeighborhood[index]) - sBar);
        }

        if(tNeighborhood.size() != sNeighborhood.size() || tNeighborhood.size() == 0)
        {
            std::cout << "Neighborhood size mismatch!" << std::endl;
            std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
            std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
            mismatches++;
            s << vec3d(0.0);
            continue;
        }

        if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
        {
            // std::cout << "Small neighborhood size!" << std::endl;
            std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
            std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
            bad++;
            s << vec3d(0.0);
            continue;
        }

        // Calc force
        vec3d fw(0.0, 0.0, 0.0);
        for(int index = 0; index < tNeighborhood.size(); ++index)
        {
            double tDiff = tmap.value(tNeighborhood[index]) - tBar;
            double sDiff = smap.value(sNeighborhood[index]) - sBar;
            vec3d G = smap.gradient(sNeighborhood[index]);

            fw += G /(2*tSigma*sSigma + eps) * (tDiff - C_ts*sDiff / (sSigma * sSigma + eps));
        }

        // if(fw.norm() > 1)
        // {
        //     std::cout << "Large NCC force detected: " << fw.norm() << std::endl;
        //     std::cout << "Node index: " << i << std::endl;
        //     std::cout << "tSigma: " << tSigma << ", sSigma: " << sSigma << std::endl;
        //     std::cout << "tBar: " << tBar << ", sBar: " << sBar << std::endl;
        //     std::cout << "C_ts: " << C_ts << std::endl;
        //     std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
        //     std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
        // }

        s << fw;
	}

    std::cout << "Number of bad nodes: " << bad << " out of " << N << std::endl;
    std::cout << "Number of neighborhood size mismatches: " << mismatches << " out of " << N << std::endl;
	return true;
}

bool FEPlotNCCForce2::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

    double varianceThreshold = 1e-3; // tune this
    double eps = 1e-12;

    FEOctreeSearch octree(&fem.GetMesh());
    octree.Init();

	int N = m.Nodes();
    int bad = 0;
    int mismatches = 0;
    #pragma omp parallel for
	for (int i=0; i<N; ++i) 
	{

		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

        // Build neighborhoods
        int radius = 2;
        vector<vec3d> tNeighborhood;
        vector<vec3d> sNeighborhood;
        for(int dz = -radius; dz <= radius; ++dz)
        {
            for(int dy = -radius; dy <= radius; ++dy)
            {
                for(int dx = -radius; dx <= radius; ++dx)
                {
                    vec3d tTemp = r0 + vec3d(dx * tmap.dx()*SCALAR, dy * tmap.dy()*SCALAR, dz * tmap.dz()*SCALAR);
                    vec3d sTemp = rt + vec3d(dx * smap.dx()*SCALAR, dy * smap.dy()*SCALAR, dz * smap.dz()*SCALAR);

                    if(tmap.valid(tTemp) && smap.valid(sTemp))
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
            std::cout << "Neighborhood size mismatch!" << std::endl;
            std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
            std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
            mismatches++;
            s << vec3d(0.0);
            continue;
        }

        if(tNeighborhood.size() < 50 || sNeighborhood.size() < 50)
        {
            // std::cout << "Small neighborhood size!" << std::endl;
            std::cout << "tNeighborhood size: " << tNeighborhood.size() << std::endl;
            std::cout << "sNeighborhood size: " << sNeighborhood.size() << std::endl;
            bad++;
            s << vec3d(0.0);
            continue;
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
            sBar += smap.value(pos);
        }
        sBar /= sNeighborhood.size();

        vec3d gradSBar(0.0, 0.0, 0.0);
        for(vec3d& pos : sNeighborhood)
        {
            gradSBar += smap.gradient(pos);
        }
        gradSBar /= sNeighborhood.size();

        double tSigma = 0.0;
        for(vec3d& pos : tNeighborhood)
        {
            double diff = tmap.value(pos) - tBar;
            tSigma += diff * diff;
        }
        tSigma = sqrt(tSigma + eps);

        double sSigma = 0.0;
        for(vec3d& pos : sNeighborhood)
        {
            double diff = smap.value(pos) - sBar;
            sSigma += diff * diff;
        }
        sSigma = sqrt(sSigma + eps);

        if (tSigma < varianceThreshold || sSigma < varianceThreshold)
        {
            // No image-based force in featureless regions
            // Rely purely on elastic regularization
            s << vec3d(0, 0, 0);
            continue;
        }

        // Calc others
        double C = 0.0;
        vec3d dC(0.0, 0.0, 0.0);
        for(int index = 0; index < tNeighborhood.size(); ++index)
        {
            double tDiff = tmap.value(tNeighborhood[index]) - tBar;
            double sDiff = smap.value(sNeighborhood[index]) - sBar;

            vec3d gradDiff = smap.gradient(sNeighborhood[index]) - gradSBar;

            C += tDiff * sDiff;

            dC += gradDiff * tDiff;
        }

        // double sDiffSum = 0.0;
        vec3d B(0.0, 0.0, 0.0);
        for(int index = 0; index < tNeighborhood.size(); ++index)
        {
            B += (smap.gradient(sNeighborhood[index]) - gradSBar)*(smap.value(sNeighborhood[index]) - sBar);
        }

        s << (dC - B * C / (sSigma * sSigma + eps)) / (2 * tSigma * sSigma + eps);
	}

    std::cout << "Number of bad nodes: " << bad << " out of " << N << std::endl;
    std::cout << "Number of neighborhood size mismatches: " << mismatches << " out of " << N << std::endl;
	return true;
}

bool FEPlotDiff::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);

			s << T - S;
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotRawDiff::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap0 = pc->GetRawTemplateMap();
	ImageMap& smap0 = pc->GetRawTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap0.valid(r0) && smap0.valid(rt))
		{
			double T = tmap0.value(r0);
			double S = smap0.value(rt);

			s << S - T;
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotGradT::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;

		if (tmap.valid(r0))
		{
			vec3d g = tmap.gradient(r0);
			s << g;
		}
		else
			s << vec3d(0.0);
	}
	return true;
}

bool FEPlotGradS::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d rt = m.Node(i).m_rt;

		if (smap.valid(rt))
		{
			vec3d g = smap.gradient(rt);
			s << g;
		}
		else
			s << vec3d(0.0);
	}
	return true;
}