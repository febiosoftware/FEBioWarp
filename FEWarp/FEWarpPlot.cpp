#include "StdAfx.h"
#include "FECore/FESolidDomain.h"
#include "FEWarpPlot.h"
#include "ImageMap.h"
#include "FEWarpConstraint.h"
#include "FEWarpImageConstraint.h"
#include <FECore/FEModel.h>

extern ImageMap ds1map;
extern ImageMap ds2map;
extern ImageMap ds3map;

bool FEPlotTemplate::Save(FEMesh &m, std::vector<float> &a)
{
	FEModel& fem = *m_pfem;
	a.assign(m.Nodes(), 0.f);
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);
		FEWarpImageConstraint* pci = dynamic_cast<FEWarpImageConstraint*>(pc);
		if (pci) { return SaveWarpImage(m, pci, a); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs) { return SaveWarpMesh(m, pcs, a); }
	}
	return true;
}

bool FEPlotTemplate::SaveWarpImage(FEMesh& m, FEWarpImageConstraint* pc, vector<float>& a)
{
	// get the template image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i) a[i] = (float) tmap.value(m.Node(i).m_r0);
	return true;
}

bool FEPlotTemplate::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, vector<float>& a)
{
	FEWarpSurface* ps = pc->GetTemplate();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		a[i] = (float) ps->value(m.Node(i).m_r0);
	}
	return true;
}

bool FEPlotTarget::Save(FEMesh &m, std::vector<float> &a)
{
	FEModel& fem = *m_pfem;
	a.assign(m.Nodes(), 0.f);
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);
		FEWarpImageConstraint* pci = dynamic_cast<FEWarpImageConstraint*>(pc);
		if (pci) { return SaveWarpImage(m, pci, a); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs) { return SaveWarpMesh(m, pcs, a); }
	}
	return true;
}


bool FEPlotTarget::SaveWarpImage(FEMesh& m, FEWarpImageConstraint* pc, vector<float>& a)
{
	// get the target image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i) a[i] = (float) smap.value(m.Node(i).m_rt);
	return true;
}

bool FEPlotTarget::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, vector<float>& a)
{
	FEWarpSurface* ps = pc->GetTarget();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		a[i] = (float) ps->value(m.Node(i).m_rt);
	}
	return true;
}

bool FEPlotEnergy::Save(FEMesh &m, std::vector<float> &a)
{
	// find the warping constraint
	FEModel& fem = *m_pfem;
	FEWarpImageConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpImageConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		double T = tmap.value(m.Node(i).m_r0);
		double S = smap.value(m.Node(i).m_rt);
		a.push_back((float) (0.5*(T-S)*(T-S)));
	}
	return true;
}

bool FEPlotForce::Save(FEMesh &m, std::vector<float> &a)
{
	// find the warping constraint
	FEModel& fem = *m_pfem;
	FEWarpImageConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpImageConstraint*>(fem.NonlinearConstraint(i));
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

		double T = tmap.value(r0);

		ImageMap::POINT ps = smap.map(rt);

		double S = smap.value(ps);

		double dx = smap.dx();
		double dy = smap.dy();
		double dz = smap.dz();

		double dSx = ds1map.value(ps)/dx;
		double dSy = ds2map.value(ps)/dy;
		double dSz = ds3map.value(ps)/dz;

		vec3d fw = (vec3d(dSx, dSy, dSz))*((T - S));

		a.push_back((float) fw.x);
		a.push_back((float) fw.y);
		a.push_back((float) fw.z);
	}
	return true;
}
