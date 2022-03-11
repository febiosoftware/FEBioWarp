#include "stdafx.h"
#include "FEWarpPlot.h"
#include "ImageMap.h"
#include "FEWarpConstraint.h"
#include "FEWarpImageConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FEDataStream.h>
#include <FECore/FEAnalysis.h>

bool FEPlotTemplate::Save(FEMesh &m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* step = fem.GetCurrentStep();
	if (step == nullptr) return false;
	for (int i=0; i<step->StepComponents(); ++i)
	{
		FEWarpImageConstraint* pci = dynamic_cast<FEWarpImageConstraint*>(step->GetStepComponent(i));
		if (pci) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(step->GetStepComponent(i));
		if (pcs) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}

bool FEPlotTemplate::SaveWarpImage(FEMesh& m, FEWarpImageConstraint* pc, FEDataStream& s)
{
	// get the template image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i) s << tmap.value(m.Node(i).m_r0);
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
	FEAnalysis* step = fem.GetCurrentStep();
	if (step == nullptr) return false;
	for (int i = 0; i < step->StepComponents(); ++i)
	{
		FEWarpImageConstraint* pci = dynamic_cast<FEWarpImageConstraint*>(step->GetStepComponent(i));
		if (pci) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(step->GetStepComponent(i));
		if (pcs) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}


bool FEPlotTarget::SaveWarpImage(FEMesh& m, FEWarpImageConstraint* pc, FEDataStream& s)
{
	// get the target image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i) s << smap.value(m.Node(i).m_rt);
	return true;
}

bool FEPlotTarget::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s)
{
	FEWarpSurface* ps = pc->GetTarget();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		s << ps->value(m.Node(i).m_rt);
	}
	return true;
}

bool FEPlotEnergy::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
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

		s << (0.5*(T - S)*(T - S));
	}
	return true;
}

bool FEPlotForce::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
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
		double S = smap.value(rt);
		vec3d G = smap.gradient(rt);
		vec3d fw = G*((T - S));

		s << fw;
	}
	return true;
}
