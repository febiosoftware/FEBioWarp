#include "stdafx.h"
#include "FEWarpPlot.h"
#include "FEWarpConstraint.h"
#include "FEWarpVolumeConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FEDataStream.h>
#include <FECore/FEAnalysis.h>

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

bool FEPlotForce::Save(FEMesh &m, FEDataStream& s)
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
			vec3d G = smap.gradient(rt);
			vec3d fw = G * ((T - S));

			s << fw;
		}
		else 
			s << vec3d(0.0);
	}
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

			s << S - T;
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