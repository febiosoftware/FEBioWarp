// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FECore/febio.h"
#include "FECore/FEBioFactory.h"
#include "FEWarpImageConstraint.h"
#include "FEWarpSurfaceConstraint.h"
#include "FEWarpPlot.h"

#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

FEBioKernel* pFEBio;

//-----------------------------------------------------------------------------
class FEWarpImageConstraintFactory : public FEBioFactory_T<FENLConstraint>
{
public:
	FEWarpImageConstraintFactory() : FEBioFactory_T<FENLConstraint>("warp-image"){}
	bool IsType(FENLConstraint* pf) { return (dynamic_cast<FEWarpImageConstraint*>(pf) != 0); }
	FENLConstraint* Create(FEModel* pfem) { return new FEWarpImageConstraint(pfem); }
};

FEWarpImageConstraintFactory warp_image_factory;

//-----------------------------------------------------------------------------
class FEWarpMeshConstraintFactory : public FEBioFactory_T<FENLConstraint>
{
public:
	FEWarpMeshConstraintFactory() : FEBioFactory_T<FENLConstraint>("warp-mesh"){}
	bool IsType(FENLConstraint* pf) { return (dynamic_cast<FEWarpSurfaceConstraint*>(pf) != 0); }
	FENLConstraint* Create(FEModel* pfem) { return new FEWarpSurfaceConstraint(pfem); }
};

FEWarpMeshConstraintFactory warp_mesh_factory;


//-----------------------------------------------------------------------------
class FEPlotTemplateFactory : public FEBioFactory_T<FEPlotData>
{
public:
	FEPlotTemplateFactory() : FEBioFactory_T<FEPlotData>("template"){}
	bool IsType(FEPlotData* po) { return (dynamic_cast<FEPlotTemplate*>(po) != 0); }
	FEPlotData* Create(FEModel* pfem) { return new FEPlotTemplate(pfem); }
};

FEPlotTemplateFactory plot_template_factory;

//-----------------------------------------------------------------------------
class FEPlotTargetFactory : public FEBioFactory_T<FEPlotData>
{
public:
	FEPlotTargetFactory() : FEBioFactory_T<FEPlotData>("target"){}
	bool IsType(FEPlotData* po) { return (dynamic_cast<FEPlotTarget*>(po) != 0); }
	FEPlotData* Create(FEModel* pfem) { return new FEPlotTarget(pfem); }
};

FEPlotTargetFactory plot_target_factory;

//-----------------------------------------------------------------------------
class FEPlotEnergyFactory : public FEBioFactory_T<FEPlotData>
{
public:
	FEPlotEnergyFactory() : FEBioFactory_T<FEPlotData>("energy"){}
	bool IsType(FEPlotData* po) { return (dynamic_cast<FEPlotEnergy*>(po) != 0); }
	FEPlotData* Create(FEModel* pfem) { return new FEPlotEnergy(pfem); }
};

FEPlotEnergyFactory plot_energy_factory;

//-----------------------------------------------------------------------------
class FEPlotForceFactory : public FEBioFactory_T<FEPlotData>
{
public:
	FEPlotForceFactory() : FEBioFactory_T<FEPlotData>("warp force"){}
	bool IsType(FEPlotData* po) { return (dynamic_cast<FEPlotForce*>(po) != 0); }
	FEPlotData* Create(FEModel* pfem) { return new FEPlotForce(pfem); }
};

FEPlotForceFactory plot_force_factory;

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT FEBioFactory* RegisterPlugin(FEBioKernel& febio)
{
	static int n = -1;
	n++;

	switch (n)
	{
	case 0: return &warp_image_factory;
	case 1: return &warp_mesh_factory;
	case 2: return &plot_template_factory;
	case 3: return &plot_target_factory;
	case 4: return &plot_energy_factory;
	case 5: return &plot_force_factory;
	default:
		pFEBio = &febio;
		return 0;
	}
}
