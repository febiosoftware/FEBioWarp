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
class FEWarpImageConstraintFactory : public FEBioFactory
{
public:
	FEWarpImageConstraintFactory() : FEBioFactory(FENLCONSTRAINT_ID, "warp-image"){}
	void* Create(FEModel* pfem) { return new FEWarpImageConstraint(pfem); }
};

FEWarpImageConstraintFactory warp_image_factory;

//-----------------------------------------------------------------------------
class FEWarpMeshConstraintFactory : public FEBioFactory
{
public:
	FEWarpMeshConstraintFactory() : FEBioFactory(FENLCONSTRAINT_ID, "warp-mesh"){}
	void* Create(FEModel* pfem) { return new FEWarpSurfaceConstraint(pfem); }
};

FEWarpMeshConstraintFactory warp_mesh_factory;

//-----------------------------------------------------------------------------
class FEPlotTemplateFactory : public FEBioFactory
{
public:
	FEPlotTemplateFactory() : FEBioFactory(FEPLOTDATA_ID, "template"){}
	void* Create(FEModel* pfem) { return new FEPlotTemplate(pfem); }
};

FEPlotTemplateFactory plot_template_factory;

//-----------------------------------------------------------------------------
class FEPlotTargetFactory : public FEBioFactory
{
public:
	FEPlotTargetFactory() : FEBioFactory(FEPLOTDATA_ID, "target"){}
	void* Create(FEModel* pfem) { return new FEPlotTarget(pfem); }
};

FEPlotTargetFactory plot_target_factory;

//-----------------------------------------------------------------------------
class FEPlotEnergyFactory : public FEBioFactory
{
public:
	FEPlotEnergyFactory() : FEBioFactory(FEPLOTDATA_ID, "energy"){}
	void* Create(FEModel* pfem) { return new FEPlotEnergy(pfem); }
};

FEPlotEnergyFactory plot_energy_factory;

//-----------------------------------------------------------------------------
class FEPlotForceFactory : public FEBioFactory
{
public:
	FEPlotForceFactory() : FEBioFactory(FEPLOTDATA_ID, "warp force"){}
	void* Create(FEModel* pfem) { return new FEPlotForce(pfem); }
};

FEPlotForceFactory plot_force_factory;

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT void PluginInitialize(FEBioKernel& febio)
{

}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT void PluginCleanup()
{

}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT int PluginNumClasses()
{
	return 6;
}

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT FEBioFactory* PluginGetFactory(int i)
{
	switch (i)
	{
	case 0: return &warp_image_factory;
	case 1: return &warp_mesh_factory;
	case 2: return &plot_template_factory;
	case 3: return &plot_target_factory;
	case 4: return &plot_energy_factory;
	case 5: return &plot_force_factory;
	default:
		return 0;
	}
}
