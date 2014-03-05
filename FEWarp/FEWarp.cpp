// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FECoreFactory.h"
#include "FEWarpImageConstraint.h"
#include "FEWarpSurfaceConstraint.h"
#include "FEWarpPlot.h"

#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

FECoreKernel* pFEBio;

//-----------------------------------------------------------------------------
class FEWarpImageConstraintFactory : public FECoreFactory
{
public:
	FEWarpImageConstraintFactory() : FECoreFactory(FENLCONSTRAINT_ID, "warp-image"){}
	void* Create(FEModel* pfem) { return new FEWarpImageConstraint(pfem); }
};

FEWarpImageConstraintFactory warp_image_factory;

//-----------------------------------------------------------------------------
class FEWarpMeshConstraintFactory : public FECoreFactory
{
public:
	FEWarpMeshConstraintFactory() : FECoreFactory(FENLCONSTRAINT_ID, "warp-mesh"){}
	void* Create(FEModel* pfem) { return new FEWarpSurfaceConstraint(pfem); }
};

FEWarpMeshConstraintFactory warp_mesh_factory;

//-----------------------------------------------------------------------------
class FEPlotTemplateFactory : public FECoreFactory
{
public:
	FEPlotTemplateFactory() : FECoreFactory(FEPLOTDATA_ID, "template"){}
	void* Create(FEModel* pfem) { return new FEPlotTemplate(pfem); }
};

FEPlotTemplateFactory plot_template_factory;

//-----------------------------------------------------------------------------
class FEPlotTargetFactory : public FECoreFactory
{
public:
	FEPlotTargetFactory() : FECoreFactory(FEPLOTDATA_ID, "target"){}
	void* Create(FEModel* pfem) { return new FEPlotTarget(pfem); }
};

FEPlotTargetFactory plot_target_factory;

//-----------------------------------------------------------------------------
class FEPlotEnergyFactory : public FECoreFactory
{
public:
	FEPlotEnergyFactory() : FECoreFactory(FEPLOTDATA_ID, "energy"){}
	void* Create(FEModel* pfem) { return new FEPlotEnergy(pfem); }
};

FEPlotEnergyFactory plot_energy_factory;

//-----------------------------------------------------------------------------
class FEPlotForceFactory : public FECoreFactory
{
public:
	FEPlotForceFactory() : FECoreFactory(FEPLOTDATA_ID, "warp force"){}
	void* Create(FEModel* pfem) { return new FEPlotForce(pfem); }
};

FEPlotForceFactory plot_force_factory;

//-----------------------------------------------------------------------------
extern "C" DLL_EXPORT void PluginInitialize(FECoreKernel& febio)
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
extern "C" DLL_EXPORT FECoreFactory* PluginGetFactory(int i)
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
