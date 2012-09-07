// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FECore/febio.h"
#include "FECore/FEBioFactory.h"
#include "FEWarpForce.h"
#include "FEWarpPlot.h"

FEBioKernel* pFEBio;

//-----------------------------------------------------------------------------
class FEWarpForceFactory : public FEBioFactory_T<FEBodyForce>
{
public:
	FEWarpForceFactory() : FEBioFactory_T<FEBodyForce>("warp"){}
	bool IsType(FEBodyForce* pf) { return (dynamic_cast<FEWarpForceFactory*>(pf) != 0); }
	FEBodyForce* Create(FEModel* pfem) { return new FEWarpForce(pfem); }
};

FEWarpForceFactory warp_force_factory;

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
extern "C" __declspec(dllexport) void RegisterPlugin(FEBioKernel& febio)
{
	febio.RegisterClass(&warp_force_factory);
	febio.RegisterClass(&plot_template_factory);
	febio.RegisterClass(&plot_target_factory);
	febio.RegisterClass(&plot_energy_factory);
	febio.RegisterClass(&plot_force_factory);
	pFEBio = &febio;
}
