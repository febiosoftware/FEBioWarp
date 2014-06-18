// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FECore/FECoreKernel.h"
#include "FEWarpImageConstraint.h"
#include "FEWarpSurfaceConstraint.h"
#include "FEWarpPlot.h"

FECoreKernel* pFEBio;

//-----------------------------------------------------------------------------
//! Factory classes for this plugin.
FEPluginFactory_T<FEWarpImageConstraint  , FENLCONSTRAINT_ID> warp_image_factory   ("warp-image");
FEPluginFactory_T<FEWarpSurfaceConstraint, FENLCONSTRAINT_ID> warp_mesh_factory    ("warp-mesh" );
FEPluginFactory_T<FEPlotTemplate         , FEPLOTDATA_ID    > plot_template_factory("warp-template");
FEPluginFactory_T<FEPlotTarget           , FEPLOTDATA_ID    > plot_target_factory  ("warp-target"  );
FEPluginFactory_T<FEPlotEnergy           , FEPLOTDATA_ID    > plot_energy_factory  ("warp-energy"  );
FEPluginFactory_T<FEPlotForce            , FEPLOTDATA_ID    > plot_force_factory   ("warp-force"   );

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{

}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginCleanup()
{

}

//-----------------------------------------------------------------------------
FECORE_EXPORT int PluginNumClasses()
{
	return 6;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT FECoreFactory* PluginGetFactory(int i)
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
