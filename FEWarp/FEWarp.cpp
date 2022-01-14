// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "FECore/FECoreKernel.h"
#include "FEWarpImageConstraint.h"
#include "FEWarpSurfaceConstraint.h"
#include "FEWarpPlot.h"

//-----------------------------------------------------------------------------
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);

	// we're extending the solid module
	febio.SetActiveModule("solid");

	// constraints
	REGISTER_FECORE_CLASS(FEWarpImageConstraint  , "warp-image");
	REGISTER_FECORE_CLASS(FEWarpSurfaceConstraint, "warp-mesh" );

	REGISTER_FECORE_CLASS(FEPlotTemplate, "warp-template");
	REGISTER_FECORE_CLASS(FEPlotTarget  , "warp-target"  );
	REGISTER_FECORE_CLASS(FEPlotEnergy  , "warp-energy"  );
	REGISTER_FECORE_CLASS(FEPlotForce   , "warp-force"   );
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void GetPluginVersion(int& maj, int& min, int& patch)
{
	maj = 2;
	min = 0;
	patch = 0;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginCleanup()
{

}
