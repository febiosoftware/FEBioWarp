// FEWarp.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModelUpdate.h>
#include "FEWarpImageConstraint.h"
#include "FEWarpImageConstraint2.h"
#include "FEWarpSurfaceConstraint.h"
#include "FEWarpPlot.h"
#include "FEWarpLog.h"

//-----------------------------------------------------------------------------
void UpdateWarpModel(FECoreBase* pc, FEModelUpdate& fem)
{
	if (dynamic_cast<FEWarpImageConstraint*>(pc))
	{
		fem.AddPlotVariable("warp-template");
	}
}

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
//	REGISTER_FECORE_CLASS(FEWarpImageConstraint2  , "warp-image2");
	REGISTER_FECORE_CLASS(FEWarpSurfaceConstraint, "warp-mesh" );

	REGISTER_FECORE_CLASS(FEPlotTemplate, "warp-template");
	REGISTER_FECORE_CLASS(FEPlotTarget  , "warp-target"  );
	REGISTER_FECORE_CLASS(FEPlotEnergy  , "warp-energy"  );
	REGISTER_FECORE_CLASS(FEPlotForce   , "warp-force"   );

	REGISTER_FECORE_CLASS(FELogWarpTemplate, "warp-template");
	REGISTER_FECORE_CLASS(FELogWarpTarget  , "warp-target");
	REGISTER_FECORE_CLASS(FELogWarpEnergy  , "warp-energy");
	REGISTER_FECORE_CLASS(FELogWarpForceX  , "warp-force-x");
	REGISTER_FECORE_CLASS(FELogWarpForceY  , "warp-force-y");
	REGISTER_FECORE_CLASS(FELogWarpForceZ  , "warp-force-z");

	// model updates
	febio.OnCreateEvent(UpdateModelWhenCreating<FEWarpImageConstraint>([](FEModelUpdate& fem) {
		fem.AddPlotVariable("warp-template");
		fem.AddPlotVariable("warp-target");
		fem.AddPlotVariable("warp-energy");
		fem.AddPlotVariable("warp-force");
		})
	);
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
