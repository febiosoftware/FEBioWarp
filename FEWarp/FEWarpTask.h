#pragma once
#include <FECore/FECoreTask.h>

class FEWarpTask : public FECoreTask
{
public:
	FEWarpTask(FEModel* fem);

	bool Init(const char* szfile) override;

	bool Run() override;
};
