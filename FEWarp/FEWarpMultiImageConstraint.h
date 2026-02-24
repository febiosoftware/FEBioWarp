#pragma once
#include "FEWarpImageConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>

//-----------------------------------------------------------------------------
class FEWarpMultiImageConstraint : public FEWarpImageConstraint
{
public:
	FEWarpMultiImageConstraint(FEModel* pfem);

	// initialization
	bool Init() override;

	// update
	void Update() override;

protected:
	FEImageSource* m_tmpReader;
	std::vector<FEImageSource*> m_trgReader;
	int		m_ntrg;		//!< current target

public:
	DECLARE_FECORE_CLASS();
};
