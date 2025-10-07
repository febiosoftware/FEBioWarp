#pragma once
#include "FEWarpVolumeConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>
#include <FEImgLib/ImageFilter.h>
#include <map>

//-----------------------------------------------------------------------------
class FEWarpFilteredImageConstraint : public FEWarpVolumeConstraint
{
public:
	FEWarpFilteredImageConstraint(FEModel* pfem);
	~FEWarpFilteredImageConstraint();

	// initialization
	bool Init();

	// update
	void Update();

protected:
	double		m_pt;	//!< previous time
	double		m_pr;   //!< previous blur radius	

	ImageFilter* m_filt = nullptr;
};

//-----------------------------------------------------------------------------
class FEWarpSingleFilteredImageConstraint : public FEWarpFilteredImageConstraint
{
public:
	FEWarpSingleFilteredImageConstraint(FEModel* fem);

	// initialization
	bool Init() override;

private:
	FEImageSource* m_tmpReader;
	FEImageSource* m_trgReader;
	
	DECLARE_FECORE_CLASS();
};
