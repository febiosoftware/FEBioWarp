#pragma once
#include "FEWarpVolumeConstraint.h"
#include <FEImgLib/FEImageSource.h>
#include <FEImgLib/ImageMap.h>

//-----------------------------------------------------------------------------
class FEWarpImageConstraint : public FEWarpVolumeConstraint
{
public:
	FEWarpImageConstraint(FEModel* pfem);
	~FEWarpImageConstraint();

	// initialization
	bool Init();

	// update
	void Update();

protected:
	double	m_blur;		//!< blur factor
	
	double		m_blur_cur;	//!< current blur factor
	int 		m_blur_method;
};

//-----------------------------------------------------------------------------
class FEWarpSingleImageConstraint : public FEWarpImageConstraint
{
public:
	FEWarpSingleImageConstraint(FEModel* fem);

	// initialization
	bool Init() override;

private:
	FEImageSource* m_tmpReader;
	FEImageSource* m_trgReader;

	DECLARE_FECORE_CLASS();
};
