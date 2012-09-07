#pragma once
#include "FECore/FEBodyForce.h"

class FEWarpForce : public FEBodyForce
{
public:
	FEWarpForce(FEModel* pfem);
	~FEWarpForce(void);

	vec3d force(FEMaterialPoint& pt);

	mat3ds stiffness(FEMaterialPoint& pt);

	bool Init();

public:
	char	m_sztmp[256];	// template image
	char	m_sztrg[256];	// target image

	double	m_k;

	DECLARE_PARAMETER_LIST();
};
