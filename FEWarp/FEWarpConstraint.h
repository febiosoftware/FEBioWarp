#pragma once
#include <FECore/FENLConstraint.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FESolidDomain.h>
#include <FECore/Image.h>
#include "ImageMap.h"

//-----------------------------------------------------------------------------
// This class implements a warping constraint. This is a image-based constraint
// where the goal is to map a template image onto a target image.
class FEWarpConstraint : public FENLConstraint
{
public:
	FEWarpConstraint(FEModel* pfem);
	~FEWarpConstraint(void);

public: // inherited members
	bool Init();
	void Residual(FEGlobalVector& R, const FETimePoint& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimePoint& tp);
	bool Augment(int naug, const FETimePoint& tp);
	void Serialize(DumpFile& ar);
	void Update();
	void BuildMatrixProfile(FEGlobalMatrix& M);

protected:
	void ElementWarpForce    (FESolidDomain& dom, FESolidElement& el, vector<double>& fe, double dens);
	void ElementWarpStiffness(FESolidDomain& dom, FESolidElement& el, matrix& ke, double dens);

	//! Calculate the force at a material point
	virtual vec3d wrpForce(FEMaterialPoint& pt) = 0;

	//! calculate the stiffness at a material point
	virtual mat3ds wrpStiffness(FEMaterialPoint& pt) = 0;

protected:
	bool		m_blaugon;	//!< augmented lagrangian flag
	double		m_altol;	//!< augmentation tolerance
	double		m_k;		//!< penalty parameter

	vector<int>	m_dom;		//!< list of domains to warp

	vector<vec3d>	m_Lm;	//!< Lagrange multipliers
	int				m_nint;	//!< running counter for multipliers
};
