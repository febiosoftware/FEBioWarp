#pragma once
#include <FECore/FEBodyConstraint.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
// This class implements a warping constraint. This is a image-based constraint
// where the goal is to map a template image onto a target image.
class FEWarpConstraint : public FEBodyConstraint
{
public:
	FEWarpConstraint(FEModel* pfem);
	~FEWarpConstraint(void);

public: // inherited members
	bool Init();
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void Serialize(DumpStream& ar) override;
	void Update() override;
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

protected:
	void ElementWarpForce    (FESolidDomain& dom, FESolidElement& el, vector<double>& fe);
	void ElementWarpStiffness(FESolidDomain& dom, FESolidElement& el, matrix& ke);

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
