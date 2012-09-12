#pragma once
#include <FECore/FENLConstraint.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FEBodyForce.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
// This class implements the actual image force
class FEWarpForce : public FEBodyForce
{
public:
	FEWarpForce(FEModel* pfem);

	//! Initialize force data
	bool Init();

	//! Calculate the force at a material point
	vec3d force(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	mat3ds stiffness(FEMaterialPoint& pt);

public:
	char	m_sztmp[256];	// template image
	char	m_sztrg[256];	// target image

	double	m_k;	// penalty parameter
};

//-----------------------------------------------------------------------------
// This class implements a warping constraint. This is a image-based constraint
// where the goal is to map a template image onto a target image.
class FEWarpConstraint : public FENLConstraint
{
public:
	FEWarpConstraint(FEModel* pfem);
	~FEWarpConstraint(void);

public: // inherited members
	void Init();
	void Residual(FENLSolver* psolver, vector<double>& R);
	void StiffnessMatrix(FENLSolver* psolver);
	bool Augment(int naug);
	void Serialize(DumpFile& ar);
	void Update();

protected:
	void ElementWarpForce    (FESolidDomain& dom, FESolidElement& el, vector<double>& fe, double dens);
	void ElementWarpStiffness(FESolidDomain& dom, FESolidElement& el, matrix& ke, double dens);

protected:
	FEWarpForce	m_wrp;		//!< warp force class
	bool		m_blaugon;	//!< augmented lagrangian flag
	double		m_altol;	//!< augmentation tolerance

public:
	DECLARE_PARAMETER_LIST();
};
