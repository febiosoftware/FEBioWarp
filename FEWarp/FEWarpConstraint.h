#pragma once
#include <FECore/FENLConstraint.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FEBodyForce.h>
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
	void Init();
	void Residual(FEGlobalVector& R);
	void StiffnessMatrix(FENLSolver* psolver);
	bool Augment(int naug);
	void Serialize(DumpFile& ar);
	void Update();

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

	vector<vec3d>	m_Lm;	//!< Lagrange multipliers
	int				m_nint;	//!< running counter for multipliers
};

//-----------------------------------------------------------------------------
class FEWarpImageConstraint : public FEWarpConstraint
{
public:
	FEWarpImageConstraint(FEModel* pfem);
	~FEWarpImageConstraint() {}

	// initialization
	void Init();

public:
	ImageMap& GetTemplateMap() { return m_tmap; }
	ImageMap& GetTargetMap  () { return m_smap; }

public:
	//! Calculate the force at a material point
	vec3d wrpForce(FEMaterialPoint& pt);

	//! calculate the stiffness at a material point
	mat3ds wrpStiffness(FEMaterialPoint& pt);

protected:
	double	m_r0[3];	//!< minimum range
	double	m_r1[3];	//!< maximum range

	Image	m_tmp;	//!< template image
	Image	m_trg;	//!< target image

	ImageMap m_tmap;
	ImageMap m_smap;

public:
	DECLARE_PARAMETER_LIST();
};