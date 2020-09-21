#include "stdafx.h"
#include "FEWarpConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>
#include <iostream>
using namespace std;

extern FECoreKernel* pFEBio;

//=============================================================================

FEWarpConstraint::FEWarpConstraint(FEModel* pfem) : FENLConstraint(pfem)
{
	m_blaugon = false;
	m_altol = 0.1;
	m_k = 0;
}

//-----------------------------------------------------------------------------
FEWarpConstraint::~FEWarpConstraint(void)
{

}

//-----------------------------------------------------------------------------
bool FEWarpConstraint::Init()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// if no domains are selected, add all domains
	if (m_dom.empty()) 
	{
		for (int i=0; i<mesh.Domains(); ++i) m_dom.push_back(i);
	}

	// figure out how many integration points we need
	int nint = 0;
	int ND = m_dom.size();
	for (int i=0; i<ND; ++i)
	{
		FESolidDomain* dom = dynamic_cast<FESolidDomain*>(&mesh.Domain(m_dom[i]));
		if (dom)
		{
			int NE = dom->Elements();
			for (int j=0; j<NE; ++j)
			{
				FESolidElement& el = dom->Element(j);
				nint += el.GaussPoints();
			}
		}
	}

	// allocate storage for Lagrange multipliers
	m_Lm.assign(nint, vec3d(0,0,0));

	return true;
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	vector<double> fe;
	vector<int> lm;

	// reset multiplier counter
	m_nint = 0;

	// loop over all domains
	int NDOM = m_dom.size();
	for (int n=0; n<NDOM; ++n)
	{
		FESolidDomain* dom = dynamic_cast<FESolidDomain*>(&mesh.Domain(m_dom[n]));
		if (dom)
		{
			int NEL = dom->Elements();
			for (int i=0; i<NEL; ++i)
			{
				FESolidElement& el = dom->Element(i);

				// get the element force vector and initialize it to zero
				int ndof = 3*el.Nodes();
				fe.assign(ndof, 0);

				// apply body forces
				ElementWarpForce(*dom, el, fe);

				// get the element's LM vector
				dom->UnpackLM(el, lm);

				// assemble element 'fe'-vector into global R vector
				R.Assemble(el.m_node, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::ElementWarpForce(FESolidDomain& dom, FESolidElement& el, vector<double>& fe)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	vec3d f;

	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i)
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
	}

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);

		detJ = dom.detJ0(el, n)*gw[n];

		// get the force
		f = m_Lm[m_nint++] + wrpForce(mp);

		H = el.H(n);

		// don't forget to multiply with the density
		double dens = pme->Density(mp);

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] -= H[i]*dens*f.x*detJ;
			fe[3*i+1] -= H[i]*dens*f.y*detJ;
			fe[3*i+2] -= H[i]*dens*f.z*detJ;
		}						
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// element stiffness matrix
	vector<int> lm;

	// loop over all domains
	int NDOM = m_dom.size();
	for (int n=0; n<NDOM; ++n)
	{
		FESolidDomain* dom = dynamic_cast<FESolidDomain*>(&mesh.Domain(m_dom[n]));
		if (dom)
		{
			// repeat over all solid elements
			int NE = dom->Elements();
			for (int iel=0; iel<NE; ++iel)
			{
				FESolidElement& el = dom->Element(iel);

				// create the element's stiffness matrix
				FEElementMatrix ke;
				int ndof = 3*el.Nodes();
				ke.resize(ndof, ndof);
				ke.zero();

				// calculate inertial stiffness
				ElementWarpStiffness(*dom, el, ke);

				// get the element's LM vector
				dom->UnpackLM(el, lm);

				// assemble element matrix in global stiffness matrix
				ke.SetIndices(lm);
				ke.SetNodes(el.m_node);
				LS.Assemble(ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::ElementWarpStiffness(FESolidDomain& dom, FESolidElement& el, matrix& ke)
{
	int neln = el.Nodes();
	int ndof = ke.columns()/neln;

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	mat3ds K;

	FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		detJ = dom.detJ0(el, n)*gw[n];

		// get the stiffness
		K = wrpStiffness(mp);

		H = el.H(n);

		double dens = pme->Density(mp);

		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				ke[ndof*i  ][ndof*j  ] -= H[i]*H[j]*dens*K(0,0)*detJ;
				ke[ndof*i  ][ndof*j+1] -= H[i]*H[j]*dens*K(0,1)*detJ;
				ke[ndof*i  ][ndof*j+2] -= H[i]*H[j]*dens*K(0,2)*detJ;

				ke[ndof*i+1][ndof*j  ] -= H[i]*H[j]*dens*K(1,0)*detJ;
				ke[ndof*i+1][ndof*j+1] -= H[i]*H[j]*dens*K(1,1)*detJ;
				ke[ndof*i+1][ndof*j+2] -= H[i]*H[j]*dens*K(1,2)*detJ;

				ke[ndof*i+2][ndof*j  ] -= H[i]*H[j]*dens*K(2,0)*detJ;
				ke[ndof*i+2][ndof*j+1] -= H[i]*H[j]*dens*K(2,1)*detJ;
				ke[ndof*i+2][ndof*j+2] -= H[i]*H[j]*dens*K(2,2)*detJ;
			}
	}	
}

//-----------------------------------------------------------------------------
bool FEWarpConstraint::Augment(int naug, const FETimeInfo& tp)
{
	if (m_blaugon == false) return true;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	vector<vec3d> L0(m_Lm);
	vector<vec3d> L1(m_Lm);

	m_nint = 0;
	int NDOM = m_dom.size();
	for (int i=0; i<NDOM; ++i)
	{
		FESolidDomain* dom = dynamic_cast<FESolidDomain*>(&mesh.Domain(m_dom[i]));
		if (dom)
		{
			int NE = dom->Elements();
			for (int j=0; j<NE; ++j)
			{
				FESolidElement& el = dom->Element(j);
				int nint = el.GaussPoints();
				int neln = el.Nodes();

				// nodal coordinates
				vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
				for (int i=0; i<neln; ++i)
				{
					r0[i] = mesh.Node(el.m_node[i]).m_r0;
					rt[i] = mesh.Node(el.m_node[i]).m_rt;
				}

				for (int n=0; n<nint; ++n, ++m_nint)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
					pt.m_r0 = el.Evaluate(r0, n);
					pt.m_rt = el.Evaluate(rt, n);

					L1[m_nint] = L0[m_nint] + wrpForce(pt);
				}
			}
		}
	}

	// calculate the norm
	double normL0 = 0, normL1 = 0;
	for (int i=0; i<m_nint; ++i) 
	{
		normL0 += L0[i]*L0[i];
		normL1 += L1[i]*L1[i];
	}

	double Lerr = fabs((normL1 - normL0)/normL1);

	feLog("warping norm: %lg %lg\n", Lerr, m_altol);

	if (Lerr >= m_altol)
	{
		// update multipliers
		m_Lm = L1;
	}

	// check convergence
	return (Lerr < m_altol);
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Update()
{

}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);

	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving())
		{
			ar << m_dom;
			ar << m_Lm;
			ar << m_nint;
		}
		else
		{
			ar >> m_dom;
			ar >> m_Lm;
			ar >> m_nint;
		}
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	// nothing to do here.
}
