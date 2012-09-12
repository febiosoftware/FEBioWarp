#include "StdAfx.h"
#include "FEWarpConstraint.h"
#include "Image.h"
#include "ImageMap.h"
#include <FECore/FEElasticMaterial.h>
#include <FECore/FEModel.h>
#include <iostream>
using namespace std;

BEGIN_PARAMETER_LIST(FEWarpConstraint, FENLConstraint);
	ADD_PARAMETER(m_wrp.m_sztmp, FE_PARAM_STRING, "template");
	ADD_PARAMETER(m_wrp.m_sztrg, FE_PARAM_STRING, "target"  );
	ADD_PARAMETER(m_wrp.m_k    , FE_PARAM_DOUBLE, "penalty" );
	ADD_PARAMETER(m_blaugon    , FE_PARAM_BOOL  , "laugon"  );
	ADD_PARAMETER(m_altol      , FE_PARAM_DOUBLE, "altol"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
Image tmp;
Image trg;

Image ds[3];	// first derivate of trg
Image dds[6];	// second derivate of trg (assumed symmetric: xx, yy, zz, xy, yz, xz)

ImageMap tmap(tmp);
ImageMap smap(trg);

ImageMap ds1map(ds[0]);
ImageMap ds2map(ds[1]);
ImageMap ds3map(ds[2]);

ImageMap dds1map(dds[0]);
ImageMap dds2map(dds[1]);
ImageMap dds3map(dds[2]);
ImageMap dds4map(dds[3]);
ImageMap dds5map(dds[4]);
ImageMap dds6map(dds[5]);

//-----------------------------------------------------------------------------
FEWarpForce::FEWarpForce(FEModel* pfem) : FEBodyForce(pfem)
{
	m_k = 0;
}

//-----------------------------------------------------------------------------
bool FEWarpForce::Init()
{
	// allocate storage for image data
	tmp.Create(64, 64, 1);
	trg.Create(64, 64, 1);
	ds[0].Create(64, 64, 1);
	ds[1].Create(64, 64, 1);
	ds[2].Create(64, 64, 1);

	dds[0].Create(64, 64, 1);
	dds[1].Create(64, 64, 1);
	dds[2].Create(64, 64, 1);
	dds[3].Create(64, 64, 1);
	dds[4].Create(64, 64, 1);
	dds[5].Create(64, 64, 1);

	tmap.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	smap.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	ds1map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	ds2map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	ds3map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));

	dds1map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	dds2map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	dds3map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	dds4map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	dds5map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));
	dds6map.SetRange(vec3d(0,0,0), vec3d(64, 64, 0));

	// try to load the image data
	if (tmp.Load(m_sztmp) == false) { cerr << "Failed reading template data " << m_sztmp << endl; return false; }
	else cout << "Template loaded successfully\n";

	if (trg.Load(m_sztrg) == false) { cerr << "Failed reading target data " << m_sztrg << endl; return false; }
	else cout << "Target loaded successfully\n";

	image_derive_x(trg, ds[0]);
	image_derive_y(trg, ds[1]);
	image_derive_z(trg, ds[2]);

	image_derive_x(ds[0], dds[0]);
	image_derive_y(ds[1], dds[1]);
	image_derive_z(ds[2], dds[2]);
	image_derive_y(ds[0], dds[3]);
	image_derive_z(ds[1], dds[4]);
	image_derive_z(ds[0], dds[5]);

	return true;
}

//-----------------------------------------------------------------------------
vec3d FEWarpForce::force(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double T = tmap.value(pt.r0);

	ImageMap::POINT ps = smap.map(pt.rt);

	double S = smap.value(ps);

	double dx = smap.dx();
	double dy = smap.dy();
	double dz = smap.dz();

	double dSx = ds1map.value(ps)/dx;
	double dSy = ds2map.value(ps)/dy;
	double dSz = ds3map.value(ps)/dz;

	vec3d Fw = (vec3d(dSx, dSy, dSz))*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpForce::stiffness(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double T = tmap.value(pt.r0);

	ImageMap::POINT ps = smap.map(pt.rt);

	double S = smap.value(ps);

	double dx = smap.dx();
	double dy = smap.dy();
	double dz = smap.dz();

	double dSx = ds1map.value(ps)/dx;
	double dSy = ds2map.value(ps)/dy;
	double dSz = ds3map.value(ps)/dz;

	vec3d dS = vec3d(dSx, dSy, dSz);

	mat3ds H;
	H.xx() = dds1map.value(ps)/(dx*dx);
	H.yy() = dds2map.value(ps)/(dy*dy);
	H.zz() = dds3map.value(ps)/(dz*dz);
	H.xy() = dds4map.value(ps)/(dx*dy);
	H.yz() = dds5map.value(ps)/(dy*dz);
	H.xz() = dds6map.value(ps)/(dx*dz);

	return H*((T - S)*m_k) - dyad(dS)*m_k;
}

//=============================================================================

FEWarpConstraint::FEWarpConstraint(FEModel* pfem) : FENLConstraint(pfem), m_wrp(pfem)
{
}

//-----------------------------------------------------------------------------
FEWarpConstraint::~FEWarpConstraint(void)
{
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Init()
{
	m_wrp.Init();
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Residual(FENLSolver* psolver, vector<double>& R)
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	vector<double> fe;
	vector<int> lm;

	// loop over all domains
	int NDOM = mesh.Domains();
	for (int n=0; n<NDOM; ++n)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(n));

		// don't forget to multiply with the density
		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
		double dens = pme->Density();

		int NEL = dom.Elements();
		for (int i=0; i<NEL; ++i)
		{
			FESolidElement& el = dom.Element(i);

			// get the element force vector and initialize it to zero
			int ndof = 3*el.Nodes();
			fe.assign(ndof, 0);

			// apply body forces
			ElementWarpForce(dom, el, fe, dens);

			// get the element's LM vector
			dom.UnpackLM(el, lm);

			// assemble element 'fe'-vector into global R vector
			psolver->AssembleResidual(el.m_node, lm, fe, R);
		}
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::ElementWarpForce(FESolidDomain& dom, FESolidElement& el, vector<double>& fe, double dens)
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	vec3d f;

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
		FEMaterialPoint& mp = *el.m_State[n];
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		pt.r0 = el.Evaluate(r0, n);
		pt.rt = el.Evaluate(rt, n);

		detJ = dom.detJ0(el, n)*gw[n];

		// get the force
		f = m_wrp.force(mp);

		H = el.H(n);

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] -= H[i]*dens*f.x*detJ;
			fe[3*i+1] -= H[i]*dens*f.y*detJ;
			fe[3*i+2] -= H[i]*dens*f.z*detJ;
		}						
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::StiffnessMatrix(FENLSolver* psolver)
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	// element stiffness matrix
	matrix ke;
	vector<int> lm;

	// loop over all domains
	int NDOM = mesh.Domains();
	for (int n=0; n<NDOM; ++n)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(n));

		// don't forget to multiply with the density
		FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(dom.GetMaterial());
		double dens = pme->Density();

		// repeat over all solid elements
		int NE = dom.Elements();
		for (int iel=0; iel<NE; ++iel)
		{
			FESolidElement& el = dom.Element(iel);

			// create the element's stiffness matrix
			int ndof = 3*el.Nodes();
			ke.resize(ndof, ndof);
			ke.zero();

			// calculate inertial stiffness
			ElementWarpStiffness(dom, el, ke, dens);

			// get the element's LM vector
			dom.UnpackLM(el, lm);

			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::ElementWarpStiffness(FESolidDomain& dom, FESolidElement& el, matrix& ke, double dens)
{
	int neln = el.Nodes();
	int ndof = ke.columns()/neln;

	// jacobian
	double detJ;
	double *H;
	double* gw = el.GaussWeights();
	mat3ds K;

	// loop over integration points
	int nint = el.GaussPoints();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.m_State[n];
		detJ = dom.detJ0(el, n)*gw[n];

		// get the stiffness
		K = m_wrp.stiffness(mp);

		H = el.H(n);

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
bool FEWarpConstraint::Augment(int naug)
{
	return true;
}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Update()
{

}

//-----------------------------------------------------------------------------
void FEWarpConstraint::Serialize(DumpFile& ar)
{

}
