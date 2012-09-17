#include "StdAfx.h"
#include "FEWarpConstraint.h"
#include <FECore/FEElasticMaterial.h>
#include <FECore/FEModel.h>
#include <iostream>
using namespace std;

//=============================================================================

FEWarpConstraint::FEWarpConstraint(FEModel* pfem) : FENLConstraint(pfem)
{
	m_blaugon = false;
	m_altol = 0.01;
	m_k = 0;
}

//-----------------------------------------------------------------------------
FEWarpConstraint::~FEWarpConstraint(void)
{

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
		f = wrpForce(mp);

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
		K = wrpStiffness(mp);

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

//=============================================================================
// FEWarpImageConstraint
//-----------------------------------------------------------------------------

Image ds[3];	// first derivate of trg
Image dds[6];	// second derivate of trg (assumed symmetric: xx, yy, zz, xy, yz, xz)

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
BEGIN_PARAMETER_LIST(FEWarpImageConstraint, FEWarpConstraint);
	ADD_PARAMETER(m_tmp    , FE_PARAM_IMAGE_3D, "template");
	ADD_PARAMETER(m_trg    , FE_PARAM_IMAGE_3D, "target"  );
	ADD_PARAMETER(m_k      , FE_PARAM_DOUBLE  , "penalty" );
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL    , "laugon"  );
	ADD_PARAMETER(m_altol  , FE_PARAM_DOUBLE  , "altol"   );
	ADD_PARAMETERV(m_r0    , FE_PARAM_DOUBLEV, 3, "range_min");
	ADD_PARAMETERV(m_r1    , FE_PARAM_DOUBLEV, 3, "range_max");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEWarpImageConstraint::FEWarpImageConstraint(FEModel* pfem) : FEWarpConstraint(pfem), m_tmap(m_tmp), m_smap(m_trg)
{

}

//-----------------------------------------------------------------------------
void FEWarpImageConstraint::Init()
{
	int nx = m_tmp.width ();
	int ny = m_tmp.height();
	int nz = m_tmp.depth ();

	// allocate storage for image data
	ds[0].Create(nx, ny, nz);
	ds[1].Create(nx, ny, nz);
	ds[2].Create(nx, ny, nz);

	dds[0].Create(nx, ny, nz);
	dds[1].Create(nx, ny, nz);
	dds[2].Create(nx, ny, nz);
	dds[3].Create(nx, ny, nz);
	dds[4].Create(nx, ny, nz);
	dds[5].Create(nx, ny, nz);

	vec3d r0(m_r0[0], m_r0[1], m_r0[2]);
	vec3d r1(m_r1[0], m_r1[1], m_r1[2]);

	m_tmap.SetRange(r0, r1);
	m_smap.SetRange(r0, r1);
	ds1map.SetRange(r0, r1);
	ds2map.SetRange(r0, r1);
	ds3map.SetRange(r0, r1);

	dds1map.SetRange(r0, r1);
	dds2map.SetRange(r0, r1);
	dds3map.SetRange(r0, r1);
	dds4map.SetRange(r0, r1);
	dds5map.SetRange(r0, r1);
	dds6map.SetRange(r0, r1);

	image_derive_x(m_trg, ds[0]);
	image_derive_y(m_trg, ds[1]);
	image_derive_z(m_trg, ds[2]);

	image_derive_x(ds[0], dds[0]);
	image_derive_y(ds[1], dds[1]);
	image_derive_z(ds[2], dds[2]);
	image_derive_y(ds[0], dds[3]);
	image_derive_z(ds[1], dds[4]);
	image_derive_z(ds[0], dds[5]);
}

//-----------------------------------------------------------------------------
vec3d FEWarpImageConstraint::wrpForce(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double T = m_tmap.value(pt.r0);

	ImageMap::POINT ps = m_smap.map(pt.rt);

	double S = m_smap.value(ps);

	double dx = m_smap.dx();
	double dy = m_smap.dy();
	double dz = m_smap.dz();

	double dSx = ds1map.value(ps)/dx;
	double dSy = ds2map.value(ps)/dy;
	double dSz = ds3map.value(ps)/dz;

	vec3d Fw = (vec3d(dSx, dSy, dSz))*((S - T)*m_k);

	return Fw;
}

//-----------------------------------------------------------------------------
mat3ds FEWarpImageConstraint::wrpStiffness(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double T = m_tmap.value(pt.r0);

	ImageMap::POINT ps = m_smap.map(pt.rt);

	double S = m_smap.value(ps);

	double dx = m_smap.dx();
	double dy = m_smap.dy();
	double dz = m_smap.dz();

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
