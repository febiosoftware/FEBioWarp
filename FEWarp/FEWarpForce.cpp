#include "StdAfx.h"
#include "FEWarpForce.h"
#include "Image.h"
#include "ImageMap.h"
#include <FECore/FEElasticMaterial.h>
#include <iostream>
using namespace std;

BEGIN_PARAMETER_LIST(FEWarpForce, FEBodyForce);
	ADD_PARAMETER(m_sztmp, FE_PARAM_STRING, "template");
	ADD_PARAMETER(m_sztrg, FE_PARAM_STRING, "target"  );
	ADD_PARAMETER(m_k    , FE_PARAM_DOUBLE, "penalty" );
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
FEWarpForce::~FEWarpForce(void)
{
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
