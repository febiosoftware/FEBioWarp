#include "stdafx.h"
#include "FEWarpImageConstraint.h"
#include <FEBioMech/FEElasticMaterial.h>

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
bool FEWarpImageConstraint::Init()
{
	FEWarpConstraint::Init();

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

	return true;
}

//-----------------------------------------------------------------------------
vec3d FEWarpImageConstraint::wrpForce(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double T = m_tmap.value(pt.m_r0);

	ImageMap::POINT ps = m_smap.map(pt.m_rt);

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

	double T = m_tmap.value(pt.m_r0);

	ImageMap::POINT ps = m_smap.map(pt.m_rt);

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
