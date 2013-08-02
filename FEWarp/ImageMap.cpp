#include "stdafx.h"
#include "ImageMap.h"
#include "math.h"

ImageMap::ImageMap(Image& img) : m_img(img)
{
	m_r0 = vec3d(0,0,0);
	m_r1 = vec3d(0,0,0);
}

ImageMap::~ImageMap(void)
{
}

void ImageMap::SetRange(vec3d r0, vec3d r1)
{
	m_r0 = r0;
	m_r1 = r1;
}

ImageMap::POINT ImageMap::map(vec3d p)
{
	double x = (p.x - m_r0.x)/(m_r1.x - m_r0.x);
	double y = (p.y - m_r0.y)/(m_r1.y - m_r0.y);
	double z = (p.z - m_r0.z)/(m_r1.z - m_r0.z);

	int nx = m_img.width ();
	int ny = m_img.height();
	int nz = m_img.depth ();

	double dx = 1.0 / (double) (nx - 1);
	double dy = 1.0 / (double) (ny - 1);
	double dz = (nz>1?1.0 / (double) (nz - 1):1.0);

	int i = (int) floor(x*(nx - 1));
	int j = (int) floor(y*(ny - 1));
	int k = (int) floor(z*(nz - 1));

	if (i == nx-1) i--;
	if (j == ny-1) j--;
	if ((k == nz-1)&&(nz>1)) k--;

	double r = (x - i*dx)/dx;
	double s = (y - j*dy)/dy;
	double t = (z - k*dz)/dz;

	POINT pt;
	pt.i = i;
	pt.j = j;
	pt.k = k;

	if (nz == 1)
	{
		pt.h[0] = (1-r)*(1-s);
		pt.h[1] = r*(1-s);
		pt.h[2] = r*s;
		pt.h[3] = s*(1-r);
	}
	else
	{
		pt.h[0] = (1-r)*(1-s)*(1-t);
		pt.h[1] = r*(1-s)*(1-t);
		pt.h[2] = r*s*(1-t);
		pt.h[3] = s*(1-r)*(1-t);
		pt.h[4] = (1-r)*(1-s)*t;
		pt.h[5] = r*(1-s)*t;
		pt.h[6] = r*s*t;
		pt.h[7] = s*(1-r)*t;
	}

	return pt;
}

double ImageMap::value(POINT& p)
{
	int nx = m_img.width ();
	int ny = m_img.height();
	int nz = m_img.depth ();

	if (nz == 1)
	{
		if ((p.i<0) || (p.i >= nx-1)) return 0.0;
		if ((p.j<0) || (p.j >= ny-1)) return 0.0;

		double v[4];
		v[0] = m_img.value(p.i  , p.j  , 0);
		v[1] = m_img.value(p.i+1, p.j  , 0);
		v[2] = m_img.value(p.i+1, p.j+1, 0);
		v[3] = m_img.value(p.i  , p.j+1, 0);

		return (p.h[0]*v[0] + p.h[1]*v[1] + p.h[2]*v[2] + p.h[3]*v[3]);
	}
	else
	{
		if ((p.i<0) || (p.i >= nx-1)) return 0.0;
		if ((p.j<0) || (p.j >= ny-1)) return 0.0;
		if ((p.k<0) || (p.k >= nz-1)) return 0.0;

		double v[8];
		v[0] = m_img.value(p.i  , p.j  , p.k  );
		v[1] = m_img.value(p.i+1, p.j  , p.k  );
		v[2] = m_img.value(p.i+1, p.j+1, p.k  );
		v[3] = m_img.value(p.i  , p.j+1, p.k  );
		v[4] = m_img.value(p.i  , p.j  , p.k+1);
		v[5] = m_img.value(p.i+1, p.j  , p.k+1);
		v[6] = m_img.value(p.i+1, p.j+1, p.k+1);
		v[7] = m_img.value(p.i  , p.j+1, p.k+1);

		return (p.h[0]*v[0] + p.h[1]*v[1] + p.h[2]*v[2] + p.h[3]*v[3] + p.h[4]*v[4] + p.h[5]*v[5] + p.h[6]*v[6] + p.h[7]*v[7]);
	}
}
