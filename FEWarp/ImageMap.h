#pragma once
#include "Image.h"
#include "FECore/vec3d.h"

class ImageMap
{
public:
	struct POINT
	{
		int		i, j, k;
		double	h[8];
	};

public:
	ImageMap(Image& img);
	~ImageMap(void);

	void SetRange(vec3d r0, vec3d r1);

	POINT map(vec3d p);

	double value(POINT& p);
	double value(vec3d r) { return value(map(r)); }

	double dx() { return (m_r1.x - m_r0.x)/(double) (m_img.width () - 1); }
	double dy() { return (m_r1.y - m_r0.y)/(double) (m_img.height() - 1); }
	double dz() { int nz = m_img.depth(); if (nz == 1) return 1.0; else return (m_r1.z - m_r0.z)/(double) (m_img.depth () - 1); }

protected:
	Image&	m_img;
	vec3d	m_r0;
	vec3d	m_r1;	
};
