#pragma once

class Image
{
public:
	Image(void);
	~Image(void);

	void Create(int nx, int ny, int nz);
	bool Load(const char* szfile);

	int width () { return m_nx; }
	int height() { return m_ny; }
	int depth () { return m_nz; }

	float& value(int x, int y, int z) { return m_pf[(z*m_ny + y)*m_nx+x]; }

	void zero();

protected:
	float*	m_pf;
	int		m_nx, m_ny, m_nz;
};

void image_derive_x(Image& s, Image& d);
void image_derive_y(Image& s, Image& d);
void image_derive_z(Image& s, Image& d);
