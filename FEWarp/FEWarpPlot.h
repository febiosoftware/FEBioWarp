#pragma once

#include "FECore/FEPlotData.h"
#include "FEWarpImageConstraint.h"
#include "FEWarpSurfaceConstraint.h"

class FEPlotTemplate : public FENodeData
{
public:
	FEPlotTemplate(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE), m_pfem(pfem) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);

protected:
	bool SaveWarpImage(FEMesh& m, FEWarpImageConstraint*   pc, FEDataStream& s);
	bool SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s);

protected:
	FEModel*	m_pfem;
};

class FEPlotTarget : public FENodeData
{
public:
	FEPlotTarget(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE), m_pfem(pfem) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);

protected:
	bool SaveWarpImage(FEMesh& m, FEWarpImageConstraint*   pc, FEDataStream& s);
	bool SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s);

protected:
	FEModel*	m_pfem;
};

class FEPlotEnergy : public FENodeData
{
public:
	FEPlotEnergy(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE), m_pfem(pfem) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
protected:
	FEModel*	m_pfem;
};

class FEPlotForce : public FENodeData
{
public:
	FEPlotForce(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE), m_pfem(pfem) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
protected:
	FEModel*	m_pfem;
};
