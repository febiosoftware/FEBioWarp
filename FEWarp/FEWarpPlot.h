#pragma once

#include <FECore/FEPlotData.h>
#include "FEWarpVolumeConstraint.h"
#include "FEWarpSurfaceConstraint.h"

class FEPlotTemplate : public FEPlotNodeData
{
public:
	FEPlotTemplate(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);

protected:
	bool SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint*   pc, FEDataStream& s);
	bool SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s);
};

class FEPlotTarget : public FEPlotNodeData
{
public:
	FEPlotTarget(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);

protected:
	bool SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint*   pc, FEDataStream& s);
	bool SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s);
};

class FEPlotEnergy : public FEPlotNodeData
{
public:
	FEPlotEnergy(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotNCCEnergy : public FEPlotNodeData
{
public:
	FEPlotNCCEnergy(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotForce : public FEPlotNodeData
{
public:
	FEPlotForce(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotNCCForce : public FEPlotNodeData
{
public:
	FEPlotNCCForce(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotNCCForce2 : public FEPlotNodeData
{
public:
	FEPlotNCCForce2(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotDiff : public FEPlotNodeData
{
public:
	FEPlotDiff(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotRawDiff : public FEPlotNodeData
{
public:
	FEPlotRawDiff(FEModel* pfem) : FEPlotNodeData(pfem, PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotGradT : public FEPlotNodeData
{
public:
	FEPlotGradT(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};

class FEPlotGradS : public FEPlotNodeData
{
public:
	FEPlotGradS(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, FEDataStream& s);
};