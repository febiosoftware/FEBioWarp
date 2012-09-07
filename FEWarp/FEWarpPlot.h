#pragma once

#include "FECore/FEPlotData.h"

class FEPlotTemplate : public FENodeData
{
public:
	FEPlotTemplate(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, vector<float>& a);
};

class FEPlotTarget : public FENodeData
{
public:
	FEPlotTarget(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, vector<float>& a);
};

class FEPlotEnergy : public FENodeData
{
public:
	FEPlotEnergy(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_NODE) {}
	virtual bool Save(FEMesh& m, vector<float>& a);
};

class FEPlotForce : public FENodeData
{
public:
	FEPlotForce(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_NODE) {}
	virtual bool Save(FEMesh& m, vector<float>& a);
};
