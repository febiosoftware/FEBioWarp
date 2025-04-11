#pragma once
#include <FECore/ElementDataRecord.h>

class FEWarpImageConstraint;

class FELogWarp : public FELogElemData
{
public:
	FELogWarp(FEModel* fem);

protected:
	FEWarpImageConstraint* m_wrp;
};

class FELogWarpTemplate : public FELogWarp
{
public:
	FELogWarpTemplate(FEModel* fem) : FELogWarp(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpTarget : public FELogWarp
{
public:
	FELogWarpTarget(FEModel* fem) : FELogWarp(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpEnergy : public FELogWarp
{
public:
	FELogWarpEnergy(FEModel* fem) : FELogWarp(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpForce_ : public FELogWarp
{
public: 
	FELogWarpForce_(FEModel* fem) : FELogWarp(fem) {}
	vec3d force(FEElement& el);
};

class FELogWarpForceX : public FELogWarpForce_
{
public:
	FELogWarpForceX(FEModel* fem) : FELogWarpForce_(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpForceY : public FELogWarpForce_
{
public:
	FELogWarpForceY(FEModel* fem) : FELogWarpForce_(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpForceZ : public FELogWarpForce_
{
public:
	FELogWarpForceZ(FEModel* fem) : FELogWarpForce_(fem) {}
	double value(FEElement& el) override;
};

class FELogWarpDiff : public FELogWarp
{
public:
	FELogWarpDiff(FEModel* fem) : FELogWarp(fem) {}
	double value(FEElement& el) override;
};
