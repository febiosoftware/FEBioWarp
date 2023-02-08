#include "stdafx.h"
#include "FEWarpLog.h"
#include "FEWarpImageConstraint.h"
#include <FECore/FEModel.h>

FELogWarp::FELogWarp(FEModel* fem) : FELogElemData(fem) 
{
	// find the warping constraint
	m_wrp = nullptr;
	for (int i = 0; i < fem->NonlinearConstraints(); ++i)
	{
		m_wrp = dynamic_cast<FEWarpImageConstraint*>(fem->NonlinearConstraint(i));
		if (m_wrp) break;
	}
}

double FELogWarpTemplate::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	double T = 0.0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		T += tmap.value(r0);
	}
	T /= (double)el.Nodes();
	return T;
}

double FELogWarpTarget::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& smap = m_wrp->GetTargetMap();
	double S = 0.0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		S += smap.value(rt);
	}
	S /= (double)el.Nodes();
	return S;
}

double FELogWarpEnergy::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	ImageMap& smap = m_wrp->GetTargetMap();
	double e = 0.0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		double T = tmap.value(r0);
		double S = smap.value(rt);
		e += (0.5 * (T - S) * (T - S));
	}
	e /= (double)el.Nodes();
	return e;
}

vec3d FELogWarpForce_::force(FEElement& el)
{
	if (m_wrp == nullptr) return vec3d(0,0,0);
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	ImageMap& smap = m_wrp->GetTargetMap();
	vec3d fw(0, 0, 0);
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		double T = tmap.value(r0);
		double S = smap.value(rt);
		vec3d G = smap.gradient(rt);
		fw += G * ((T - S));
	}
	fw /= (double)el.Nodes();
	return fw;
}

double FELogWarpForceX::value(FEElement& el)
{
	vec3d fw = force(el);
	return fw.x;
}

double FELogWarpForceY::value(FEElement& el)
{
	vec3d fw = force(el);
	return fw.y;
}

double FELogWarpForceZ::value(FEElement& el)
{
	vec3d fw = force(el);
	return fw.z;
}
