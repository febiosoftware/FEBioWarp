#include "stdafx.h"
#include "FEWarpLog.h"
#include "FEWarpVolumeConstraint.h"
#include <FECore/FEModel.h>

FELogWarp::FELogWarp(FEModel* fem) : FELogElemData(fem) 
{
	// find the warping constraint
	m_wrp = nullptr;
	for (int i = 0; i < fem->NonlinearConstraints(); ++i)
	{
		m_wrp = dynamic_cast<FEWarpVolumeConstraint*>(fem->NonlinearConstraint(i));
		if (m_wrp) break;
	}
}

double FELogWarpTemplate::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	double T = 0.0;
	int n = 0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		T += tmap.valid(r0) ? tmap.value(r0) : 0.0;
		n += tmap.valid(r0) ? 1 : 0;
	}
	T /= (double)n;
	return T;
}

double FELogWarpTarget::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& smap = m_wrp->GetTargetMap();
	double S = 0.0;
	int n = 0;

	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		S += smap.valid(rt) ? smap.value(rt) : 0.0;
		n += smap.valid(rt) ? 1 : 0;
	}
	S /= (double)n;
	return S;
}

double FELogWarpEnergy::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	ImageMap& smap = m_wrp->GetTargetMap();
	double e = 0.0;
	int n = 0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		double T = tmap.valid(r0) ? tmap.value(r0) : 0.0;
		double S = tmap.valid(r0) ? smap.value(rt) : 0.0;
		e += (tmap.valid(r0) || smap.valid(rt)) ? (0.5 * (T - S) * (T - S)) : 0.0;
		n += (tmap.valid(r0) || smap.valid(rt)) ? 1 : 0;
	}
	e /= (double)n;
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
		double T = tmap.valid(r0) ? tmap.value(r0) : 0.0;
		double S = smap.valid(rt) ? smap.value(rt) : 0.0;
		vec3d G = smap.valid(rt) ? smap.gradient(rt) : vec3d(0.0);
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

double FELogWarpDiff::value(FEElement& el)
{
	if (m_wrp == nullptr) return 0.0;
	FEMesh& mesh = GetFEModel()->GetMesh();
	ImageMap& tmap = m_wrp->GetTemplateMap();
	ImageMap& smap = m_wrp->GetTargetMap();
	double e = 0.0;
	int n = 0;
	for (int i = 0; i < el.Nodes(); ++i)
	{
		vec3d r0 = mesh.Node(el.m_node[i]).m_r0;
		vec3d rt = mesh.Node(el.m_node[i]).m_rt;
		double T = tmap.valid(r0) ? tmap.value(r0) : 0.0;
		double S = smap.valid(rt) ? smap.value(rt) : 0.0;
		e += (tmap.valid(r0) || smap.valid(rt)) ? S - T : 0.0;
		n += (tmap.valid(r0) || smap.valid(rt)) ? 1 : 0.0;
	}
	e /= (double)el.Nodes();
	return e;
}
