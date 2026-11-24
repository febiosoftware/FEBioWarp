#include "stdafx.h"
#include "FEWarpPlot.h"
#include "FEWarpConstraint.h"
#include "FEWarpVolumeConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FEDataStream.h>
#include <FECore/FEAnalysis.h>
#include <cmath>
#include <vector>

static FEWarpVolumeConstraint* FindWarpImageConstraint(FEModel& fem)
{
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		auto* pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc != nullptr) return pc;
	}

	return nullptr;
}

static bool GetImageMeanAndStd(ImageMap& map, double& mean, double& sigma)
{
	Image& image = map.GetImage();

	const int nx = image.width();
	const int ny = image.height();
	const int nz = image.depth();
	const int n = nx * ny * nz;

	if (n <= 0) return false;

	mean = 0.0;
	for (int k = 0; k < nz; ++k)
		for (int j = 0; j < ny; ++j)
			for (int i = 0; i < nx; ++i)
				mean += image.value(i, j, k);
	mean /= n;

	double sumSquaredDifference = 0.0;
	for (int k = 0; k < nz; ++k)
		for (int j = 0; j < ny; ++j)
			for (int i = 0; i < nx; ++i)
			{
				const double difference = image.value(i, j, k) - mean;
				sumSquaredDifference += difference * difference;
			}

	sigma = sqrt(sumSquaredDifference / n);
	return (sigma > 1e-12);
}

static bool SaveLocalNCC(FEMesh& m, ImageMap& tmap,	ImageMap& smap,	FEDataStream& s)
{
	constexpr int radius = 2;
	constexpr int minimumSamples = 50;
	constexpr double eps = 1e-12;

	const int nodeCount = m.Nodes();

	for (int i = 0; i < nodeCount; ++i)
	{
		const vec3d r0 = m.Node(i).m_r0;
		const vec3d rt = m.Node(i).m_rt;

		std::vector<vec3d> templatePoints;
		std::vector<vec3d> targetPoints;

		for (int dz = -radius; dz <= radius; ++dz)
			for (int dy = -radius; dy <= radius; ++dy)
				for (int dx = -radius; dx <= radius; ++dx)
				{
					const vec3d templatePoint = r0 + vec3d(dx * tmap.dx(), dy * tmap.dy(), dz * tmap.dz());
					const vec3d targetPoint = rt + vec3d(dx * smap.dx(), dy * smap.dy(), dz * smap.dz());

					if (tmap.valid(templatePoint) && smap.valid(targetPoint))
					{
						templatePoints.push_back(templatePoint);
						targetPoints.push_back(targetPoint);
					}
				}

		const int sampleCount = static_cast<int>(templatePoints.size());

		if (sampleCount < minimumSamples)
		{
			s << 0.0;
			continue;
		}

		double templateMean = 0.0;
		double targetMean = 0.0;

		for (int j = 0; j < sampleCount; ++j)
		{
			templateMean += tmap.value(templatePoints[j]);
			targetMean += smap.value(targetPoints[j]);
		}

		templateMean /= sampleCount;
		targetMean /= sampleCount;

		double numerator = 0.0;
		double templateNorm2 = 0.0;
		double targetNorm2 = 0.0;

		for (int j = 0; j < sampleCount; ++j)
		{
			const double templateDifference = tmap.value(templatePoints[j]) - templateMean;

			const double targetDifference = smap.value(targetPoints[j]) - targetMean;

			numerator += templateDifference * targetDifference;
			templateNorm2 += templateDifference * templateDifference;
			targetNorm2 += targetDifference * targetDifference;
		}

		const double denominator = std::sqrt(templateNorm2 * targetNorm2);

		s << (denominator > eps ? numerator / denominator : 0.0);
	}

	return true;
}

bool FEPlotTemplate::Save(FEMesh &m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);

		FEWarpVolumeConstraint* pci = dynamic_cast<FEWarpVolumeConstraint*>(pc);
		if (pci && pci->IsActive()) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs && pcs->IsActive()) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}

bool FEPlotTemplate::SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint* pc, FEDataStream& s)
{
	// get the template image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		if (tmap.valid(r0))
			s << tmap.value(r0);
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotTemplate::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s)
{
	FEWarpSurface* ps = pc->GetTemplate();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		s << ps->value(m.Node(i).m_r0);
	}
	return true;
}

bool FEPlotTarget::Save(FEMesh &m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* pc = fem.NonlinearConstraint(i);

		FEWarpVolumeConstraint* pci = dynamic_cast<FEWarpVolumeConstraint*>(pc);
		if (pci && pci->IsActive()) { return SaveWarpImage(m, pci, s); }

		FEWarpSurfaceConstraint* pcs = dynamic_cast<FEWarpSurfaceConstraint*>(pc);
		if (pcs && pcs->IsActive()) { return SaveWarpMesh(m, pcs, s); }
	}
	return true;
}


bool FEPlotTarget::SaveWarpImage(FEMesh& m, FEWarpVolumeConstraint* pc, FEDataStream& s)
{
	// get the target image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d rt = m.Node(i).m_rt;
		if (smap.valid(rt))
			s << smap.value(rt);
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotTarget::SaveWarpMesh(FEMesh& m, FEWarpSurfaceConstraint* pc, FEDataStream& s)
{
	FEWarpSurface* ps = pc->GetTarget();
	ps->Update();
	int N = m.Nodes();
	for (int i=0; i<N; ++i)
		s << ps->value(m.Node(i).m_rt);
	return true;
}

bool FEPlotEnergy::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i=0; i<N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);
			s << (0.5 * (T - S) * (T - S));
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotForce::Save(FEMesh &m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i=0; i<fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	double lam = ppc->GetPenalty();
	int N = m.Nodes();
	for (int i=0; i<N; ++i) 
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);
			vec3d G = smap.gradient(rt);
			vec3d fw = - G * ((T - S)) * lam;

			s << fw;
		}
		else 
			s << vec3d(0.0);
	}
	return true;
}

bool FEPlotNCCForce::Save(FEMesh& m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();

	FEWarpVolumeConstraint* warp = nullptr;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		warp = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (warp && warp->IsActive()) break;
	}

	if (warp == nullptr) return false;

	ImageMap& tmap = warp->GetTemplateMap();
	ImageMap& smap = warp->GetTargetMap();

	double tMean = 0.0;
	double tSigma = 0.0;
	double sMean = 0.0;
	double sSigma = 0.0;

	const bool validTemplate = GetImageMeanAndStd(tmap, tMean, tSigma);
	const bool validTarget = GetImageMeanAndStd(smap, sMean, sSigma);

	if ((validTemplate == false) || (validTarget == false))
	{
		for (int i = 0; i < m.Nodes(); ++i)
			s << vec3d(0.0);
		return true;
	}

	const double penalty = warp->GetPenalty();
	const int nodes = m.Nodes();

	for (int i = 0; i < nodes; ++i)
	{
		const vec3d r0 = m.Node(i).m_r0;
		const vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			const double tHat = (tmap.value(r0) - tMean) / tSigma;
			const double sHat = (smap.value(rt) - sMean) / sSigma;

			const vec3d force = -smap.gradient(rt) * (penalty / sSigma) * (tHat - sHat);

			s << force;
		}
		else
			s << vec3d(0.0);
	}

	return true;
}

bool FEPlotDiff::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap.valid(r0) && smap.valid(rt))
		{
			double T = tmap.value(r0);
			double S = smap.value(rt);

			s << T - S;
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotRawDiff::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		if (pc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap0 = pc->GetRawTemplateMap();
	ImageMap& smap0 = pc->GetRawTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;
		vec3d rt = m.Node(i).m_rt;

		if (tmap0.valid(r0) && smap0.valid(rt))
		{
			double T = tmap0.value(r0);
			double S = smap0.value(rt);

			s << T - S;
		}
		else
			s << 0.0;
	}
	return true;
}

bool FEPlotNCC::Save(FEMesh& m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();

	FEWarpVolumeConstraint* pc = FindWarpImageConstraint(fem);

	if (pc == nullptr) return false;

	return SaveLocalNCC(m, pc->GetTemplateMap(), pc->GetTargetMap(), s);
}

bool FEPlotRawNCC::Save(FEMesh& m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();

	FEWarpVolumeConstraint* pc = FindWarpImageConstraint(fem);

	if (pc == nullptr) return false;

	return SaveLocalNCC(m, pc->GetRawTemplateMap(),pc->GetRawTargetMap(), s);
}

static void CalculateImageStatistics(
	ImageMap& map,
	double& mean,
	double& sigma)
{
	Image& image = map.GetImage();

	const int nx = image.width();
	const int ny = image.height();
	const int nz = image.depth();
	const int voxelCount = nx * ny * nz;

	mean = 0.0;
	sigma = 0.0;

	for (int k = 0; k < nz; ++k)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int i = 0; i < nx; ++i)
			{
				mean += image.value(i, j, k);
			}
		}
	}

	mean /= voxelCount;

	for (int k = 0; k < nz; ++k)
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int i = 0; i < nx; ++i)
			{
				const double difference =
					image.value(i, j, k) - mean;

				sigma += difference * difference;
			}
		}
	}

	sigma = std::sqrt(sigma / voxelCount);
}

bool FEPlotNCCEnergy::Save(FEMesh& m, FEDataStream& s)
{
	FEModel& fem = *GetFEModel();

	FEWarpVolumeConstraint* pc = FindWarpImageConstraint(fem);

	if (pc == nullptr) return false;

	ImageMap& tmap = pc->GetTemplateMap();
	ImageMap& smap = pc->GetTargetMap();

	double templateMean = 0.0;
	double templateSigma = 0.0;
	double targetMean = 0.0;
	double targetSigma = 0.0;

	CalculateImageStatistics(tmap, templateMean, templateSigma);
	CalculateImageStatistics(smap, targetMean, targetSigma);

	constexpr double eps = 1e-12;

	for (int i = 0; i < m.Nodes(); ++i)
	{
		const vec3d r0 = m.Node(i).m_r0;
		const vec3d rt = m.Node(i).m_rt;

		if (!tmap.valid(r0) || !smap.valid(rt))
		{
			s << 0.0;
			continue;
		}

		const double templateNormalized = (tmap.value(r0) - templateMean) / (templateSigma + eps);

		const double targetNormalized = (smap.value(rt) - targetMean) / (targetSigma + eps);
		const double difference = templateNormalized - targetNormalized;

		s << 0.5 * difference * difference;
	}

	return true;
}

bool FEPlotGradT::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& tmap = pc->GetTemplateMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d r0 = m.Node(i).m_r0;

		if (tmap.valid(r0))
		{
			vec3d g = tmap.gradient(r0);
			s << g;
		}
		else
			s << vec3d(0.0);
	}
	return true;
}

bool FEPlotGradS::Save(FEMesh& m, FEDataStream& s)
{
	// find the warping constraint
	FEModel& fem = *GetFEModel();
	FEWarpVolumeConstraint* pc = 0;
	FEWarpConstraint* ppc = 0;
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		pc = dynamic_cast<FEWarpVolumeConstraint*>(fem.NonlinearConstraint(i));
		ppc = dynamic_cast<FEWarpConstraint*>(fem.NonlinearConstraint(i));
		if (pc && ppc) break;
	}
	if (pc == 0) return false;

	// get the image map
	ImageMap& smap = pc->GetTargetMap();

	int N = m.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d rt = m.Node(i).m_rt;

		if (smap.valid(rt))
		{
			vec3d g = smap.gradient(rt);
			s << g;
		}
		else
			s << vec3d(0.0);
	}
	return true;
}