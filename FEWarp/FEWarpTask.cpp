#include "stdafx.h"
#include "FEWarpTask.h"
#include <FECore/FEModel.h>
#include <FECore/XMLReader.h>
#include "FEWarpConstraint.h"
#include "FEWarpImageConstraint.h"
#include <FECore/xmltool.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FELoadCurve.h>
#include <FECore/FEPlotDataStore.h>
#include <iostream>
#include "FEWarpPlot.h"

FEWarpTask::FEWarpTask(FEModel* fem) : FECoreTask(fem)
{

}

FEWarpSingleImageConstraint* parseFile(FEModel* fem, XMLTag& tag)
{
	FEWarpSingleImageConstraint* pwc = fecore_alloc(FEWarpSingleImageConstraint, fem);

	++tag;
	while (!tag.isend())
	{
		if (tag == "constraint")
		{
			pwc->BuildParamList();
			if (fexml::readParameterList(tag, pwc) == false)
			{
				throw XMLReader::InvalidTag(tag);
			}

			fem->AddNonlinearConstraint(pwc);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}

	return pwc;
}

bool FEWarpTask::Init(const char* szfile)
{
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;

	XMLReader xml;
	if (xml.Open(szfile) == false) return false;
	XMLTag tag;
	if (!xml.FindTag("fewarp_spec", tag))
	{
		std::cerr << "This is not a valid FEWarp specification file." << std::endl;
		return false;
	}

	FEWarpSingleImageConstraint* pwc = nullptr;

	try {
		pwc = parseFile(fem, tag);
	}
	catch (...)
	{
		std::cerr << "Something went wrong.";
		return false;
	}

	xml.Close();

	// Add load controllers for the blur and penalty parameters
	int nlc = fem->LoadControllers();
	FELoadCurve* plcBlur = fecore_alloc(FELoadCurve, fem);
	plcBlur->SetID(nlc++);
	plcBlur->Clear();
	plcBlur->Add(0, 1);
	plcBlur->Add(1, 0);
	fem->AddLoadController(plcBlur);
	FEParam* blurParam = pwc->GetParameter("blur"); assert(blurParam != nullptr);
	fem->AttachLoadController(blurParam, plcBlur);

	FELoadCurve* plcPenalty = fecore_alloc(FELoadCurve, fem);
	plcPenalty->SetID(nlc++);
	plcPenalty->Clear();
	plcPenalty->Add(0, 0);
	plcPenalty->Add(1, 1);
	fem->AddLoadController(plcPenalty);
	FEParam* penaltyParam = pwc->GetParameter("penalty"); assert(penaltyParam != nullptr);
	fem->AttachLoadController(penaltyParam, plcPenalty);

	// Add the warping plot variables
	vector<int> dummy;
	FEPlotDataStore& plotData = fem->GetPlotDataStore();
	plotData.AddPlotVariable("warp-template", dummy);
	plotData.AddPlotVariable("warp-target", dummy);
	plotData.AddPlotVariable("warp-energy", dummy);

	return fem->Init();
}

bool FEWarpTask::Run()
{
	FEModel* fem = GetFEModel();
	if (fem == nullptr) return false;
	return fem->Solve();
}
