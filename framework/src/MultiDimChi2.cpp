#include <algorithm>
#include <iomanip>
#include "TMath.h"
#include "MultiDimChi2.h"
MultiDimChi2::MultiDimChi2(const std::vector<PDFWithData*>& _allObjects, vector<std::string> wantedObservables) : allObjects(_allObjects)
{
	std::cout << "MultiDimChi2 will fill " << allObjects.size() << " " << wantedObservables.size() << "D histogram(s)\n";
	// Loop through each fit and bin each dataset
	int obj_counter = 0; // Just to give each histogram a unique name
	for(const auto PDFAndData: allObjects)
	{
		obj_counter++;
		IDataSet* thisDataSet = PDFAndData->GetDataSet();
		PhaseSpaceBoundary* thisBound = thisDataSet->GetBoundary();
		// Get the boundaries of each observable and construct the binning scheme
		std::vector<double> x_min, x_max;
		std::vector<int> x_bins;
		for(const auto& observable: wantedObservables)
		{
			IConstraint* thisConstraint = thisBound->GetConstraint(ObservableRef(observable));
			if( !thisConstraint->IsDiscrete() )
			{
				x_min.push_back(thisConstraint->GetMinimum());
				x_max.push_back(thisConstraint->GetMaximum());
				x_bins.push_back(5); // TODO: read from config
				std::cout << x_bins.back() << " bins for observable " << observable << "\n";
			}
			else
			{
				std::cerr << "Warning: discrete constraint on " << observable << " in this PhaseSpaceBoundary. Removing it." << std::endl;
				wantedObservables.erase(std::remove(wantedObservables.begin(), wantedObservables.end(), observable), wantedObservables.end());
			}
		}
		// Construct the histogram and add it to the container
		std::string histname = "BinnedData_"+std::to_string(obj_counter);
		BinnedData.push_back(std::unique_ptr<THnD>(new THnD(histname.c_str(), "", x_bins.size(), x_bins.data(), x_min.data(), x_max.data())));
		// Get the weight name, if it exists
		bool weighted = thisDataSet->GetWeightsWereUsed();
		ObservableRef weightName;
		if(weighted) weightName = ObservableRef(thisDataSet->GetWeightName());
		// Fill the histogram
		std::cout << "Filling histogram " << obj_counter << " with " << thisDataSet->GetDataNumber() << " points\n";
		for(int i = 0; i < thisDataSet->GetDataNumber(); i++)
		{
			DataPoint* thisPoint = thisDataSet->GetDataPoint(i);
			std::vector<double> values;
			for(const auto& observable: wantedObservables)
				values.push_back(thisPoint->GetObservable(observable)->GetValue());
			double weight = weighted ? thisPoint->GetObservable(weightName)->GetValue() : 1.;
			BinnedData.back()->Fill(values.data(), weight);
		}
	}
	std::cout << "Filled histograms\n";
	// Store the wanted observables locally
	for(const auto& observable: wantedObservables)
		Observables.push_back(ObservableRef(observable));
	std::cout << "MultiDimChi2 successfully constructed" << std::endl;
}

void MultiDimChi2::PerformMuiltDimTest() const
{
	int obj_counter = 0;
	for(const auto PDFAndData: allObjects)
	{
		// Retrieve the binned data hist
		const THnD& DataHist = *BinnedData[obj_counter];
		obj_counter++;
		// Get the PDF and the boundary
		IPDF* thisPDF = PDFAndData->GetPDF();
		IDataSet* thisDataSet = PDFAndData->GetDataSet();
		PhaseSpaceBoundary* thisBound = thisDataSet->GetBoundary();
		// Get the observed and expected number of events per bin
		std::vector<double> observed_events, error_events, expected_events;
		std::cout << "Looping over " << DataHist.GetNbins() << " bins in histogram " << obj_counter << "\n";
		for(unsigned binNum = 1; binNum <= DataHist.GetNbins(); binNum++)
		{
			observed_events.push_back(DataHist.GetBinContent(binNum));
			error_events.push_back(DataHist.GetBinError(binNum));
			expected_events.push_back(CalculateExpected(*thisPDF, *thisBound, *thisDataSet, DataHist, binNum));
		}
		// Calculate the chi2 by summing over all bins
		double TotalChi2 = CalcChi2(expected_events, observed_events, error_events);
		// Get the number of degrees of freedom
		double nDoF = DataHist.GetNbins() - thisPDF->GetPhysicsParameters()->GetAllFloatNames().size();
		// Print the result
		std::cout << DataHist.GetNdimensions() << "D chi2 result";
		if(allObjects.size() > 1) std::cout << " for fit #" << obj_counter;
		std::cout << "\n";
		std::cout << std::setw(12) << "chi2: " << TotalChi2 << "\n";
		std::cout << std::setw(12) << "nDoF: " << nDoF << "\n";
		std::cout << std::setw(12) << "chi2/nDoF: " << TotalChi2 / nDoF << "\n";
		std::cout << std::setw(12) << "p-value: " << TMath::Prob( TotalChi2, nDoF ) << "\n";
		std::cout << std::endl;
	}
	return;
}
double MultiDimChi2::CalcChi2(const std::vector<double>& expected_events, const std::vector<double>& observed_events, const std::vector<double>& errors) const
{
	(void)errors; // XXX why?
	double Chi2Value = 0;
	for(unsigned i=0; i < expected_events.size(); i++)
	{
		double thisChi2 = expected_events[i] - observed_events[i] + observed_events[i] * log( observed_events[i] / expected_events[i] );
		if( !std::isnan(thisChi2) && !std::isinf(thisChi2) )
			Chi2Value += thisChi2;
	}
	return 2.*Chi2Value;
}
double MultiDimChi2::CalculateExpected(IPDF& thisPDF, PhaseSpaceBoundary& fullPhaseSpace, const IDataSet& thisDataSet, const THnD& DataHist, int bindex) const
{
	double ExpectedEvents = 0;
	for(auto& combination: fullPhaseSpace.GetDiscreteCombinations())
	{
		// Get the integral over the full phase space
		double TotalIntegral = thisPDF.GetPDFIntegrator()->Integral(combination, &fullPhaseSpace);
		// Set the observable constraints to the bin boundaries
		PhaseSpaceBoundary binPhaseSpace(fullPhaseSpace);
		for(unsigned iobs = 0; iobs < Observables.size(); iobs++)
		{
			double newMin = DataHist.GetAxis(iobs)->GetBinLowEdge(bindex);
			double newMax = DataHist.GetAxis(iobs)->GetBinUpEdge(bindex);
			binPhaseSpace.SetConstraint(Observables[iobs], newMin, newMax, "noUnits_Chi2");
		}
		// Get the integral over the bin
		double BinIntegral = thisPDF.GetPDFIntegrator()->NumericallyIntegrateDataPoint(combination, &binPhaseSpace, thisPDF.GetDoNotIntegrateList());
		// Multiply by the size of the sample to get the expected events in this bin
		double SampleSize = thisDataSet.GetDataNumber(combination);
		ExpectedEvents += SampleSize*BinIntegral/TotalIntegral;
	}
	return ExpectedEvents;
}

