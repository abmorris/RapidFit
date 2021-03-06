#include <algorithm>
#include <iomanip>
#include "TMath.h"
#include "MultiDimChi2.h"
MultiDimChi2::MultiDimChi2(const std::vector<PDFWithData*>& _allObjects, std::vector<std::pair<int, std::string>> wantedObservables) : allObjects(_allObjects)
{
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
			IConstraint* thisConstraint = thisBound->GetConstraint(ObservableRef(observable.second));
			if( !thisConstraint->IsDiscrete() )
			{
				x_min.push_back(thisConstraint->GetMinimum());
				x_max.push_back(thisConstraint->GetMaximum());
				x_bins.push_back(observable.first);
				std::cout << observable.first << " bins for observable " << observable.second << " with range [" << x_min.back() << "," << x_max.back() << "]" << "\n";
			}
			else
			{
				std::cerr << "Warning: discrete constraint on " << observable.second << " in this PhaseSpaceBoundary. Removing it." << std::endl;
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
		for(int i = 0; i < thisDataSet->GetDataNumber(); i++)
		{
			DataPoint* thisPoint = thisDataSet->GetDataPoint(i);
			std::vector<double> values;
			for(const auto& observable: wantedObservables)
				values.push_back(thisPoint->GetObservable(observable.second)->GetValue());
			double weight = weighted ? thisPoint->GetObservable(weightName)->GetValue() : 1.;
			BinnedData.back()->Fill(values.data(), weight);
		}
	}
	// Store the wanted observables locally
	for(const auto& observable: wantedObservables)
		Observables.push_back(ObservableRef(observable.second));
}

std::vector<MultiDimChi2::Result> MultiDimChi2::PerformMuiltDimTest(bool poisson) const
{
	std::vector<MultiDimChi2::Result> results;
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
		std::vector<double> observed_events, expected_events;
		unsigned nbins = 1;
		for(unsigned iobs = 0; iobs < Observables.size(); iobs++)
			nbins *= DataHist.GetAxis(iobs)->GetNbins();
		for(unsigned binNum = 0; binNum < nbins; binNum++)
		{
			std::vector<int> indices = GetIndices(binNum, DataHist);
			observed_events.push_back(DataHist.GetBinContent(indices.data()));
			expected_events.push_back(CalculateExpected(*thisPDF, *thisBound, thisBound->GetDiscreteCombinations(), *thisDataSet, Observables, GetBinBoundaries(DataHist, indices)));
		}
		// Calculate the chi2 by summing over all bins
		double TotalChi2 = CalcChi2(expected_events, observed_events, poisson);
		// Get the number of degrees of freedom
		double nDoF = nbins - thisPDF->GetPhysicsParameters()->GetAllFloatNames().size();
		// Print the result
		std::cout << DataHist.GetNdimensions() << "D chi2 result";
		if(allObjects.size() > 1) std::cout << " for fit #" << obj_counter;
		std::cout << "\n";
		std::cout << std::setw(12) << "chi2: " << TotalChi2 << "\n";
		std::cout << std::setw(12) << "nDoF: " << nDoF << "\n";
		std::cout << std::setw(12) << "chi2/nDoF: " << TotalChi2 / nDoF << "\n";
		std::cout << std::setw(12) << "p-value: " << TMath::Prob( TotalChi2, nDoF ) << "\n";
		std::cout << std::endl;
		results.push_back(Result(TotalChi2, nDoF));
	}
	return results;
}
std::vector<int> MultiDimChi2::GetIndices(unsigned binNum, const THnD& DataHist) const
{
	std::vector<int> indices;
	for(unsigned iobs = 0; iobs < Observables.size(); iobs++)
	{
		int nbins = DataHist.GetAxis(iobs)->GetNbins();
		indices.push_back((binNum % nbins) + 1);
		binNum /= nbins;
	}
	return indices;
}
std::vector<std::pair<double, double>> MultiDimChi2::GetBinBoundaries(const THnD& DataHist, const std::vector<int>& indices) const
{
	std::vector<std::pair<double, double>> binboundaries;
	for(unsigned iobs = 0; iobs < Observables.size(); iobs++)
	{
		double xMin = DataHist.GetAxis(iobs)->GetBinLowEdge(indices[iobs]);
		double xMax = DataHist.GetAxis(iobs)->GetBinUpEdge(indices[iobs]);
		binboundaries.push_back({xMin, xMax});
	}
	return binboundaries;
}
double MultiDimChi2::CalcChi2(const std::vector<double>& expected_events, const std::vector<double>& observed_events, const bool poisson)
{
	double Chi2Value = 0;
	for(unsigned i = 0; i < expected_events.size() && i < observed_events.size(); i++)
	{
		double thisChi2;
		if(poisson)
			thisChi2 = 2.0 * (expected_events[i] + observed_events[i] * (std::log(observed_events[i]/expected_events[i]) - 1));
		else
			thisChi2 = std::pow(expected_events[i] - observed_events[i],2) / expected_events[i];
		if( !std::isnan(thisChi2) && !std::isinf(thisChi2) )
			Chi2Value += thisChi2;
	}
	return Chi2Value;
}
double MultiDimChi2::CalculateExpected(IPDF& thisPDF, PhaseSpaceBoundary& fullPhaseSpace, const std::vector<DataPoint*>& combinations, const IDataSet& thisDataSet, const std::vector<ObservableRef>& theseObservables, const std::vector<std::pair<double, double>>& boundaries)
{
	if(theseObservables.size() != boundaries.size())
		std::cout << "Number of bin boundaries (" << boundaries.size() << ") does not match number of observables (" << theseObservables.size() << "). Proceeding to loop over the smallest number." << std::endl;
	double ExpectedEvents = 0;
	for(auto& combination: combinations)
	{
		// Get the integral over the full phase space
		double TotalIntegral = thisPDF.GetPDFIntegrator()->Integral(combination, &fullPhaseSpace);
		// Set the observable constraints to the bin boundaries
		PhaseSpaceBoundary binPhaseSpace(fullPhaseSpace);
		for(unsigned iobs = 0; iobs < theseObservables.size() && iobs < boundaries.size(); iobs++)
			binPhaseSpace.SetConstraint(theseObservables[iobs], boundaries[iobs].first, boundaries[iobs].second, "noUnits_Chi2");
		// Get the integral over the bin
		double BinIntegral = thisPDF.GetPDFIntegrator()->NumericallyIntegrateDataPoint(combination, &binPhaseSpace, thisPDF.GetDoNotIntegrateList());
		// Multiply by the size of the sample to get the expected events in this bin
		double SampleSize = thisDataSet.GetDataNumber(combination);
		ExpectedEvents += SampleSize*BinIntegral/TotalIntegral;
	}
	return ExpectedEvents;
}

