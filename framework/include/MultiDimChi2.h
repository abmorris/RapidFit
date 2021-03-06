#ifndef MULTIDIMCHI2_H
#define MULTIDIMCHI2_H
// Standard C++
#include <string>
#include <vector>
#include <memory>
// ROOT
#include "THn.h"
// RapidFit
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "PDFWithData.h"

class MultiDimChi2
{
	public:
		MultiDimChi2(const std::vector<PDFWithData*>& _allObjects, std::vector<std::pair<int, std::string>> wantedObservables);
		struct Result
		{
			Result(double _chi2, int _ndof) : chi2(_chi2), ndof(_ndof) {}
			double chi2;
			int ndof;
		};
		std::vector<MultiDimChi2::Result> PerformMuiltDimTest(bool poisson) const;
		static double CalcChi2(const std::vector<double>& expected_events, const std::vector<double>& observed_events, const bool poisson);
		static double CalculateExpected(IPDF& thisPDF, PhaseSpaceBoundary& fullPhaseSpace, const std::vector<DataPoint*>& combinations, const IDataSet& thisDataSet, const std::vector<ObservableRef>& theseObservables, const std::vector<std::pair<double, double>>& boundaries);
	private:
		// Stuff to store locally
		vector<ObservableRef> Observables;
		std::vector<std::unique_ptr<THnD>> BinnedData; // THnD is neither copyable nor moveable so I have to resort to this nonsense. Thanks ROOT!
		std::vector<PDFWithData*> allObjects;
		// Helper functions
		std::vector<std::pair<double, double>> GetBinBoundaries(const THnD& DataHist, const std::vector<int>& indices) const;
		std::vector<int> GetIndices(unsigned binNum, const THnD& DataHist) const;
};

#endif

