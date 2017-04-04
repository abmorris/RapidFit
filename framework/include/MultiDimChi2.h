
#ifndef MULTIDIMCHI2_H
#define MULTIDIMCHI2_H

#include "THn.h"
#include <vector>
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"
#include "PDFWithData.h"
#include <string>
#include "IPDF.h"

using std::vector;
using std::string;

struct ThisObsBinning
{
	string ObservableName;
	vector<double> binCenters;
	double thisMin;
	double thisMax;
	double thisStepSize;
	unsigned int theseBins;
};

class MultiDimChi2
{
	public:
		MultiDimChi2( vector<PDFWithData*> allObjects, PhaseSpaceBoundary* thisBound, vector<string> wantedObservables );
		void PerformMuiltDimTest();
	private:
		// Initialisation stuff
		void ConstructAllCoordinates();
		void AddCoordinates( unsigned int thisDim );
		void ConstructBinCenters();
		void ConstructBoundaries( PhaseSpaceBoundary* totalPhaseSpace, vector<string> );
		void ConstructInternalHisto( vector<string> wantedObservables, PhaseSpaceBoundary* thisBound );
		void populateAllObjects( vector<PDFWithData*> allObjects );
		// Calculation helpers
		double CalcChi2( vector<double> expected_events, vector<double> observed_events, vector<double> );
		double CalculateTotalExpected( vector<double> thisBinCenter );
		double CalculateRange( PhaseSpaceBoundary* thisBound );
		void ConstructIntegralsRatios( vector<string> wantedObservables );
		double CorrectYield( IDataSet* thisSet, DataPoint* thisPoint );
		double PDF2DataNormalisation( unsigned int PDFNum, const unsigned int combinationIndex, DataPoint* thisDataPoint );
		// Member variables
		vector<double> x_min;
		vector<double> x_max;
		vector<ObservableRef> goodObservables;
		vector<int> x_bins;
		std::unique_ptr<THnD> internalHisto; // Because THnD isn't copyable for no apparent reason
		vector<ThisObsBinning> theseDimensions;
		unsigned int nDim;
		vector<vector<double>> allBinCenters;
		vector<IPDF*> allPDFs;
		vector<IDataSet*> allDataSets;
		vector<PDFWithData*> allPDFData;
		vector<PhaseSpaceBoundary*> allBoundaries;
		unsigned int data_binning;
		vector<vector<double> > ratioOfIntegrals;
		vector<vector<double> > combinationIntegrals;
		vector<double> weightNorms;
};

#endif


