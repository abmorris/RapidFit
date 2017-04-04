
#include "TH1D.h"
#include "TCanvas.h"
#include "THn.h"
#include "MultiDimChi2.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"
#include <vector>
#include "StringProcessing.h"
#include "IPDF.h"
#include <string>
#include <cmath>
#include "ClassLookUp.h"
#include "PDFWithData.h"
#include "RapidFitIntegrator.h"
#include "RapidFitIntegratorConfig.h"
#include "DebugClass.h"
#include "TMath.h"
MultiDimChi2::MultiDimChi2( vector<PDFWithData*> allObjects, PhaseSpaceBoundary* thisBound, vector<string> wantedObservables )
{
	//	Populate PDF, DataSets and PhaseSpaceBoundaries to be used in this analysis
	cout << "MultiDimChi2 Hello!" << endl << endl;
	cout << "initializing objects" << endl;
	this->populateAllObjects( allObjects );
	this->ConstructBoundaries( thisBound, wantedObservables );
	cout << "Deciding binning and constructing THnD" << endl;
	//	Construct the Histogram object which contains all of the Data of Interest
	this->ConstructInternalHisto( wantedObservables, thisBound );
	cout << "Calculating Bin Centers" << endl;
	//      Construct the Central coordinates of every bin in each axis, this will be the basis for generating the corrdinates to run the Chi2 test over
	this->ConstructBinCenters();
	cout << "Calculating the Numerical vs Analytical Ratios" << endl;
	//	Construct the Ratios of the Numerical vs analytical Integrals (condensing the differenced between the PDF numerator and denominator to be a single number but what else can we do...)
	this->ConstructIntegralsRatios( wantedObservables );
	cout << "Constructing all of the Coordinates to perform a Chi2 test over" << endl;
	//	Construct the full set of coordinates that the Chi2 test is to be run over
	this->ConstructAllCoordinates();
	cout << "Finished Initializing!" << endl << endl;
}
void MultiDimChi2::populateAllObjects( vector<PDFWithData*> allObjects )
{
	allPDFData = allObjects;	//	Copy Reference, don't take ownership!
	for( unsigned int i=0; i< allObjects.size(); ++i )
	{
		allDataSets.push_back( allObjects[i]->GetDataSet() );	//      Copy Reference, don't take ownership!
		allPDFs.push_back( ClassLookUp::CopyPDF( allObjects[i]->GetPDF() ) );	// 	More than likely going to cause something to change so lets make it on our copy, just incase
		RapidFitIntegratorConfig* thisConf = allPDFs.back()->GetPDFIntegrator()->GetIntegratorConfig();
		thisConf->useGSLIntegrator = true;
		thisConf->numThreads = 16; // XXX: why?
		allPDFs.back()->GetPDFIntegrator()->SetUpIntegrator( thisConf );
	}
}
void MultiDimChi2::ConstructIntegralsRatios( vector<string> wanted_params )
{
	(void) wanted_params;
	ratioOfIntegrals.clear();
	combinationIntegrals.clear();
	weightNorms.clear();
	for( unsigned int PDFNum=0; PDFNum< allPDFs.size(); ++PDFNum )
	{
		cout << "\tCalculating Ratios for PDF: " << PDFNum+1 << endl;
		IPDF* thisPDF = allPDFs[PDFNum];
		RapidFitIntegrator* pdfIntegrator = thisPDF->GetPDFIntegrator();
		IDataSet* thisData = allDataSets[PDFNum];
		vector<double> thisratioOfIntegrals;
		vector<double> thisCombination_integral;
		PhaseSpaceBoundary* thisBound = allBoundaries[PDFNum];
		thisPDF->ChangePhaseSpace( thisBound );
		vector<DataPoint*> allCombinations = thisBound->GetDiscreteCombinations();
		thisratioOfIntegrals.assign(allCombinations.size(),1.);
		thisCombination_integral.assign(allCombinations.size(),1.);
		if( thisratioOfIntegrals.size() == 0 || thisratioOfIntegrals.empty() ) thisratioOfIntegrals =vector<double> ( 1, 1. );
		ratioOfIntegrals.push_back( thisratioOfIntegrals );
		double thisWeight_norm = 0.;
		double weight_sum = thisData->GetSumWeights();
		thisWeight_norm = weight_sum / thisData->GetDataNumber(NULL);
		weightNorms.push_back( thisWeight_norm );
		cout << "Alpha: " << thisData->GetAlpha() << endl;
		pdfIntegrator->ForceTestStatus( true );
	}
}
void MultiDimChi2::ConstructInternalHisto( vector<string> wantedObservables, PhaseSpaceBoundary* thisBound )
{
	nDim=0;
	x_min.clear(); x_max.clear(); goodObservables.clear(); x_bins.clear();
	cout << "\tDecicing Binning" << endl;
	data_binning = 1;
	unsigned int x_binning=0;
	for( unsigned int i=0; i< wantedObservables.size(); ++i )
	{
		IConstraint* thisConstraint = thisBound->GetConstraint( wantedObservables[i] );
		if( !thisConstraint->IsDiscrete() )
		{
			++nDim;
			x_min.push_back( thisConstraint->GetMinimum() );
			x_max.push_back( thisConstraint->GetMaximum() );
			goodObservables.push_back( wantedObservables[i] );
			if( wantedObservables[i] == "time" ) x_binning = 50;
			else x_binning = 4;
			//x_binning = 2;
			x_bins.push_back( x_binning );
			data_binning *= x_binning;
			cout << wantedObservables[i] << "  " << x_min.back() << "<->" << x_max.back() << "  /  " << x_binning << endl;
			ThisObsBinning thisDimension;
			thisDimension.ObservableName = wantedObservables[i];
			thisDimension.thisMin = x_min.back();
			thisDimension.thisMax = x_max.back();
			thisDimension.theseBins = x_binning;
			thisDimension.thisStepSize = (thisDimension.thisMax-thisDimension.thisMin)/((double)thisDimension.theseBins);
			theseDimensions.push_back( thisDimension );
		}
	}
	cout << "\tConstructing THnD" << endl;
	cout << "\t\tinternal_Chi2nDim\t" << nDim << "\t" << x_bins[0] << "\t" << x_min[0] << "\t" << x_max[0] << endl;
	internalHisto = std::unique_ptr<THnD>(new THnD( "internal_Chi2nDim", "internal_Chi2nDim", (int)nDim, &(x_bins[0]), &(x_min[0]), &(x_max[0]) ));
	vector<double> thisValues( goodObservables.size(), 0. );
	double thisWeight = 0.;
	cout << "\tPopulating THnD" << endl;
	ObservableRef weightName;
	for( unsigned int i=0; i< allDataSets.size(); ++i )
	{
		weightName = ObservableRef( allDataSets[i]->GetWeightName() );
		for( unsigned int j=0; j< (unsigned)allDataSets[i]->GetDataNumber(); ++j )
		{
			DataPoint* thisPoint = allDataSets[i]->GetDataPoint( j );
			cout << "data:" << endl;
			for( unsigned int k=0; k< goodObservables.size(); ++k )
			{
				thisValues[k] = thisPoint->GetObservable( goodObservables[k] )->GetValue();
			}
			cout << "weights:" << endl;
			thisWeight = 1.;//thisPoint->GetObservable( weightName )->GetValue();
			internalHisto->Fill( &(thisValues[0]), thisWeight );
		}
	}
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		TH1D* h1 = internalHisto->Projection( i, "E" );
		TString name("c");
		name+=i;
		TCanvas c1( name, name );
		h1->Draw("");
		c1.Print("Output.pdf");
	}
}
void MultiDimChi2::ConstructBinCenters()
{
	double thisMin=0.;
	double thisMax=0.;
	int theseBins=0;
	double thisCenter=0.;
	double thisStepSize=0.;
	double numSteps=0.;
	for( unsigned int i=0; i< nDim; ++i )
	{
		vector<double> thisObservableBins;
		thisMin = x_min[i];
		thisMax = x_max[i];
		theseBins = x_bins[i];
		thisStepSize = (thisMax-thisMin)/((double)theseBins);
		for( unsigned int j=1; j<= (unsigned)theseBins; ++j )
		{
			numSteps=(double)j; numSteps-=0.5;
			thisCenter = thisMin+numSteps*thisStepSize;
			thisObservableBins.push_back( thisCenter );
		}
		theseDimensions[i].binCenters = thisObservableBins;
	}
}
void MultiDimChi2::PerformMuiltDimTest()
{
	cout << "MultiDimChi2: About to Perform Chi1 Calculation" << endl;
	PhaseSpaceBoundary* thisBoundary=NULL;
	(void) thisBoundary;
	vector<double> expected_events( allBinCenters.size(), 0. );
	vector<double> observed_events( allBinCenters.size(), 0. );
	vector<double> error_events( allBinCenters.size(), 0. );
	unsigned int histo_binNum=0;
	vector<double> thisBinCenter( (allBinCenters)[0].size(), 0. );
	cout << "Looping over: " << allBinCenters.size()+1 << " coordianates!" << endl;
	for( unsigned int binNum=0; binNum< allBinCenters.size(); ++binNum )
	{
		thisBinCenter = allBinCenters.at( binNum );
		histo_binNum = (unsigned)internalHisto->GetBin( &(thisBinCenter[0]) );
		cout << "Coordinate: " << binNum+1 << " of: " << allBinCenters.size() << endl;
		cout << "Determining Number of Observed Events in Bin: " << histo_binNum << endl;
		observed_events[binNum] = internalHisto->GetBinContent( histo_binNum );
		cout << "Calculating Number of Expected Events in Bin: " << histo_binNum << endl;
		expected_events[binNum] = this->CalculateTotalExpected( thisBinCenter );
		error_events[binNum] = internalHisto->GetBinError( histo_binNum );
		cout << "O: " << observed_events[binNum] << "  E: " << expected_events[binNum] << "  Err: " << error_events[binNum] << endl;
	}
	cout << "Finished Looping over all coordinates!" << endl;
	double TotalChi2 = this->CalcChi2( expected_events, observed_events, error_events );
	DebugClass::Dump2TTree( "file.root", expected_events, "TTree", "expected" );
	DebugClass::Dump2TTree( "file2.root", observed_events, "TTree", "observed" );
	DebugClass::Dump2TTree( "file3.root", error_events, "TTree", "error" );
	vector<double> pull; for( unsigned int i=0; i< expected_events.size(); ++i ) pull.push_back( (observed_events[i] - expected_events[i])/error_events[i] );
	DebugClass::Dump2TTree( "file4.root", pull, "TTree", "pull" );
	unsigned int freeNum = 28;//this->CountnDoF( finalSet ); // XXX: stupid
	double binz=1.;
	for( unsigned int i=0; i< theseDimensions.size(); ++i ) binz *= (double)theseDimensions[i].theseBins;
	double nDoF = binz - (double)freeNum;
	cout << "nDoF: " << nDoF << endl;
	cout << "chi2/nDoF: " << TotalChi2 / nDoF << endl;
	cout << "p-value: " << TMath::Prob( TotalChi2, nDoF ) << endl;
	cout << "The Total Chi2 = " << TotalChi2 << endl << endl;
	return;
}
double MultiDimChi2::CalcChi2( vector<double> expected_events, vector<double> observed_events, vector<double> errors )
{
	double Chi2Value=0.;
	double thisChi2=0.;
	for( unsigned int i=0; i< expected_events.size(); ++i )
	{
		cout << "O: " << observed_events[i] << "  E: " << expected_events[i] << "  Err: " << errors[i] << "  P: " << (observed_events[i]-expected_events[i])/errors[i] << endl;
		thisChi2 = expected_events[i] - observed_events[i] + observed_events[i] * log( observed_events[i] / expected_events[i] );
		if( !std::isnan(thisChi2) && !std::isinf(thisChi2) )
		{
			Chi2Value += thisChi2;
		}
	}
	return 2.*Chi2Value;
}
double MultiDimChi2::CalculateTotalExpected( vector<double> thisBinCenter )
{
	double result_AllPDFs = 0.;
	cout << "Calculating Expected Number of Events for this Bin" << endl;
	IPDF* thisPDF=NULL;
	IDataSet* thisDataSet=NULL;
	vector<ObservableRef> theseDim;
	vector<string> theseDimName;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		theseDim.push_back( ObservableRef( theseDimensions[i].ObservableName ) );
		theseDimName.push_back( theseDimensions[i].ObservableName );
	}
	for( unsigned int PDFNum=0; PDFNum< allPDFs.size(); ++PDFNum )
	{
		thisPDF = allPDFs[PDFNum];
		thisDataSet = allDataSets[PDFNum];
		RapidFitIntegrator* thisPDFIntegrator = thisPDF->GetPDFIntegrator();
		thisPDFIntegrator->SetPDF( thisPDF );
		thisPDFIntegrator->ForceTestStatus( false );
		vector<string> doNotIntegrate = thisPDF->GetDoNotIntegrateList();
		doNotIntegrate = StringProcessing::CombineUniques( doNotIntegrate, theseDimName );
		PhaseSpaceBoundary thisPhaseSpace( *allBoundaries[PDFNum] );
		PhaseSpaceBoundary thisPhaseSpace2( *allBoundaries[PDFNum] );
		vector<DataPoint> theseDataPoints;
		vector<DataPoint*> tempPoints = thisPhaseSpace.GetDiscreteCombinations();
		for( unsigned int i=0; i< tempPoints.size(); ++i ) /*if( i == 5 )*/ theseDataPoints.push_back( DataPoint( *tempPoints[i] ) );
		double thisResult=0.;
		cout << "\tI have " << theseDataPoints.size() << " integrals to perform for PhaseSpace: " << PDFNum+1 << endl;
		for( unsigned int combinationNum = 0; combinationNum< theseDataPoints.size(); ++combinationNum )
		{
			DataPoint thisDataPoint( theseDataPoints[combinationNum] );
			for( unsigned int i=0; i< theseDim.size(); ++i )
			{
				Observable newObs( string(theseDim[i]), thisBinCenter[i], "noUnits_Chi2" );
				thisDataPoint.SetObservable( theseDim[i], &newObs );
				double newMin = thisBinCenter[i]-0.5* theseDimensions[i].thisStepSize;
				double newMax = thisBinCenter[i]+0.5* theseDimensions[i].thisStepSize;
				thisPhaseSpace.SetConstraint( string(theseDim[i]), newMin, newMax, "noUnits_Chi2" );
				//cout << string(theseDim[i]) << "  " << newMin << "  " << newMax << endl;
			}
			vector<string> discNames = thisPhaseSpace.GetDiscreteNames();
			for( unsigned int i=0; i< discNames.size(); ++i )
			{
				double value = thisDataPoint.GetObservable( discNames[i] )->GetValue();
				thisPhaseSpace.SetConstraint( discNames[i], vector<double>(1,value), "noUnits_Chi2" );
				thisPhaseSpace2.SetConstraint( discNames[i], vector<double>(1,value), "noUnits_Chi2" );
			}
			cout << "\tCalculating Integral: " << combinationNum +1 << endl;
			thisPDFIntegrator->ForceTestStatus( true );
			thisDataPoint.SetPhaseSpaceBoundary( &thisPhaseSpace );
			double Total = thisPDF->Evaluate( &thisDataPoint );// thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace2, thisPDF->GetDoNotIntegrateList() );
			double Integral = 1.;//thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace, thisPDF->GetDoNotIntegrateList() );
			cout << "\tIntegral = " << Integral << endl;//" : " << Integral2 << " : " << Integral3 << endl;
			cout << "\tTotal = " << Total << endl;
			cout << "\tScaling Integral to DataSet." << endl;
			double PDF2DataNorm = this->CorrectYield( thisDataSet, &thisDataPoint );
			cout << "\tScale = " << PDF2DataNorm << endl;
			double thisYield = (Integral/Total) * PDF2DataNorm;
			cout << "\tExpected Yield for thisCombination = " << thisYield << endl;
			thisResult += thisYield;
		}
		result_AllPDFs += thisResult;
	}
	cout << "\tReturning Expected Number of Events for this Coordinate" << endl;
	return result_AllPDFs;
}
double MultiDimChi2::CorrectYield( IDataSet* thisSet, DataPoint* thisPoint )
{
	double total_yield = 0.;
	vector<DataPoint> thesePoints = thisSet->GetDiscreteSubSet( NULL );
	if( thisSet->GetWeightsWereUsed() )
	{
		ObservableRef thisRef( thisSet->GetWeightName() );
		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			bool isInRange = true;
			if( isInRange )
			{
				total_yield += thesePoints[i].GetObservable( thisRef )->GetValue();//GetEventWeight();//GetObservable( thisRef )->GetValue();
			}
		}
	}
	else
	{
		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			bool isInRange = true;
			for( unsigned int j=0; j< theseDimensions.size(); ++j )
			{
				double center = thisPoint->GetObservable( theseDimensions[j].ObservableName )->GetValue();
				double StepSize = theseDimensions[j].thisStepSize;
				double new_min = center-0.5*StepSize;
				double new_max = center+0.5*StepSize;
				double thisPointVal = thesePoints[i].GetObservable( theseDimensions[j].ObservableName )->GetValue();
				if( thisPointVal > new_max || thisPointVal < new_min )
				{
					isInRange = false;
					break;
				}
			}
			if( isInRange )
			{
				++total_yield;
			}
		}
	}
	return total_yield;
}
void MultiDimChi2::ConstructAllCoordinates()
{
	vector<double> thisBinCenter( theseDimensions.size(), 0. );
	unsigned int total_points = 1;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		total_points*=theseDimensions[i].theseBins;
	}
	allBinCenters.assign( total_points, thisBinCenter );
	for( unsigned int i=0; i< nDim; ++i )
	{
		this->AddCoordinates( i );
	}
}
void MultiDimChi2::AddCoordinates( unsigned int thisDim )
{
	cout << "Adding Coordinated for dimension: " << thisDim << endl;
	unsigned int number_of_set_repeats = 1;
	for( unsigned int i=0; i< thisDim; ++i )
	{
		number_of_set_repeats *= theseDimensions[i].theseBins;
	}
	cout << "There are: " << (thisDim-0) << " dimensions 'outside' this one." << endl;
	unsigned int number_of_individual_repeats = 1;
	for( unsigned int i=thisDim+1; i< nDim; ++i )
	{
		cout << i << "  " << theseDimensions[i].theseBins << endl;
		number_of_individual_repeats *= theseDimensions[i].theseBins;
	}
	cout << "There are: " << ((nDim-1)-thisDim) << " dimensions 'inside' this one." << endl;
	cout << number_of_set_repeats << "  x  " << number_of_individual_repeats << endl;
	unsigned int global_count=0;
	for( unsigned int i=0; i< number_of_set_repeats; ++i )
	{
		for( unsigned int j=0; j< theseDimensions[thisDim].theseBins; ++j )
		{
			for( unsigned int k=0; k< number_of_individual_repeats; ++k )
			{
				allBinCenters[global_count][thisDim] = theseDimensions[thisDim].binCenters[j];
				++global_count;
			}
		}
	}
}
void MultiDimChi2::ConstructBoundaries( PhaseSpaceBoundary* totalPhaseSpace, vector<string> wanted_observables )
{
	(void)wanted_observables;
	for( unsigned int i=0; i< allPDFs.size(); ++i )
	{
		vector<string> allDescribedObservables = allPDFs[i]->GetPrototypeDataPoint();
		allBoundaries.push_back( new PhaseSpaceBoundary( allDescribedObservables ) );
		PhaseSpaceBoundary* thisPhaseSpace = allBoundaries.back();
		for( unsigned int j=0; j< allDescribedObservables.size(); ++j )
		{
			IConstraint* thisConst = ClassLookUp::CopyConstraint( totalPhaseSpace->GetConstraint( allDescribedObservables[j] ) );
			thisPhaseSpace->SetConstraint( allDescribedObservables[j], thisConst );
			delete thisConst;
		}
		allPDFs[i]->ChangePhaseSpace( thisPhaseSpace );
		allDataSets[i]->SetBoundary( thisPhaseSpace );
	}
}
double MultiDimChi2::CalculateRange( PhaseSpaceBoundary* thisBound )
{
	(void)thisBound;
	double range = 1.;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		double this_min = x_min[i];//thisConst->GetMinimum();
		double this_max = x_max[i];//thisConst->GetMaximum();
		range *= fabs( this_max - this_min );
	}
	return range;
}
double MultiDimChi2::PDF2DataNormalisation( unsigned int PDFDataNum, const unsigned int combinationIndex, DataPoint* thisDataPoint )
{
	IPDF* thisPDF = allPDFs[ PDFDataNum ];
	IDataSet* thisDataSet = allDataSets[ PDFDataNum ];
	PhaseSpaceBoundary* thisBound = allBoundaries[ PDFDataNum ];
	double normalisation=1.;
	PhaseSpaceBoundary thisBound3( *thisBound );
	vector<string> wantedParams;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		wantedParams.push_back( theseDimensions[i].ObservableName );
		double value = thisDataPoint->GetObservable( wantedParams.back() )->GetValue();
		thisBound3.SetConstraint( theseDimensions[i].ObservableName, vector<double>(1, value), "unitless" );
	}
	vector<string> fixed_param = thisBound3.GetDiscreteNames();
	for( unsigned int i=0; i< fixed_param.size(); ++i )
	{
		double value = thisDataPoint->GetObservable( fixed_param[i] )->GetValue();
		thisBound3.SetConstraint( fixed_param[i], vector<double>(1, value), "unitless" );
	}
	vector<string> doNotList = thisPDF->GetDoNotIntegrateList();
	doNotList = StringProcessing::CombineUniques( doNotList, wantedParams );
	RapidFitIntegrator* pdfIntegrator = thisPDF->GetPDFIntegrator();
	pdfIntegrator->ForceTestStatus( false );
	double thisRatio =  pdfIntegrator->TestIntegral( thisDataPoint, &thisBound3, doNotList );
	cout << "\t\tScaling Based on Numerical/Analytical Ratio: " << 1./thisRatio << endl;
	normalisation /= thisRatio;                    //      Attempt to correct for Numerical != analytical due to any constant factor due to numerical inaccuracy
	//      (some constant close to 1. exactly 1. for numerical PDFs)
	vector<DataPoint*> allCombinations = thisBound->GetDiscreteCombinations();
	double total_yield=0.;
	if( thisDataSet->GetWeightsWereUsed() )
	{
		vector<DataPoint> thesePoints = thisDataSet->GetDiscreteSubSet( allCombinations[combinationIndex] );
		cout << thesePoints.size() << endl;
		ObservableRef thisWeight( thisDataSet->GetWeightName() );
		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			total_yield += thesePoints[i].GetObservable( thisWeight )->GetValue();
		}
	}
	else
	{
		total_yield = thisDataSet->GetDataNumber( allCombinations[combinationIndex] );
	}
	cout << total_yield << endl;
	normalisation *= total_yield;
	PhaseSpaceBoundary thisBound2( *thisBound );
	for( unsigned int AxisNum=0; AxisNum< theseDimensions.size(); ++AxisNum )
	{
		double ThisStep = theseDimensions[AxisNum].thisStepSize;
		double centralValue = thisDataPoint->GetObservable( theseDimensions[AxisNum].ObservableName )->GetValue();
		double newMax = centralValue+0.5*ThisStep;
		double newMin = centralValue-0.5*ThisStep;
		cout << "\t\tConstructing New Constraint for: " << theseDimensions[AxisNum].ObservableName << "\t" << newMin << " : " << newMax << endl;
		thisBound2.SetConstraint( theseDimensions[AxisNum].ObservableName, newMin, newMax, "noUnit_MultiDimChi2" );
	}
	cout << "\t\tCalculating Bin Integral" << endl;
	return normalisation;
}

