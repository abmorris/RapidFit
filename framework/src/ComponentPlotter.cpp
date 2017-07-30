/*!
 * @class ComponentPlotter
 *
 * A class for plotting PDF projections onto histograms
 *
 * @author Robert Currie rcurrie@cern.ch
 */
///	ROOT Headers
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TFolder.h"
#include "TTree.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TPad.h"
///	RapidFit Headers
#include "ComponentPlotter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"
#include "PhaseSpaceBoundary.h"
#include "ClassLookUp.h"
#include "ComponentRef.h"
#include "ObservableDiscreteConstraint.h"
#include "RapidFitRandom.h"
#include "MultiDimChi2.h"
#include "NormalisedSumPDF.h"
///	System Headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <cstdlib>
#include <algorithm>
#ifdef __USE_VALGRIND
#include <valgrind/callgrind.h>
#endif
#define DOUBLE_TOLERANCE 1E-6
using std::stringstream;
using std::setw;
using std::setprecision;
using std::flush;
using std::left;
//Constructor with correct arguments
ComponentPlotter::ComponentPlotter( IPDF * NewPDF, IDataSet * NewDataSet, TString PDFStr, TFile* filename, string ObservableName, CompPlotter_config* config, int PDF_Num ) :
	observableName( ObservableName ), weightName(), plotPDF( ClassLookUp::CopyPDF(NewPDF) ), plotData( NewDataSet ),
	pdfIntegrator( NULL ), weightsWereUsed(false), weight_norm(1.), discreteNames(), continuousNames(), full_boundary( NULL ), PlotFile( filename ),
	total_points( (config!=NULL)?config->PDF_points:128 ), data_binning( (config!=NULL)?config->data_bins:100 ), pdfStr( PDFStr ),
	logY( (config!=NULL)?config->logY:false ), logX( (config!=NULL)?config->logX:false ), this_config( config ), boundary_min( -99999 ), boundary_max( -99999 ), step_size( -99999 ),
	onlyZero( (config!=NULL)?config->OnlyZero:false ), combination_integral(vector<double>()), ratioOfIntegrals(1,1.), wanted_weights(), format(),
	data_subsets(), allCombinations(), combinationWeights(), combinationDescriptions(), observableValues(), binned_data(), total_components(),
	chi2(), N(), allPullData(), PDFNum(PDF_Num), initialBoundary(NULL)
{
	initialBoundary = new PhaseSpaceBoundary( *plotData->GetBoundary() );
	vector<string> allDescribedObservables = plotPDF->GetPrototypeDataPoint();
	full_boundary = new PhaseSpaceBoundary( allDescribedObservables );
	for( unsigned int i=0; i< allDescribedObservables.size(); ++i )
	{
		IConstraint* thisConst = initialBoundary->GetConstraint( allDescribedObservables[i] );
		full_boundary->SetConstraint( allDescribedObservables[i], thisConst );
	}
	plotData->SetBoundary( full_boundary );
	plotPDF->ChangePhaseSpace( full_boundary );
	TH1::SetDefaultSumw2(true);
	plotPDF->TurnCachingOff();
	plotPDF->SetComponentStatus( true );
	pdfIntegrator = new RapidFitIntegrator( plotPDF );
	pdfIntegrator->SetPDF( plotPDF );
	if( pdfIntegrator->GetUseGSLIntegrator() ) cout << "Using GSL for projections" << endl;
	RapidFitIntegratorConfig* projectionIntegratorConfig = NULL;
	if( config != NULL ) projectionIntegratorConfig = config->integratorConfig;
	if( projectionIntegratorConfig != NULL )
	{
		RapidFitIntegratorConfig* thisIntConfig = new RapidFitIntegratorConfig( *projectionIntegratorConfig );
		pdfIntegrator->SetUpIntegrator( thisIntConfig );
		plotPDF->SetUpIntegrator( thisIntConfig );
		delete thisIntConfig;
	}
	//Make the histogram of this observable
	format = new EdStyle();
	format->SetStyle();
	gStyle->SetPadLeftMargin( (Float_t)0.15 );
	gStyle->SetTitleOffset((Float_t)0.9,"Y");
	boundary_min = full_boundary->GetConstraint( observableName )->GetMinimum();
	boundary_max = full_boundary->GetConstraint( observableName )->GetMaximum();
	if(config->xmax > -99999) boundary_max = config->xmax;
	if(config->xmin > -99999) boundary_min = config->xmin;
	step_size = ( boundary_max - boundary_min ) / (double)( total_points - 1 );
	//Work out what to plot
	vector<DataPoint*> allCombinations_input;
	if( !config->plotAllCombinations )
	{
		for( const auto& DiscreteObs: full_boundary->GetDiscreteNames() )
		{
			if( full_boundary->GetConstraint( DiscreteObs )->GetValues().size() > 1 )
			{
				full_boundary->RemoveConstraint( DiscreteObs );
				ObservableDiscreteConstraint thisConstraint( DiscreteObs, {config->defaultCombinationValue}, " ", "" );
				full_boundary->AddConstraint( DiscreteObs, &thisConstraint );
			}
		}
		plotPDF->ChangePhaseSpace( full_boundary );
	}
	allCombinations_input = full_boundary->GetDiscreteCombinations();
	if( config->ForceCombinationNumber != -1 )
	{
		cout << "Requested ONLY to use Combination Number: " << config->ForceCombinationNumber << endl;
		if( config->ForceCombinationNumber > (int)allCombinations_input.size() )
		{
			cout << "CANNOT USE Combination Number: " << config->ForceCombinationNumber << " Ignoring!" << endl;
		}
		allCombinations.push_back(*allCombinations_input[config->ForceCombinationNumber]);
	}
	else
	{
		for( auto thisCombination : allCombinations_input)
		{
			double thisNum = plotData->GetDataNumber( thisCombination );
			if( fabs(thisNum) > 1E-5 )
				allCombinations.push_back( *thisCombination );
		}
	}
	cout << endl << "All Combinations with Data:" << endl;
	cout << allCombinations.size() << endl << endl;
	SetupCombinationDescriptions();
	if( config != NULL )
	{
		if( !config->combination_names.empty() ) combinationDescriptions = config->combination_names;
	}
	// Hack to get the internal cached values for NormalisedSumPDF to work
	if(plotPDF->GetName()=="NormalisedSumPDF")
	{
		int combination_counter = 0;
		for(auto& combination: allCombinations)
		{
			double integral = plotPDF->GetPDFIntegrator()->Integral( &combination, full_boundary );
			std::array<double,2> CachedIntegrals = ((NormalisedSumPDF*)plotPDF)->GetCachedIntegrals( &combination, full_boundary );
			std::cout << "Combination " << combination_counter << "\n\ttotal: " << integral << ", PDF 1: " << CachedIntegrals[0] << ", PDF 2: " << CachedIntegrals[1] << endl;
			combination_counter++;
		}
		
	}
	for( unsigned int i=0; i< allCombinations.size(); ++i )
	{
		pdfIntegrator->ForceTestStatus( false );
		allCombinations[i].SetPhaseSpaceBoundary( full_boundary );
		double thisIntegral = 0.;
		try
		{
			allCombinations[i].SetPhaseSpaceBoundary( full_boundary );
			thisIntegral = pdfIntegrator->NumericallyIntegrateDataPoint( &allCombinations[i], full_boundary, plotPDF->GetDoNotIntegrateList() );
		}
		catch(...)
		{
			thisIntegral = 1.;
			cout << endl << "CANNOT PROPERLY NORMALISE WHOLE PDF, THIS WILL LEAD TO NORMALISATION ISSUES OVER THE WHOLE PDF" << endl << endl;
		}
		combination_integral.push_back( thisIntegral );
		if( plotPDF->GetNumericalNormalisation() == true )
			ratioOfIntegrals.push_back( 1. );
		else
		{
			if( config->ScaleNumerical )
				ratioOfIntegrals.push_back( pdfIntegrator->GetRatioOfIntegrals() );
			else
				ratioOfIntegrals.push_back( 1./pdfIntegrator->GetRatioOfIntegrals() );
		}
	}
	if( ratioOfIntegrals.size() == 0 || ratioOfIntegrals.empty() )
		ratioOfIntegrals =vector<double> ( 1, 1. );
	vector<double> minimum, maximum;
	vector<int> binNumber;
	(void) minimum; (void) maximum; (void) binNumber;
	if( config != NULL )
	{
		if( config->component_names.empty() )
		{
			if( plotPDF->GetName() == plotPDF->GetLabel() )
				config->component_names.push_back( plotPDF->GetName() );
			else
				config->component_names.push_back( "Total" );
			vector<string> pdfComponents = plotPDF->PDFComponents();
			pdfComponents = StringProcessing::MoveElementToStart( pdfComponents, "0" );
			for( unsigned int i=0; i< pdfComponents.size(); ++i )
			{
				ComponentRef thisRef( pdfComponents[i], observableName );
				if( pdfComponents[i] != "0" )
					config->component_names.push_back( plotPDF->GetComponentName( &thisRef ) );
			}
		}
	}
	TH1::SetDefaultSumw2( true );
	pdfIntegrator->ForceTestStatus( true );
	plotData->SetBoundary( initialBoundary );
}
//Destructor
ComponentPlotter::~ComponentPlotter()
{
	// XXX many/all of these objects are copyable. what's with all the pointer abuse?
	if( plotPDF != NULL )	delete plotPDF;
	if( pdfIntegrator != NULL ) delete pdfIntegrator;
	delete format;
	if( initialBoundary != NULL ) delete initialBoundary;
	if( full_boundary != NULL ) delete full_boundary;
}
//Create a root file containing a projection plot over one observable
void ComponentPlotter::ProjectObservable()
{
#ifdef __USE_VALGRIND
	CALLGRIND_START_INSTRUMENTATION;
#endif
	plotData->SetBoundary( full_boundary );
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();
	plotPDF->UnsetCache();
	//Check the observable can be plotted
	bool continuous = !( plotData->GetBoundary()->GetConstraint( observableName )->IsDiscrete() );
	bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &(observableName) ) == -1 );
	if( continuous && doIntegrate )
	{
		//Do the projecting of the pdfs
		cout << "Projecting " << observableName << endl;
		GenerateProjectionData();
	}
	else
	{
		cerr << "CANNOT PROJECT: " << observableName << endl;
		return;
	}
	plotData->SetBoundary( initialBoundary );
#ifdef __USE_VALGRIND
	CALLGRIND_STOP_INSTRUMENTATION;
	CALLGRIND_DUMP_STATS;
#endif
}
//	Returns a vector of arrays of dimention		AllCombinations (+ 0 combination) * total_points
//	This data corresponds to the X axis of the projected functions from the PDF
//	This is a simple loop of adding numbers onto a running total and recording the running total between minima and maxima in a vector
vector<vector<double>> ComponentPlotter::MakeXProjectionData( unsigned int num_combinations )
{
	vector<vector<double>> new_dataarray;
	unsigned int total = num_combinations;
	if( total > 1 ) ++total;
	for (unsigned int combinationIndex = 0; combinationIndex < total; ++combinationIndex )
	{
		vector<double> observableValueArray;
		for(int pointIndex = 0; pointIndex < total_points; ++pointIndex )
			//	Start at minimum of observable in phase-space
			//	move a step equal to the size of the interval that you should take for this subset of data
			//	due to the fact that the stepsize corresponding to the tag=-1 (combinationIndex=0)
			//	may be different ot the stepsize corresponding to the tag=1 (combinationIndex=1)
			observableValueArray.push_back(boundary_min + ( step_size * (pointIndex) ));
		new_dataarray.push_back( observableValueArray );
	}
	return new_dataarray;
}
//      This returns the Y component of the projections for all combinations and for all components for each projection
//      This calls ComponentPlotter::ProjectObservableComponent for each combination
//      Combination 0 is constricted as the total of the individual components once all components for all sub component have been calculated
vector<vector<double>> ComponentPlotter::MakeYProjectionData( string component_str )
{
	vector<vector<double>> new_dataarray;
	//	Loop over all discrete combinations for this PDF configuration
	for( unsigned int combinationIndex = 0; combinationIndex < allCombinations.size(); ++combinationIndex )
	{
		cout << "Constructing PDF Integral of: " << observableName << " Combination: " << combinationIndex+1 << " of " << allCombinations.size() << ". For component: " << component_str <<" .\t";
		allCombinations[combinationIndex].SetPhaseSpaceBoundary( full_boundary );
		cout << observableName << ": " << boundary_min << " <-> " << boundary_min + ( step_size * (total_points-1) ) << endl;
		cout << "Starting Projection: " << combinationIndex+1 << " of " << allCombinations.size() <<"."<< endl;
		//Get the raw probability distribution for this component from the PDF
		//Calculate the projection for this combination
		vector<double> projectionValueArray =
			//	DataPoint with all information on this configuration,
			ProjectObservableComponent( &allCombinations[combinationIndex],
					//	minimum of range in Observable,	num of points in range,	plot_interval,	component of interest
					observableName, boundary_min, total_points, step_size, component_str );
		//Normalise the PDF such that the final decision of bin number is the only variable missing.
		//Update the data average values, and make the projection graph arrays
		double PDFNormalisation = PDF2DataNormalisation( combinationIndex );
		//	Perform Normalisation
		std::transform(projectionValueArray.begin(),projectionValueArray.end(),projectionValueArray.begin(),std::bind1st(std::multiplies<double>(), PDFNormalisation));
		cout << "Finished Combination: " << combinationIndex+1 << " of " << allCombinations.size() <<"."<< endl;
		new_dataarray.push_back( projectionValueArray );
	}
//	When we have more than 1 discrete component we need to create component 0 which contains the total PDF result at this coordinate
	if( new_dataarray.size() > 1 )
	{
		//      Construct 0th combination and fill it as the sum of all sub-combinations of the dataset
		vector<double> zero_component;
		for( unsigned int i=0; i< (unsigned)total_points; ++i )
		{
			//	Zeroth component is defined as the total of each of the sub components at this point
			zero_component.push_back(0);
			for( unsigned int j=0; j< new_dataarray.size(); ++j )
				zero_component[i]+= new_dataarray[j][i];
		}
		new_dataarray.insert(new_dataarray.begin(), zero_component );
	}
	return new_dataarray;
}
//	Generate Projection Data for a given Observable
//	When this function starts it collects statistical information and infromation on the boundary of the observable
//	From here it loops over all components and combinations for this PDFWithData object
//	For configurations with more than 1 combination I create combination 0 at the end to contain the total sum of all of the unique combination datasets and projections
void ComponentPlotter::GenerateProjectionData()
{
	//	Now we have:
	//			plot phase-space
	//	Lets plot the projection
	//	Check if the 0'th component has been defined and add it if it's missing.
	vector<string> PDF_Components;
	if( onlyZero == true )
		PDF_Components = vector<string>(1,"0");
	else
	{
		PDF_Components = plotPDF->PDFComponents();
		if( StringProcessing::VectorContains( PDF_Components, string("0") ) == -1 )
			PDF_Components.insert(PDF_Components.begin(),"0");
	}
	cout << endl << "Components: " << PDF_Components.size() << endl;
	// This is still quite bad, but you should see the old version
	vector<vector<vector<double>>> X_values;
	vector<vector<vector<double>>> Y_values;
	//	Loop Over ALL components that are provided by this PDF
	for( unsigned int i=0; i< PDF_Components.size(); ++i )
	{
		cout << "\n\t\tCOMPONENT: " << i+1 << " of: " << PDF_Components.size() << "\t\tPDF: "<< plotPDF->GetLabel() << endl <<endl;
		//	Generate and store the X and Y values of the projection plot
		X_values.push_back( MakeXProjectionData( (unsigned)allCombinations.size() ) );
		Y_values.push_back( MakeYProjectionData( PDF_Components[i] ) );
	}
	//	Write the output to the output file
	TDirectory* here = gDirectory;
	TString PDFDesc; PDFDesc+=PDFNum;
	vector<string> combDescs=combinationDescriptions;
	for( unsigned int i=0; i< combDescs.size(); ++i )
	{
		combDescs[i].append("_PDF_");
		combDescs[i].append(PDFDesc.Data());
	}
	cout << endl << "Writing Combinations:" << endl;
	for( unsigned int i=0; i< combinationDescriptions.size(); ++i ) cout << combinationDescriptions[i] << endl;
	cout << endl;
	WriteOutput( X_values, Y_values, combinationDescriptions );
	here->cd();
	return;
}
//	This routine bins the data for a requested combinationNumber into a TH1, the number of bins is determined through the data_binning int
//
TH1* ComponentPlotter::FormatData( unsigned int combinationNumber )
{
	vector<DataPoint> wanted_points = data_subsets[combinationNumber];
	TString component_num_str;component_num_str+=combinationNumber;
	TH1* returnable = new TH1D( "Data_For_Component_"+component_num_str, "Data_For_Component_"+component_num_str, data_binning, boundary_min, boundary_max );
	vector<double> wanted_data, this_wanted_weights;
	ObservableRef new_ref( observableName );
	for( auto& point_i : wanted_points )
	{
		wanted_data.push_back( point_i.GetObservable( new_ref )->GetValue() );
	}
	if( wanted_weights.empty() && weightsWereUsed ) wanted_weights = vector<double>( wanted_data.size(), 1. );	//	'Should' never occur... but better than a segfault
	this_wanted_weights = wanted_weights;
	//	Fill a histogram with the weighted data, with weights normalised to 1
	//	Think of this as multiplying by a very clever version of the number 1 ;)
	//	This gives the correct normalisation
	vector<double>::iterator point_i = wanted_data.begin();
	vector<double>::iterator weight_i = this_wanted_weights.begin();
	if( weightsWereUsed )
		for( ; point_i != wanted_data.end(); ++point_i, ++weight_i )
			returnable->Fill( *point_i, *weight_i );
	else
		for( ; point_i != wanted_data.end(); ++point_i )
			returnable->Fill( *point_i, 1. );
	return returnable;
}
//	This routine plots all of the data for all combinations and stores the data int TTree objects in a .root file
void ComponentPlotter::WriteOutput( vector<vector<vector<double>>>& X_values, vector<vector<vector<double>>>& Y_values, vector<string> CombinationDescriptions )
{
	//	Make sure we start in the correct place
	PlotFile->cd();
	//	Seperate directories per PDF in the XML
	if( gDirectory->GetDirectory( pdfStr ) == 0 )	gDirectory->mkdir( pdfStr );
	gDirectory->cd( pdfStr );
	//	Seperate directories per observable in the XML
	if( gDirectory->GetDirectory( observableName.c_str() ) == 0 )	gDirectory->mkdir( observableName.c_str() );
	gDirectory->cd( observableName.c_str() );
	//	Pointer to here
	TDirectory* PlotDirectory = gDirectory;
	//	Loop over all components
	for( unsigned int componentIndex=0; componentIndex < X_values.size(); ++componentIndex )
	{
		TString componentName("Component_");componentName+=componentIndex;
		if( gDirectory->GetDirectory( componentName ) == 0 )	gDirectory->mkdir( componentName );
		gDirectory->cd( componentName );
		TDirectory* componentDir = gDirectory;
		//	Loop over all combinations for this component
		for( unsigned int combinationIndex=0; combinationIndex < X_values[componentIndex].size(); ++combinationIndex )
		{
			TString combinationName("Combination_");combinationName+=combinationIndex;
			if( gDirectory->GetDirectory( combinationName ) == 0 )	gDirectory->mkdir( combinationName );
			gDirectory->cd( combinationName );
			//	Save all of the X and Y data in the TTree
			TString TTree_name( "TTree_" ); TTree_name+=componentIndex; TTree_name.Append("_"); TTree_name+=combinationIndex;
			TString TTree_title( TTree_name );
			TTree_name.Append("_");TTree_name += RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
			string TTreeName( TTree_name.Data() );
			replace( TTreeName.begin(), TTreeName.end(), '.', '_' );
			TTree* this_data = new TTree( TTreeName.c_str(), TTree_title );
			vector<double> unnormalised;
			vector<double> normalisation;
			double Normalisation = 0.;
			if( X_values[componentIndex].size() > 1 )
			{
				if( combinationIndex == 0 )
					for( unsigned int i=0; i< X_values[componentIndex].size()-1; ++i )
						Normalisation += PDF2DataNormalisation( i );
				else
					Normalisation = PDF2DataNormalisation( combinationIndex-1 );
			}
			else
				Normalisation = PDF2DataNormalisation( 0 );
			for( const auto& num_i : Y_values[componentIndex][combinationIndex])
			{
				unnormalised.push_back( num_i / Normalisation );
				normalisation.push_back( Normalisation );
			}
			WriteBranch( this_data, "X_data", X_values[componentIndex][combinationIndex] );
			WriteBranch( this_data, "Y_data", Y_values[componentIndex][combinationIndex] );
			WriteBranch( this_data, "Y_data_unNormalised", unnormalised );
			WriteBranch( this_data, "Y_data_Normalisation", normalisation );
			this_data->Write( "", TObject::kOverwrite );
			TString Data_Name("Raw_Data_"); Data_Name+=componentIndex; Data_Name.Append("_");Data_Name+=combinationIndex;
			TString Data_Title( Data_Name );
			Data_Name.Append("_"); Data_Name += RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
			TTree* raw_data = new TTree( Data_Name, Data_Title );
			vector<double> real_raw_data;
			ObservableRef tempref( observableName );
			for( auto& point_i : data_subsets[combinationIndex] )
				real_raw_data.push_back( point_i.GetObservable( tempref )->GetValue() );
			WriteBranch( raw_data, "Value", real_raw_data );
			//	Bin the data for this combination
			TH1* data_plot = FormatData( combinationIndex );
			TString ext("_"); ext+=componentIndex; ext.Append("_"); ext+=combinationIndex;
			ext.Append("_"); ext+=RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
			string cleanExt( ext.Data() );
			replace( cleanExt.begin(), cleanExt.end(), '.', '_' );
			data_plot->SetName( data_plot->GetName() + TString(cleanExt.c_str()) );
			data_plot->SetTitle( "" );
			//	Plot the component on this combination and save the TCanvas
			TString Graph_Name("TGraph_");Graph_Name+=componentIndex;
			Graph_Name.Append("_");Graph_Name+=combinationIndex;
			Graph_Name.Append("_");Graph_Name+=RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
			string graphCleanName( Graph_Name.Data() );
			replace( graphCleanName.begin(), graphCleanName.end(), '.', '_' );
			TGraph* data_graph = new TGraph( total_points,  X_values[componentIndex][combinationIndex].data(),  Y_values[componentIndex][combinationIndex].data() );
			data_graph->SetTitle( "" ); data_graph->SetName( graphCleanName.c_str() );
			TString Canvas_Name("TCanvas_");Canvas_Name+=componentIndex;
			Canvas_Name.Append("_");Canvas_Name+=combinationIndex;
			Canvas_Name.Append("_");Canvas_Name+=RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
			string CanvasCleanName( Canvas_Name.Data() );
			replace( CanvasCleanName.begin(), CanvasCleanName.end(), '.', '_' );
			TCanvas* c1 = EdStyle::RapidFitCanvas( CanvasCleanName.c_str(), "" );
			data_plot->Draw();
			data_graph->Draw("PC SAME");
			c1->Update();
			double Y_min = -99999.;
			double Y_max = -99999.;
			if( this_config != NULL )
			{
				Y_min = this_config->ymin;
				Y_max = this_config->ymax;
			}
			c1->Update();
			if( Y_min <= -99999. ) Y_min = data_graph->GetYaxis()->GetXmin();
			if( Y_max <= -99999. ) Y_max = data_graph->GetYaxis()->GetXmax();
			data_graph->GetYaxis()->SetRangeUser( Y_min, Y_max );
			c1->Update();
			c1->Write();
			componentDir->cd();
		}
		PlotDirectory->cd();
	}
	cout << "Data Written to File, Making Graphs" << endl;
	//	Starting at the top of the file again
	PlotFile->cd();
	gDirectory->cd( pdfStr );
	if( gDirectory->GetDirectory( "overlay_graphs" ) == 0 )	gDirectory->mkdir( "overlay_graphs" );
	gDirectory->cd( "overlay_graphs" );
	//	Overlay all components for each combination
	//
	//	This is painful as your looping over the upper 2D of a 3D vector inside out...
	//	Sorry but I'm not rewriting the whole class or inverting the vector of vectors as we will just get more mistakes
	//
	//	This definitely works and isn't quite so bad as I feares but knowing which component we are using is a pain
	//
	for( unsigned int combinationIndex=0; combinationIndex < X_values[0].size(); ++combinationIndex )
	{
		vector<TGraph*> these_components;
		//	Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
		for( unsigned int componentIndex=0; componentIndex < X_values.size(); ++componentIndex )
		{
			TGraph* data_graph = new TGraph( total_points, X_values[componentIndex][combinationIndex].data(),  Y_values[componentIndex][combinationIndex].data() );
			string data_name = "data_" + std::to_string(combinationIndex) + "_" + std::to_string(componentIndex) + "_" + std::to_string(RapidFitRandom::GetFrameworkRandomFunction()->Rndm());
			replace( data_name.begin(), data_name.end(), '.', '_' );
			data_graph->SetName( data_name.c_str() ); data_graph->SetTitle( data_name.c_str() );
			data_graph->SetLineColor( (Color_t)(componentIndex+1) );
			data_graph->SetLineWidth( (Width_t)3 );
			data_graph->SetMarkerColor( (Color_t)(componentIndex+1) );
			these_components.push_back( data_graph );
		}
		total_components.push_back( these_components );
		TH1* data_plot = FormatData( combinationIndex );
		TString ext("_"); ext+=combinationIndex; ext.Append("_"); ext+=RapidFitRandom::GetFrameworkRandomFunction()->Rndm();
		string cleanExt( ext.Data() );
		replace( cleanExt.begin(), cleanExt.end(), '.', '_' );
		data_plot->SetName( data_plot->GetName() + TString(cleanExt.c_str()) );
		data_plot->SetTitle( "" );
		//	Object to hold the binned data for the life of this plot
		binned_data.push_back( new TGraphErrors( data_plot ) );
		binned_data.back()->Write("",TObject::kOverwrite);
		string desc;
		if( combinationIndex > 0 ) desc = CombinationDescriptions[combinationIndex-1];
		else if( combinationIndex == 0 ) desc = "All_Data";
		TString PDFDesc;PDFDesc+=PDFNum;
		desc.append("_PDF_"); desc.append(PDFDesc.Data());
		//	For the moment haven't decided if I should pass the global config to ALL sub plots, I shall get user input on this
		CompPlotter_config temp( *this_config );
		temp.logY = logY;
		temp.logX = logX;
		//	Static function so has to be told everything about what you want to plot!
		OutputPlot( binned_data.back(), these_components, observableName, desc, plotData->GetBoundary(), RapidFitRandom::GetFrameworkRandomFunction(), &temp );
		std::vector<DataPoint*> combs;
		if( this_config != NULL && combinationIndex == 0 )
			combs = full_boundary->GetDiscreteCombinations();
		else
			combs.push_back(&allCombinations[combinationIndex-1]);
		if( this_config->CalcChi2 == true && combinationIndex == 0 )
		{
			cout << endl;
			cout << "Calculating Chi^2:" << endl;
			std::vector<double> expected, observed;
			for( unsigned int i=0; i< (unsigned) binned_data[combinationIndex]->GetN(); ++i )
			{
				double bin_center = binned_data[combinationIndex]->GetX()[i];
				double bin_err = binned_data[combinationIndex]->GetEX()[i]; // This is half the bin width. I know this is stupid, but the data is in a TGraphErrors not a TH1 for some reason.
				expected.push_back(MultiDimChi2::CalculateExpected(*plotPDF, *full_boundary, combs, *plotData, {observableName}, {{bin_center-bin_err,bin_center+bin_err}}));
				observed.push_back(binned_data[combinationIndex]->GetY()[i]);
			}
			chi2 = MultiDimChi2::CalcChi2(expected, observed, false);
			//	dof = num Populated Bins - ndof in PDF
			N =  binned_data.back()->GetXaxis()->GetNbins();
			for( unsigned int i=0; i< (unsigned) binned_data[0]->GetN(); ++i ) if( fabs(binned_data[combinationIndex]->GetY()[i]) <= 0 ) --N;
			double n = (double) plotPDF->GetPhysicsParameters()->GetAllFloatNames().size();
			double denominator = (double)(N - n - 1. );
			cout << endl << "Chi^2/ndof :\t" << setprecision(10) << chi2 / denominator << endl;
		}
		if( this_config->DrawPull == true )
		{
			cout << endl;
			cout << "Calculating Pull Values" << endl;
			std::vector<double> localPullData;
			for( unsigned int i=0; i< (unsigned)binned_data[combinationIndex]->GetN(); ++i )
			{
				double bin_center = binned_data[combinationIndex]->GetX()[i];
				double bin_err = binned_data[combinationIndex]->GetEX()[i];
				double this_bin = MultiDimChi2::CalculateExpected(*plotPDF, *full_boundary, combs, *plotData, {observableName}, {{bin_center-bin_err,bin_center+bin_err}});// (*this)( bin_center );
				localPullData.push_back(this_bin);
			}
			TString desc_pull( desc );
			desc_pull.Append("_pull");
			cout << endl << "Making Plots" << endl;
			OutputPlot( binned_data.back(), these_components, observableName, string(desc_pull.Data()), plotData->GetBoundary(), RapidFitRandom::GetFrameworkRandomFunction(), &temp, localPullData );
			if( this_config != NULL && combinationIndex == 0 ) allPullData = localPullData;
		}
	}
	//	If there is more than 1 combination it's useful to plot the total's on the same graph with total component 0
	if( GetComponents().size() > 1 )
	{
		vector<TGraph*> allZerothComponents;
		vector<vector<TGraph*> > allGraphs = GetComponents();
		for( unsigned int i=0; i< allGraphs.size(); ++i )
		{
			allZerothComponents.push_back( allGraphs[i][0] );
		}
		string desc("_All_Combinations");
		TString PDFDesc;PDFDesc+=PDFNum;
		desc.append("_PDF_"); desc.append(PDFDesc.Data());
		vector<string> combDescs;
		if( this_config->combination_names.empty() )
		{
			combDescs.push_back( "Total");
			for( unsigned int i=0; i< combinationDescriptions.size(); ++i )
				combDescs.push_back( combinationDescriptions[i] );
		}
		CompPlotter_config temp( *this_config );
		vector<int> colorVec;
		for( unsigned int i=0; i< combinationDescriptions.size(); ++i )
			colorVec.push_back( (i+2) );
		if( temp.color_key.empty() )
			temp.color_key = colorVec;
		if( temp.combination_names.empty() )
			temp.component_names = combDescs;
		else
			temp.component_names = temp.combination_names;
		ComponentPlotter::OutputPlot( binned_data[0], allZerothComponents, observableName, desc, plotData->GetBoundary(), RapidFitRandom::GetFrameworkRandomFunction(), &temp );
		if( !allPullData.empty() )
		{
			desc.append("_wPulls");
			ComponentPlotter::OutputPlot( binned_data[0], allZerothComponents, observableName, desc, plotData->GetBoundary(), RapidFitRandom::GetFrameworkRandomFunction(), &temp, allPullData );
		}
	}
}
//	Again I won't be bothered coding this to check the input, this is the USERS problem!
vector<TGraph*> ComponentPlotter::MergeComponents( vector<vector<TGraph*> > input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( input.size() == 1 ) return input[0];
	vector<TGraph*> output_components;
	for( unsigned int component_i=0; component_i< input[0].size(); ++component_i )
	{
		//	Collect all component 0 into a vector or component 1 into a seperate vector etc...
		vector<TGraph*> this_component;
		for( unsigned int combination_i = 0; combination_i < input.size(); ++combination_i )
			this_component.push_back( input[combination_i][component_i] );
		//	Copy all of the data into a vector of vector<double>
		vector<vector<double> > X_val, Y_val;
		for( auto graph_i : this_component )
		{
			int data_num = graph_i->GetN();
			double* x_pointer = graph_i->GetX();
			double* y_pointer = graph_i->GetY();
			vector<double> this_x, this_y;
			for( int i=0; i < data_num; ++i )
			{
				this_x.push_back( x_pointer[i] );
				this_y.push_back( y_pointer[i] );
			}
			X_val.push_back( this_x );
			Y_val.push_back( this_y );
		}
		vector<double> final_X_val(X_val[0].size()), final_Y_val(Y_val[0].size());
		//      Sum the X and Y values
		for( unsigned int i=0; i< X_val[0].size(); ++i )
		{
			final_X_val[i] = X_val[0][i];
			final_Y_val[i] = 0.;
			for( unsigned int j=0; j< X_val.size(); ++j )
				final_Y_val[i] += Y_val[j][i];
		}
		//      Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
		TString TGraphName("TGraph_");TGraphName+=rand->Rndm();
		string TGRaphCleanName( TGraphName.Data() );
		replace( TGRaphCleanName.begin(), TGRaphCleanName.end(), '.', '_' );
		TGraph* output_graph = new TGraph( X_val[0].size(), final_X_val.data(), final_Y_val.data() );
		output_graph->SetName( TGRaphCleanName.c_str() );
		output_graph->SetTitle("");
		output_graph->SetLineColor( this_component[0]->GetLineColor() );
		output_graph->SetLineWidth( this_component[0]->GetLineWidth() );
		output_graph->SetLineStyle( this_component[0]->GetLineStyle() );
		output_graph->SetMarkerColor( this_component[0]->GetMarkerColor() );
		output_graph->SetMarkerStyle( this_component[0]->GetMarkerStyle() );
		output_components.push_back( output_graph );
	}
	return output_components;
}
//	I won't be bothered coding this in order to check that I have been provided with compatible input. The USER should check this!
TGraphErrors* ComponentPlotter::MergeBinnedData( vector<TGraphErrors*> input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;
	if( input.size() == 1 ) return input[0];
	vector<vector<double> > X_val, X_err, Y_val, Y_err;
	for( auto graph_i : input )
	{
		int data_num = graph_i->GetN();
		double* x_pointer = graph_i->GetX();
		double* y_pointer = graph_i->GetY();
		double* x_err_pointer = graph_i->GetEX();
		double* y_err_pointer = graph_i->GetEY();
		vector<double> this_x, this_y, this_ex, this_ey;
		for( int i=0; i < data_num; ++i )
		{
			this_x.push_back( x_pointer[i] );
			this_y.push_back( y_pointer[i] );
			this_ex.push_back( x_err_pointer[i] );
			this_ey.push_back( y_err_pointer[i] );
		}
		X_val.push_back( this_x );
		Y_val.push_back( this_y );
		X_err.push_back( this_ex );
		Y_err.push_back( this_ey );
	}
	vector<double> final_X_val(X_val[0].size()), final_X_err(X_val[0].size()), final_Y_val(X_val[0].size()), final_Y_err(X_val[0].size());
	//	Sum the X and Y values and add the errors in quadrature
	for( unsigned int i=0; i< X_err[0].size(); ++i )
	{
		double y_err_sq = 0.;
		final_X_val[i] = X_val[0][i];
		final_X_err[i] = X_err[0][i];
		final_Y_val[i] = 0.;
		for( unsigned int j=0; j< Y_err.size(); ++j )
		{
			final_Y_val[i] += Y_val[j][i];
			y_err_sq += Y_err[j][i]*Y_err[j][i];
		}
		final_Y_err[i] = sqrt( y_err_sq );
	}
	//      Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
	TString TGraphErrorsName("TGraphErrors_");TGraphErrorsName+=rand->Rndm();
	string TGraphErrorsCleanName( TGraphErrorsName.Data() );
	replace( TGraphErrorsCleanName.begin(), TGraphErrorsCleanName.end(), '.', '_' );
	TGraphErrors* output_graph = new TGraphErrors( X_val[0].size(), final_X_val.data(), final_Y_val.data(), final_X_err.data(), final_Y_err.data() );
	output_graph->SetName( TGraphErrorsCleanName.c_str() );
	output_graph->SetTitle("");
	output_graph->SetLineColor( input[0]->GetLineColor() );
	output_graph->SetLineWidth( input[0]->GetLineWidth() );
	output_graph->SetLineStyle( input[0]->GetLineStyle() );
	output_graph->SetMarkerColor( input[0]->GetMarkerColor() );
	output_graph->SetMarkerStyle( input[0]->GetMarkerStyle() );
	return output_graph;
}
//	Plot all components on this combinations and print and save the canvas
void ComponentPlotter::OutputPlot( TGraphErrors* input_data, vector<TGraph*> input_components, string observableName,
		string CombinationDescription, PhaseSpaceBoundary* total_boundary, TRandom* rand, CompPlotter_config* conf, vector<double> input_bin_theory_data )
{
	if( rand == NULL ) rand = gRandom;
	TString TCanvas_Name("Overlay_"+observableName+"_"+CombinationDescription+"_");TCanvas_Name+=rand->Rndm();
	string TCanvasCleanName( TCanvas_Name.Data() );
	replace( TCanvasCleanName.begin(), TCanvasCleanName.end(), '.', '_' );
	bool logy=false;
	bool logx=false;
	TCanvas* c1 = EdStyle::RapidFitCanvas( TCanvasCleanName.c_str(), "" );
	if( conf != NULL )
	{
		if( conf->logY )
		{
			logy=true;
			c1->SetLogy( true );
		}
		if( conf->logX )
		{
			logx=true;
			c1->SetLogx( true );
		}
	}
	TPad* pad1=NULL;	//		Holds Projection Plot
	TPad* pad2=NULL;	//		Holds Pull Plot
	if( !input_bin_theory_data.empty() )
	{
		pad1 = new TPad("pad1","pad1", 0., 0.3, 1., 1.);
		pad2 = new TPad("pad2","pad2", 0., 0.0, 1., 0.3);
		pad1->Draw();
		pad2->Draw();
		pad1->cd();
	}
	TString plotTitle;
	vector<int> Style_Key, Color_Key, Width_Key;
	vector<string> component_names;
	double X_min=-99999, X_max=-99999;
	double Y_min=-99999, Y_max=-99999;
	TString X_Title, Y_Title;
	double final_chi2=-99999;
	double legend_size=0.1;
	bool addLHCb=false;
	bool addRightLHCb=false;
	bool limitPulls=false;
	bool drawSpline=true;
	double XaxisTitleScale=1.;
	double XaxisLabelScale=1.;
	double YaxisTitleScale=1.;
	double YaxisLabelScale=1.;
	if( conf != NULL )
	{
		if( conf->logY )
		{
			logy=true;
			if( !input_bin_theory_data.empty() )
			{
				pad1->SetLogy( true );
			}
			else
			{
				c1->SetLogy( true );
			}
			//c1->Update();
		}
		if( conf->logX )
		{
			logx=true;
			if( !input_bin_theory_data.empty() )
			{
				pad1->SetLogx( true );
				pad2->SetLogx( true );
			}
			else
			{
				c1->SetLogx( true );
			}
		}
		plotTitle = conf->PlotTitle;
		Style_Key = conf->style_key;
		Color_Key = conf->color_key;
		Width_Key = conf->width_key;
		component_names = conf->component_names;
		legend_size = conf->LegendTextSize;
		addLHCb = conf->addLHCb;
		addRightLHCb = conf->addRightLHCb;
		X_min = conf->xmin;
		X_max = conf->xmax;
		Y_min = conf->ymin;
		Y_max = conf->ymax;
		X_Title = conf->xtitle;
		Y_Title = conf->ytitle;
		final_chi2 = conf->Chi2Value;
		limitPulls = conf->LimitPulls;
		drawSpline = conf->useSpline;
		XaxisTitleScale = conf->XaxisTitleScale;
		XaxisLabelScale = conf->XaxisLabelScale;
		YaxisTitleScale = conf->YaxisTitleScale;
		YaxisLabelScale = conf->YaxisLabelScale;
	}
	input_data->SetTitle( plotTitle );
	input_data->Draw("AP");
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}
	if( !input_bin_theory_data.empty() )
	{
		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}
	c1->Update();
	if( X_min <= -99999 ) X_min = total_boundary->GetConstraint( observableName )->GetMinimum();
	if( X_max <= -99999 ) X_max = total_boundary->GetConstraint( observableName )->GetMaximum();
	if( Y_min <= -99999 ) Y_min = input_data->GetYaxis()->GetXmin();//logy==true?0.5:0.;
	if( Y_max <= -99999 ) Y_max = input_data->GetYaxis()->GetXmax();
	TString Y_ext(" / ( ");
	stringstream thisStream;
	if( conf != NULL )
	{
		thisStream << setw(4) << setprecision(1) /*<< scientific*/ << fabs(X_max-X_min)/((double)(conf->data_bins));
	}
	else
	{
		thisStream << setw(4) << setprecision(1) /*<< scientific*/ << fabs(X_max-X_min)/100.;
	}
	Y_ext.Append( thisStream.str() );
	Y_ext.Append(" ");
	TString Unit; Unit.Append( EdStyle::GetParamRootUnit( observableName ) );
	if( !StringProcessing::is_empty(Unit) )
	{
		Y_ext.Append(Unit);
		Y_ext.Append(" ");
	}
	Y_ext.Append(")");
	if( StringProcessing::is_empty( X_Title ) )
	{
		X_Title = EdStyle::GetParamRootName( observableName );
		TString unit = EdStyle::GetParamRootUnit( observableName );
		if( !StringProcessing::is_empty( unit ) )
		{
			X_Title.Append( " [" );
			X_Title.Append( unit );
			X_Title.Append( "]" );
		}
	}
	if( StringProcessing::is_empty( Y_Title ) )
	{
		Y_Title = "Candidates";
		Y_Title.Append( Y_ext );
	}
	input_data->GetYaxis()->SetRangeUser( Y_min, Y_max );
	input_data->GetYaxis()->SetTitle( Y_Title );
	input_data->GetYaxis()->SetTitleSize( (Float_t)YaxisTitleScale*input_data->GetYaxis()->GetTitleSize() );
	input_data->GetYaxis()->SetLabelSize( (Float_t)YaxisLabelScale*input_data->GetYaxis()->GetLabelSize() );
	input_data->GetXaxis()->SetRangeUser( X_min, X_max );
	input_data->GetXaxis()->SetTitle( X_Title );
	input_data->GetXaxis()->SetTitleSize( (Float_t)XaxisTitleScale*input_data->GetXaxis()->GetTitleSize() );
	input_data->GetXaxis()->SetLabelSize( (Float_t)XaxisLabelScale*input_data->GetXaxis()->GetLabelSize() );
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}
	c1->Modified();
	c1->Update();
	//input_data->GetXaxis()->CenterTitle( true );
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}
	c1->Update();
	TPaveText* myLatex=NULL;
	if( addLHCb )
	{
		myLatex = EdStyle::LHCbLabel();
	}
	if( addRightLHCb )
	{
		myLatex = EdStyle::RightLHCbLabel();
	}
	if( myLatex != NULL ) myLatex->Draw();
	TLegend* leg = NULL;
	if( conf->useLegend )
	{
		if( conf->TopRightLegend ) leg = EdStyle::LHCbLegend();
		else if( conf->TopLeftLegend ) leg = EdStyle::LHCbLeftLegend();
		else if( conf->BottomRightLegend ) leg = EdStyle::LHCbBottomLegend();
		else if( conf->BottomLeftLegend ) leg = EdStyle::LHCbBottomLeftLegend();
	}
	if( leg != NULL ) leg->SetTextSize( (Float_t) legend_size );
	if( leg != NULL ) leg->AddEntry( input_data, "Data", "pl" );
	unsigned int num=0;
	for( auto comp_i : input_components )
	{
		if( !Style_Key.empty() && num < Style_Key.size() )
				comp_i->SetLineStyle( (Style_t)Style_Key[num] );
		if( !Color_Key.empty() && num < Color_Key.size() )
				comp_i->SetLineColor( (Color_t)Color_Key[num] );
		if( !Width_Key.empty() && num < Width_Key.size() )
				comp_i->SetLineWidth( (Width_t)Width_Key[num] );
		if( comp_i->GetLineWidth() != 0 )
		{
			comp_i->Draw(drawSpline? "C" : "L");
			if( !component_names.empty() && leg != NULL )
			{
				if( num < component_names.size() )
					leg->AddEntry( comp_i, TString(component_names[num]), "l" );
				else
					leg->AddEntry( comp_i, "Unnamed", "l" );
			}
		}
		num++;
	}
	if( final_chi2 > 0 )
	{
		TString Chi2Text( "#chi^{2}/ndof : " );
		stringstream chi_stream; chi_stream << setprecision(4) << final_chi2;
		Chi2Text.Append( chi_stream.str() );
		if( leg !=NULL ) leg->AddEntry( (TObject*)NULL, Chi2Text, "" );
	}
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}
	c1->Update();
	if( !component_names.empty() ) if( leg != NULL ) leg->Draw();
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}
	if( !input_bin_theory_data.empty() )
	{
		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}
	c1->Update();
	if( !input_bin_theory_data.empty() )
	{
		pad2->cd();
		vector<double> pull_value, pull_error_value,  x_values,  x_errs;
		for( unsigned int i=0; i< input_bin_theory_data.size(); ++i )		//	Explicit Assumption that theory.size() == input_data->GetN()
		{
			double theory_y = input_bin_theory_data[i];
			double data_y = input_data->GetY()[i];
			double data_err = input_data->GetErrorY( i );
			double pull = ( data_y -theory_y ) / data_err;
			/*!
			 * From what I understand in comparing the PDF to data in this way the residuals are normalised with an error of 1.
			 *
			 * I don't implicitly trust this, however I don't have a better solution to calculate a sensible error at the time of writing
			 */
			double pull_err = 0;//1.;
			if( pull >= DBL_MAX || data_err < 1E-10 )	pull = 0.;
			if( fabs(data_y) < 1E-10 ) pull_err = 0.;
			pull_value.push_back( pull );
			pull_error_value.push_back( pull_err );
			x_values.push_back( input_data->GetX()[i] );
			x_errs.push_back( input_data->GetErrorX( i ) );
		}
		TGraphErrors* pullGraph = new TGraphErrors( pull_value.size(), x_values.data(), pull_value.data(), x_errs.data(), pull_error_value.data() );
		pullGraph->Draw("AB");
		// Style Y axis
		pullGraph->GetYaxis()->SetTitle( "Pull" );
		pullGraph->GetYaxis()->SetTitleSize( input_data->GetYaxis()->GetTitleSize() * 3 );
		pullGraph->GetYaxis()->SetLabelSize( input_data->GetYaxis()->GetLabelSize() * 2. );
		// Style X axis
		pullGraph->GetXaxis()->SetTitleSize( input_data->GetXaxis()->GetTitleSize() * 3. );
		pullGraph->GetXaxis()->SetLabelSize( input_data->GetXaxis()->GetLabelSize() * 3. );
		pullGraph->GetXaxis()->SetTitleOffset(pullGraph->GetXaxis()->GetTitleOffset() / 3.);
		pullGraph->GetXaxis()->SetRangeUser( X_min, X_max );
		// Move X-axis title and unit from main plot to pull plot
		pullGraph->GetXaxis()->SetTitle(input_data->GetXaxis()->GetTitle());
		input_data->GetXaxis()->SetLabelSize(0.);
		input_data->GetXaxis()->SetTitleSize(0.);
		pad1->SetBottomMargin(0.03);
		pad2->SetTopMargin(0.02);
		pad2->SetBottomMargin(0.50);
		// Pull limits
		if( limitPulls )
		{
			pullGraph->GetYaxis()->SetNdivisions( 3 );
			pullGraph->GetYaxis()->SetRangeUser( -5., 5. );
			pullGraph->SetMaximum( 5. );
			pullGraph->SetMinimum( -5. );
		}
		pad1->Modified();
		pad1->Update();
		pad2->Modified();
		pad2->Update();
		c1->Update();
	}
	c1->Write("",TObject::kOverwrite);
	TString Clean_Description = StringProcessing::Clean( CombinationDescription.c_str() );
	c1->SaveAs( TString("Overlay_"+observableName+"_"+Clean_Description+".C") );
	c1->SaveAs( TString("Overlay_"+observableName+"_"+Clean_Description+".pdf") );
	c1->SaveAs( TString("Overlay_"+observableName+"_"+Clean_Description+".png") );
}
//      This is a slimmed down version of the AddBranch function I have written elsewhere
void ComponentPlotter::WriteBranch( TTree* input_tree, TString Branch_Name, vector<double>& branch_data )
{
	if( input_tree->GetEntries() != 0 && branch_data.size() != (unsigned)input_tree->GetEntries() )
		return;
	Double_t X_object=-1.;
	TBranch* this_X_branch = input_tree->Branch( Branch_Name, &X_object, TString(Branch_Name+"/D") );
	input_tree->SetEntries( branch_data.size() );
	for(const auto& dat_x_i : branch_data )
	{
		X_object = dat_x_i;
		this_X_branch->Fill();
	}
}
//	Initialize the Desciptions of the Relevant Combinations in the DataSet
void ComponentPlotter::SetupCombinationDescriptions()
{
	cout << "Optimally:" << endl << endl;
	//      If we have 'no discrete combinations' then get the whole dataset for this pdf
	vector<DataPoint> whole_dataset;
	cout << plotData->GetDataNumber(NULL) << endl;
	for( int i=0; i< plotData->GetDataNumber(NULL); ++i )
	{
		whole_dataset.push_back( DataPoint(*plotData->GetDataPoint( i )) );
	}
	data_subsets.push_back( whole_dataset );
	ObservableRef ObservableName( observableName );
	vector<double> discrete_observableValues_i;
	//      Get all values of the wanted parameter from this subset
	for( unsigned int j=0; j< whole_dataset.size(); ++j )
	{
		discrete_observableValues_i.push_back( whole_dataset[j].GetObservable( ObservableName )->GetValue() );
	}
	cout << "Combination: 1 " << "Min: " << StatisticsFunctions::Minimum( discrete_observableValues_i );
	cout << "\tMax: " << StatisticsFunctions::Maximum( discrete_observableValues_i ) << "\tBinNum: " << StatisticsFunctions::OptimumBinNumber( discrete_observableValues_i ) << endl;
	combinationDescriptions.push_back( "All_Data" );
	discreteNames = full_boundary->GetDiscreteNames();
	//      If we have discrete combinations then we need to loop over each set seperatley to split up the data and information about the data
	if( allCombinations.size() > 1 )
	{
		combinationDescriptions.clear();
		//      For all discrete combinations that have been calculated elsewhere
		for( unsigned int combinationIndex=0; combinationIndex< allCombinations.size(); ++combinationIndex )
		{
			cout << "Combination: " << combinationIndex+2 << " ";
			//      For this combination work out what values correspond to the discrete parameters
			vector<double> values_for_this_combination;
			for( unsigned int paramIndex=0; paramIndex< discreteNames.size(); ++paramIndex )
				values_for_this_combination.push_back( allCombinations[combinationIndex].GetObservable( discreteNames[paramIndex] )->GetValue() );
			//              Get all data which falls within this paramreter set
			data_subsets.push_back( plotData->GetDiscreteSubSet( discreteNames, values_for_this_combination ) );
			vector<double> this_discrete_observableValues_i;
			//      Get all values of the wanted parameter from this subset
			for( unsigned int j=0; j< data_subsets.back().size(); ++j )
				this_discrete_observableValues_i.push_back( (data_subsets.back())[j].GetObservable( ObservableName )->GetValue() );
			cout << "Min: " << StatisticsFunctions::Minimum( this_discrete_observableValues_i );
			cout << "\tMax: " << StatisticsFunctions::Maximum( this_discrete_observableValues_i );
			cout << "\tBinNum: " << StatisticsFunctions::OptimumBinNumber( this_discrete_observableValues_i ) << endl;
			TString ThisName;
			for( unsigned int paramIndex=0; paramIndex< discreteNames.size(); ++paramIndex )
			{
				ThisName.Append( discreteNames[paramIndex] );
				ThisName.Append( "=" );
				ThisName+= values_for_this_combination[paramIndex];
				ThisName.Append( "_" );
			}
			cout << "Adding This: " << ThisName << endl;
			combinationDescriptions.push_back( ThisName.Data() );
		}
	}
	return;
}
//Return the values tracing the PDF projection in the Observable of choice
vector<double> ComponentPlotter::ProjectObservableComponent( DataPoint* InputPoint, string ObservableName, double Minimum, int PlotNumber, double PlotInterval, string component )
{
	//	Move initializer(s) outside of for loop
	double integralvalue=-1., observableValue=-1.;
	//Find the value of the observable projection at each data point
	vector<double> pointValues;
	//	This class object has been created to speed up the communication between this class and the basePDF as it may pass through several PDF wrappers
	ComponentRef comp_obj( component, ObservableName );
	string unit = InputPoint->GetObservable( string(ObservableName) )->GetUnit();
	ObservableRef thisObservableRef( ObservableName );
	//	Step over the whole observable range
	for (int pointIndex = 0; pointIndex < PlotNumber; ++pointIndex )
	{
		//	Inform the user of how far we have got :D
		//	Value of Observable we want to evaluate for this step
		observableValue = Minimum + ( PlotInterval * pointIndex );
		//	Set the value of Observable to the new step value
		//	All other CONTINUOUS observables set to average of range
		Observable newObs( string(ObservableName), observableValue, unit );
		DataPoint thisPoint( *InputPoint );
		thisPoint.SetObservable( thisObservableRef, &newObs );
		//	perform actual evaluation of the PDF with the configuration in the InputPoint, in the whole boundary with the given name and for selected component
		integralvalue = pdfIntegrator->ProjectObservable( &thisPoint, full_boundary, ObservableName, &comp_obj );
		pointValues.push_back( integralvalue );
	}
	Sanity_Check( pointValues, component );
	return pointValues;
}
//	Check that the passed vector of points aren't all NULL if they are print some debug information
void ComponentPlotter::Sanity_Check( vector<double>& pointValues, TString component )
{
	bool result = true;
	for( const auto& inf : pointValues )
	{
		if( fabs( inf ) > 1E-6 )
			result = false;
	}
	if( result )
	{
		cerr << "Failed to get anything other than 0 for PDF!!" << endl;
		full_boundary->Print();
		cerr << "component: " << component << endl;
		cerr << "\nCheck your Plots!" << endl;
	}
}
void ComponentPlotter::SetWeightsWereUsed( string input )
{
	double sum=0.;
	for( int i=0 ; i < plotData->GetDataNumber(NULL); ++i )
	{
		sum+=plotData->GetDataPoint( i )->GetObservable( input )->GetValue();
	}
	weightsWereUsed = true;
	weightName = input;
	ObservableRef weight_ref( weightName );
	for( int i=0 ; i < plotData->GetDataNumber(NULL); ++i )
		wanted_weights.push_back( plotData->GetDataPoint( i )->GetObservable( weight_ref )->GetValue() );
	weight_norm=0.;
	double weight_sum=plotData->GetSumWeights();
	weight_norm = weight_sum / plotData->GetDataNumber(NULL);
}
//	Functions to access the internal results in this class
vector<TGraphErrors*> ComponentPlotter::GetBinnedData()
{
	return binned_data;
}
vector<vector<TGraph*> > ComponentPlotter::GetComponents()
{
	return total_components;
}
double ComponentPlotter::PDF2DataNormalisation( const unsigned combinationIndex ) const // I swear this is obfuscating something much simpler. TODO Work out what it is.
{
	double normalisation=fabs(ratioOfIntegrals[ combinationIndex ]);			//	Attempt to correct for Numerical != analytical due to any constant factor due to numerical inaccuracy
	if(combinationIndex >= allCombinations.size())
	{
		cerr << "Can't get combination " << combinationIndex << " of " << allCombinations.size() << endl;
		exit(-1);
	}
	DataPoint PrototypeCombinationDatapoint = allCombinations[combinationIndex];
	double dataNum = 0.;
	if( allCombinations.empty() || allCombinations.size() == 1 ) dataNum = plotData->GetDataNumber();
	else dataNum = plotData->GetDataNumber( &PrototypeCombinationDatapoint );
	normalisation *= dataNum / (double) data_binning;				//	Normalise to this density of events	(Num of events per bin in flatPDF)
	double range = fabs( boundary_max-boundary_min );
	normalisation *= range;								//	Correct for the range of the dataset	(absolute range of Observable being projected)
	if( !std::isnan( combination_integral[ combinationIndex ] ) )
		normalisation /= combination_integral[ combinationIndex ];		//	Total Integral of the PDF	(We're plotting prob of PDF being at this point for a non-unitary PDF Evaluate)
	normalisation *= weight_norm;							//	Correct for the effect of non-unitary weights used in the fit
	return normalisation;
}
pair<double,double> ComponentPlotter::GetChi2Numbers()
{
	return make_pair( chi2, N );
}
vector<double> ComponentPlotter::GetFunctionEval()
{
	return allPullData;
}
TGraphErrors* ComponentPlotter::PullPlot1D( vector<double> input_bin_theory_data, TGraphErrors* input_data, string observableName, string CombinationDescription, TRandom* rand, CompPlotter_config* conf )
{
	if( rand == NULL ) rand = gRandom;
	vector<double> pull_value;
	vector<double> pull_error_value;
	vector<double> x_values;
	vector<double> x_errs;
	for( unsigned int i=0; i< input_bin_theory_data.size(); ++i )		//	Explicit Assumption that theory.size() == input_data->GetN()
	{
		double theory_y = input_bin_theory_data[i];
		double data_y = input_data->GetY()[i];
		double data_err = input_data->GetErrorY( i );
		double pull = ( data_y - theory_y ) / data_err;
		/*!
		 * From what I understand in comparing the PDF to data in this way the residuals are normalised with an error of 1.
		 *
		 * I don't implicitly trust this, however I don't have a better solution to calculate a sensible error at the time of writing
		 */
		double pull_err = 1.;
		if( pull >= DBL_MAX || data_err < 1E-10 )	pull = 0.;
		if( fabs(data_y) < 1E-10 ) pull_err = 0.;
		pull_value.push_back( pull );
		pull_error_value.push_back( pull_err );
		x_values.push_back( input_data->GetX()[i] );
		x_errs.push_back( input_data->GetErrorX( i ) );
	}
	TGraphErrors* pullGraph = new TGraphErrors( pull_value.size(), x_values.data(), pull_value.data(), x_errs.data(), pull_error_value.data() );
	TString pullGraphName="PullGraph_"; pullGraphName+=rand->Rndm();
	string pullGraphCleanName( pullGraphName.Data() );
	replace( pullGraphCleanName.begin(), pullGraphCleanName.end(), '.', '_' );
	pullGraph->SetName( pullGraphCleanName.c_str() ); pullGraph->SetTitle( pullGraphCleanName.c_str() );
	TString canvas_name = "PullPlot_"; canvas_name.Append(observableName); canvas_name.Append("_"); canvas_name+=rand->Rndm();
	string canvas_clean_name( canvas_name.Data() );
	replace( canvas_clean_name.begin(), canvas_clean_name.end(), '.', '_' );
	TCanvas* c1 = EdStyle::RapidFitCanvas( canvas_clean_name.c_str(), canvas_name );
	pullGraph->Draw("AP");
	c1->Update();
	if( conf != NULL )
	{
		if( conf->LimitPulls )
		{
			pullGraph->GetYaxis()->SetRangeUser(-5.,5.);
			pullGraph->SetMinimum(-5);
			pullGraph->SetMaximum(+5);
			c1->Update();
		}
	}
	TString X_Title;
	TString Y_Title;
	if( StringProcessing::is_empty( X_Title ) )
	{
		X_Title = EdStyle::GetParamRootName( observableName );
		TString unit = EdStyle::GetParamRootUnit( observableName );
		if( !StringProcessing::is_empty( unit ) )
		{
			X_Title.Append( " " );
			X_Title.Append( unit );
		}
	}
	if( StringProcessing::is_empty( Y_Title ) ) Y_Title = "Pull";
	input_data->GetYaxis()->SetTitle( Y_Title );
	input_data->GetXaxis()->SetTitle( X_Title );
	c1->Update();
	TString Clean_Description = StringProcessing::Clean( CombinationDescription.c_str() );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".C") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".pdf") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".png") );
	TCanvas* c2 = EdStyle::RapidFitCanvas( "pull", "pull");
	TString pull_name = "pull_"; pull_name.Append(observableName);
	TH1D * pull_histograms = new TH1D(pull_name, pull_name, StatisticsFunctions::OptimumBinNumber( pull_value ), -5, 5);
	for ( unsigned i = 0; i < input_bin_theory_data.size() ; i++ ) pull_histograms->Fill(pull_value[i]);
	pull_histograms->Draw();
	gStyle->SetOptFit(1);
	pull_histograms->Fit("gaus");
	c2->Update();
	c2->Print( TString("pull"+observableName+"_"+Clean_Description+".pdf") );
	delete pull_histograms;
	return pullGraph;
}
void ComponentPlotter::WriteData( TGraphErrors* Total_BinnedData, vector<TGraph*> Total_Components, TString destination )
{
	TTree* outputTree = new TTree( destination, destination );
	unsigned int points = (unsigned)Total_BinnedData->GetN();
	vector<double> X_val;
	vector<double> X_err;
	vector<double> Y_val;
	vector<double> Y_err;
	for( unsigned int i=0; i< points; ++i )
	{
		X_val.push_back( Total_BinnedData->GetX()[i] );
		X_err.push_back( Total_BinnedData->GetEX()[i] );
		Y_val.push_back( Total_BinnedData->GetY()[i] );
		Y_err.push_back( Total_BinnedData->GetEY()[i] );
	}
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Value_X", X_val  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Error_X", X_err  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Value_Y", Y_val  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Error_Y", Y_err  );
	for( unsigned int i=0; i< Total_Components.size(); ++i )
	{
		vector<double> this_X_val;
		vector<double> this_Y_val;
		unsigned int this_points = (unsigned)Total_Components[i]->GetN();
		for( unsigned int j=0; j< this_points; ++j )
		{
			this_X_val.push_back( Total_Components[i]->GetX()[j] );
			this_Y_val.push_back( Total_Components[i]->GetY()[j] );
			TString ComponentName="Component_";
			ComponentName+=i;
			TString ComponentXName=ComponentName;ComponentXName.Append("_X");
			TString ComponentYName=ComponentName;ComponentYName.Append("_Y");
			ComponentPlotter::WriteBranch( outputTree, ComponentXName, this_X_val );
			ComponentPlotter::WriteBranch( outputTree, ComponentYName, this_Y_val );
		}
	}
	outputTree->Write("",TObject::kOverwrite);
}
string ComponentPlotter::XML( int projectionType )
{
	stringstream xml;
	string projectionOuterTag="ComponentProjection";
	if( projectionType == 2 ) projectionOuterTag = "Projection";
	xml << "<" << projectionOuterTag << ">" << endl;
	xml << "\t" << "<Name>" << "someObservable" << "</Name>" << endl;
	xml << "\t" << "<CompNames>" << "component1Name:component2Name:component3Name:..." << "</CompNames> # or <ComponentNames> refers to the PDF Components eg CP-Even:CP-Odd:Swave" << endl;
	xml << "\t" << "<CombinationNames>" << "combination1Name:combination2Name:combination3Name:..." << "</CombinationNames> # refers to the different Discrete Combinations eg tag=+1:tag=-1:untagged" << endl;
	xml << "\t" << "<Xmax>" << "XaxisMax" << "</Xmax>" << endl;
	xml << "\t" << "<Xmin>" << "XaxisMin" << "</Xmin>" << endl;
	xml << "\t" << "<Ymax>" << "YaxisMax" << "</Ymax>" << endl;
	xml << "\t" << "<Ymin>" << "YaxisMin" << "</Ymin>" << endl;
	xml << "\t" << "<XTitle>" << "XaxisTitle" << "</XTitle>" << endl;
	xml << "\t" << "<YTitle>" << "YaxisTitle" << "</YTitle>" << endl;
	xml << "\t" << "<TrustNumerical>" << "True/False" << "</TrustNumerical> # Correct for Numerical vs analytical difference" << endl;
	xml << "\t" << "<CalcChi2>" << "True/False" << "</CalcChi2> # Calculate Chi2 for the Projection" << endl;
	xml << "\t" << "<DrawPull>" << "True/False" << "</DrawPull> # Draw Residual Plot" << endl;
	xml << "\t" << "<LimitPulls>" << "True/False" << "</LimitPulls> # limit pull plot to +/- 5" << endl;
	xml << "\t" << "<AddLHCb>" << "True/False" << "</AddLHCb> # Add LHCb logo on top left" << endl;
	xml << "\t" << "<AddRightLHCb>" << "True/False" << "</AddRightLHCb> # Add LHCb logo on top right" << endl;
	xml << "\t" << "<LegendTextSize>" << "Size_t" << "</LegendTextSize> # Legend font size default: 0.05" << endl;
	xml << "\t" << "<TopRightLegend>" << "True/False" << "</TopRightLegend> # Add Legend on top right of plot" << endl;
	xml << "\t" << "<TopLeftLegend>" << "True/False" << "</TopLeftLegend> # Add Legend on top left of plot" << endl;
	xml << "\t" << "<BottomRightLegend>" << "True/False" << "</BottomRightLegend> # Add Legend on bottom right of plot" << endl;
	xml << "\t" << "<BottomLeftLegend>" << "True/False" << "</BottomLeftLegend> # Add Legend on bottom left of plot" << endl;
	xml << "\t" << "<NoLegend>" << "True/False" << "</NoLegend> # Do Not draw the Legend" << endl;
	xml << "\t" << "<UseSpline>" << "True/False" << "</UseSpline> # Use a Spline to interpolate between PDF points" << endl;
	xml << "\t" << "<Threads>" << "NumberOfThreads" << "</Threads> # Number of Threads that GSL will use when multi-threading" << endl;
	xml << "\t" << "<FixedIntegrationPoints>" << "numberOfPointsPerGSLIntegral" << "</FixedIntegrationPoints> # Set the Number of points for GSL to use to be non default" << endl;
	xml << "\t" << "<UseGSLNumericalIntegration>" << "True/False" << "</UseGSLNumericalIntegration> # Use the GSL Integrator for projections, this is multi-threaded so better" << endl;
	xml << "\t" << "<StyleKey>" << "LineStyle1:LineStyle2:LineStyle3:..." << "</StyleKey> # Styles to use for different lines" << endl;
	xml << "\t" << "<ColorKey>" << "LineColor1:LineColor2:LineColor3:..." << "</ColorKey> # Colors to use for different lines" << endl;
	xml << "\t" << "<WidthKey>" << "LineWidth1:LineWidth2:LineWidth3:..." << "</WidthKey> # Widths of lines on plot, width 0 means a line is not drawn or added to Legend" << endl;
	xml << "\t" << "<LogY>" << "True/False" << "</LogY> # Set Log on Y axis" << endl;
	xml << "\t" << "<LogX>" << "True/False" << "</LogX> # Set Log on X axis" << endl;
	xml << "\t" << "<PDFpoints>" << "numberOfPDFIntegralPoints" << "</PDFpoints> # Number of bins that the PDF will be evaluated(integrated) at in Xaxis" << endl;
	xml << "\t" << "<DataBins>" << "numberOfBinsInDataHisto" << "</DataBins> # Number of bins that will be used in constructing the histogram of data" << endl;
	xml << "\t" << "<PlotAllCombinatons>" << "True/False" << "</PlotAllCombinatons> # Do we need to project tag+/-1 separately?" << endl;
	xml << "\t" << "<DefaultDiscreteValue>" << "someValue" << "</DefaultDiscreteValue> # Advanced usage to replace a series of discrete values with a non-zero default" << endl;
	xml << "\t" << "<XaxisTitleScale>" << "someScale" << "</XaxisTitleScale> # Scale the Font Size of the X axis Title" << endl;
	xml << "\t" << "<XaxisLabelScale>" << "someScale" << "</XaxisLabelScale> # Scale the Font Size of the X axis Labels" << endl;
	xml << "\t" << "<YaxisTitleScale>" << "someScale" << "</YaxisTitleScale> # Scale the Font Size of the Y axis Title" << endl;
	xml << "\t" << "<YaxisLabelScale>" << "someScale" << "</YaxisLabelScale> # Scale the Font Size of the Y axis Labels" << endl;
	xml << "</" << projectionOuterTag << ">" << endl;
	return xml.str();
}

