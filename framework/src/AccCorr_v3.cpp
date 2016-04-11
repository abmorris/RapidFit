#include "AccCorr_v3.h"
#include "TRandom3.h"
#include <sstream>
#include <fstream>

AccCorr_v3::AccCorr_v3():
AccLoc(3)
,histo()
//,xaxis()
//,yaxis()
//,zaxis()
//,histo1D()
,Dim()
,nxbins(0),nybins(0),nzbins(0),xmin(0),xmax(0),ymin(0),
 ymax(0),zmin(0),zmax(0),deltax(0),deltay(0),deltaz(0),total_num_entries(0)
,angNorm(0.),w1(0.0),w2(0.0),w3(0.0),w4(0.0),w5(0.0),w6(0.0)
,w7(0.0),w8(0.0),w9(0.0),w10(0.0),w11(0.0),w12(0.0)
,w13(0.0),w14(0.0),w15(0.0)
{
}

AccCorr_v3::~AccCorr_v3(){
}

AccCorr_v3::AccCorr_v3( const AccCorr_v3& ACclass): 
w1(ACclass.w1),w2(ACclass.w2),w3(ACclass.w3),w4(ACclass.w4),w5(ACclass.w5),w6(ACclass.w6)
,w7(ACclass.w7),w8(ACclass.w8),w9(ACclass.w9),w10(ACclass.w10),w11(ACclass.w11),w12(ACclass.w12)
,w13(ACclass.w13),w14(ACclass.w14),w15(ACclass.w15)
,histo(ACclass.histo)
//,histo1D(ACclass.histo1D)
//,xaxis(ACclass.xaxis)
//,yaxis(ACclass.yaxis)
//,zaxis(ACclass.zaxis)
,Dim(ACclass.Dim)
,AccLoc(ACclass.AccLoc),nxbins(ACclass.nxbins),nybins(ACclass.nybins),nzbins(ACclass.nzbins),xmin(ACclass.xmin),xmax(ACclass.xmax)
,ymin(ACclass.ymin),ymax(ACclass.ymax),zmin(ACclass.zmin),zmax(ACclass.zmax),deltax(ACclass.deltax),deltay(ACclass.deltay),deltaz(ACclass.deltaz)
,total_num_entries(ACclass.total_num_entries),angNorm(ACclass.angNorm)
{
}

void AccCorr_v3::AccCorr_v3Init(PDFConfigurator* config) {
	cout << "Initialising acceptance class" << endl;	
	// Determine where we get our weights from: MC (1), HIST (2), FIT (3), File (4)
	string orig = config->getConfigurationValue( "AccOrig" ) ;
	Dim = config->getConfigurationValue( "AccDim" ) ;
	if(orig=="MC") AccLoc = 1 ;
	if(orig=="HIST") AccLoc = 2 ;
	if(orig=="FIT") AccLoc = 3 ;
	if(orig=="FILE") AccLoc = 4 ;
	string fileName;	
	if (AccLoc == 2) {
		cout << "   Acceptance histogram requested: " << fileName << endl ;
		//Find name of histogram needed to define 3-D angular acceptance distribution
		fileName = config->getConfigurationValue( "AccLocation" ) ;

		//File location
		ifstream input_file;
		input_file.open( fileName.c_str(), ifstream::in );
		input_file.close();
		
		if( !getenv("RAPIDFITROOT") && input_file.fail() )
		{
			cerr << "\n\n" << endl;
			cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
			cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "\n\n" << endl;
			exit(-987);
		}
		string fullFileName;
		
		if( getenv("RAPIDFITROOT") )
		{
			string path( getenv("RAPIDFITROOT") ) ;
			
			cout << "RAPIDFITROOT defined as: " << path << endl;
			
			fullFileName = path+"/pdfs/configdata/"+fileName ;
		} else if( !input_file.fail() )
		{
			fullFileName = fileName;
		} else {
			cerr << "Shouldn't end up Here in the code!" << endl;
			exit(-892);
		}
		//Read in histo
		cout << "About to open file" << endl;
		TFile* f = TFile::Open(fullFileName.c_str());
		cout << "Found file successfully: " << fullFileName.c_str() << endl;	
		//if(Dim=="3"){
		
		histo = new TH3F(  *(    (TH3F*)f->Get("fit3d")      )     ); //(fileName.c_str())));
		// *********************************** CONSISTENCY CHECK ********************************
		// xaxis
		TAxis* xaxis = histo->GetXaxis();
        	xmin = xaxis->GetXmin();
        	xmax = xaxis->GetXmax();
        	nxbins = histo->GetNbinsX();
        	deltax = (xmax-xmin)/nxbins;
		// yaxis
        	TAxis* yaxis = histo->GetYaxis();
        	ymin = yaxis->GetXmin();
        	ymax = yaxis->GetXmax();
        	nybins = histo->GetNbinsY();
        	deltay = (ymax-ymin)/nybins;
		// zaxis
        	TAxis* zaxis = histo->GetZaxis();
        	zmin = zaxis->GetXmin();
        	zmax = zaxis->GetXmax();
       	 	nzbins = histo->GetNbinsZ();
        	deltaz = (zmax-zmin)/nzbins;

		//method for Checking whether histogram is sensible
		total_num_entries = histo->GetEntries();
		int total_num_bins = nxbins * nybins * nzbins;
		double sum = 0.;
		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; i++){
			for (int j=1; j < nybins+1; j++){
				for (int k=1; k < nzbins+1; k++){

					double bin_content = histo->GetBinContent(i,j,k);
					//cout << "Bin content: " << bin_content << endl;
					if(bin_content<=0) { zero_bins.push_back(1);}
					//cout << " Zero bins " << zero_bins.size() << endl;}
					else if (bin_content>0){
					sum +=  bin_content;}
				}
			}
		}
		
		double average_bin_content = sum / total_num_bins;

		cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
		cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
		cout << endl;
		cout << endl;

		if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*TMath::Pi() ){
			cout << "In TimeIntegratedBkg_3Dangular_Helicity::TimeIntegratedBkg_3Dangular_Helicity: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
		// **************************************** END OF CONSISTENCY CHECK *************************************

		// Normalise histogram
		double scale=histoArea(histo)/(histo->GetXaxis()->GetNbins()*histo->GetYaxis()->GetNbins()*histo->GetZaxis()->GetNbins()-zero_bins.size());
		histo->Scale(1./scale);
		/*
		}
		else {
		
			cout << "1D histogram acceptance requested" << endl;
			histo1D = new TH1F(  *(    (TH1F*)f->Get("accep1d")      )     );
			cout << "Histogram loaded successfully" << endl;
			double area=0.0;
			for(int n=1; n<histo1D->GetEntries()+1; n++) {
				area+=histo1D->GetBinContent(n);
			}
			histo1D->Scale(histo1D->GetNbinsX()/area);
		//}
		//f.Close();
		*/
	}
	if(AccLoc==4) {
		cout << "Getting weights from file" << endl;
		string AccFileName = config->getConfigurationValue( "AccLocation" ) ;
		cout << "Looking in file : " << AccFileName << endl;
		WeightsFromFile(AccFileName);
		cout << "Weights are: " << endl;
		cout << "w1 = " << w1 << endl;
		cout << "w2 = " << w2 << endl;
		cout << "w3 = " << w3 << endl;
		cout << "w4 = " << w4 << endl;
		cout << "w5 = " << w5 << endl;
		cout << "w6 = " << w6 << endl;

		cout << "w7 = " << w7 << endl;
		cout << "w8 = " << w8 << endl;
		cout << "w9 = " << w9 << endl;
		cout << "w10 = " << w10 << endl;
		cout << "w11 = " << w11 << endl;
		cout << "w12 = " << w12 << endl;
		cout << "w13 = " << w13 << endl;
		cout << "w14 = " << w14 << endl;
		cout << "w15 = " << w15 << endl;

		return;
	}
	MakeWeights(config);
	return;
}

double AccCorr_v3::AccEval(vector<double> x) {
	double factor= 0.0;	
	if(AccLoc==3){
		double y[3];
		y[0]=x[0];y[1]=x[1];y[2]=x[2];
		Bs2PhiPhi_AccCorr FIT = Bs2PhiPhi_AccCorr();
		factor = FIT.MDF(y);
		//cout << "factor =" << factor << endl;
	}
	else if(AccLoc==2) {
		// Find the bin
		//int globalbin = histo->FindBin(x[0], x[1], x[2]);
		//double num_entries_bin = histo->GetBinContent(globalbin);
		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
    		factor = CleverInterpolate(histo, x[0], x[1], x[2]);
    		cout << "factor = " << factor << endl;
		//if(Dim=="3") factor = CleverInterpolate(histo, x[0], x[1], x[2]);
		//else if(Dim=="1") {
			//cout << "Dim = " << Dim << endl;
			//factor = histo1D->GetBinContent(histo1D->FindBin(x[0]));
			//factor *= histo1D->GetBinContent(histo1D->FindBin(x[1]));
		//cout << "factor = " << factor << endl;
		//factor = num_entries_bin ;
		//}
	}
	else {
		cout << "NO LOCATION PROVIDED" << endl;
		return 0.0;
	}
	//cout << "factor (END) =" << factor << endl;
	return factor;
}

void AccCorr_v3::MakeWeights(PDFConfigurator* config) {
	// DECIDE IF WE GET THE WEIGHTS FROM A FILE OR NOT
	stringstream fname;
	fname << "Weights_" << config->getConfigurationValue( "DATATYPE" ) << "_Iter_" << config->getConfigurationValue( "IterNum" ) << ".txt";
	if( config->getConfigurationValue( "MAKE-KKRW" ) != "TRUE"){
		WeightsFromFile(fname.str());
		return;
	}
	
	// Get the sub-XML from configdata
	vector<pair<string, string> >* XMLOverrideList = new vector<pair<string,string> >;	
	XMLConfigReader* xmlFileACC = new XMLConfigReader( config->getConfigurationValue( "ACCEPXMLname" ), XMLOverrideList);

	// Extra code to use weights
	FitFunctionConfiguration * theFunction=NULL;
	OutputConfiguration * makeOutput = NULL;
	makeOutput = xmlFileACC->GetOutputConfiguration();
	theFunction = xmlFileACC->GetFitFunctionConfiguration();

	// If weights were specified then we need to let the output plotting know
	if( theFunction->GetWeightsWereUsed() ) makeOutput->SetWeightsWereUsed( theFunction->GetWeightName() ) ;
	
	// Get the PDF and dataset from our sub-XML
	PDFWithData * pdfAndData = xmlFileACC->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( xmlFileACC->GetFitParameters() );
	IDataSet * dataSet = pdfAndData->GetDataSet();
	IPDF * PDF = pdfAndData->GetPDF();	
	PhaseSpaceBoundary * boundary = dataSet->GetBoundary();
	
	// Initialise our variables
	const int numAngularTerms = 15;
	double*  f = new double[numAngularTerms]; // the angular functions
	//f[15]={0.0};
	double xi[numAngularTerms]; // the angular weights
	double ct1, phi, ct2, time; (void) time; double weightval=0.0;
	double evalPDFraw, evalPDFnorm, val;
	int numEvents = dataSet->GetDataNumber();
	for (int i = 0; i < numAngularTerms; i++) xi[i] = 0.0;

	double Sum[numAngularTerms];
	double Sum_sq[numAngularTerms][numAngularTerms];
	double cov[numAngularTerms][numAngularTerms];
	for (int i = 0; i < numAngularTerms; i++) {
		Sum[i] = 0.0;
		for (int k = 0; k < numAngularTerms; k++) {
			cov[i][k] = 0.0;
			Sum_sq[i][k] = 0.0;
		}
	}

	// ADDED TO MAKE AN EVENT WEIGHT FOR THE ITERATIVE RE-WEIGHTING
	double weightRW,BPT,weightKK,kk1pt,kk2pt,weightHel,weightBPT;
	float chi,ann,bdt;
	XMLConfigReader* xmlFileACC_DATA = new XMLConfigReader( config->getConfigurationValue( "DATAXML_PARAM" ), XMLOverrideList);
	PDFWithData * pdfAndData_DATA = xmlFileACC_DATA->GetPDFsAndData()[0];
	//
	pdfAndData_DATA->SetPhysicsParameters( xmlFileACC_DATA->GetFitParameters() );
	IDataSet * dataSet_DATA = pdfAndData_DATA->GetDataSet();
	IPDF * PDF_DATA = pdfAndData_DATA->GetPDF();
	//
	string it_num = config->getConfigurationValue( "IterNum" );
	string data_type = config->getConfigurationValue( "DATATYPE" );
	//
	cout << "Creating phi p_T vs. phi p_T histogram" << endl;
	TH2F histKK;
	if(config->getConfigurationValue( "MAKE-KKRW" )=="TRUE") histKK = makeRWHistKK(it_num,data_type,config,dataSet_DATA,PDF_DATA,dataSet,PDF);
	else {
		stringstream ss_in;
		ss_in << "phi_vs_phi_Iter_" << it_num << "_" << data_type << ".root";
		TFile f_in(ss_in.str().c_str());
		histKK = *(TH2F*)f_in.Get("Kpt2D_ratio");
	}
	//
	cout << "*************** ITERATIVE RE-WEIGHTING *****************" << endl;
	cout << "ITERATION: " << it_num <<  endl;
	cout << "Taking data parameters from: " << config->getConfigurationValue( "DATAXML_PARAM" ) << endl;
	cout << "Taking MC parameters from: " << config->getConfigurationValue( "ACCEPXMLname" ) << endl;
	//

	int seedgen=1951;
	TRandom3* rnd = new TRandom3();
	rnd->SetSeed(seedgen);
	cout << "Starting event loop..." << endl;
	for (int e = 0; e < numEvents; e++)
	{
		if (e % 10000 == 0) cout << "Event # " << e << endl;
		
		// Get the observables for the event
		DataPoint * event = dataSet->GetDataPoint(e);
		ct1 = event->GetObservable("ctheta_1")->GetValue();
		phi      = event->GetObservable("phi")->GetValue();
		ct2   = event->GetObservable("ctheta_2")->GetValue();
		time     = event->GetObservable("time")->GetValue();
		
		// ADDED TO MAKE AN EVENT WEIGHT FOR THE ITERATIVE RE-WEIGHTING
		BPT = event->GetObservable("B_s0_LOKI_B_s0_PT")->GetValue();
		kk1pt = event->GetObservable("B_s0_LOKI_phi1_PT")->GetValue();
		kk2pt = event->GetObservable("B_s0_LOKI_phi2_PT")->GetValue();
		chi = event->GetObservable("maxTRACK_CHI2NDOF_MVA")->GetValue();
		ann = event->GetObservable("minProbK_MVA")->GetValue();
		bdt = event->GetObservable("bdt")->GetValue();
		weightKK = histKK.GetBinContent( histKK.FindBin( kk1pt , kk2pt ) );
		weightBPT = getBWeight(data_type,BPT,chi,ann,bdt);
		weightHel = getHelicityWeight(config,PDF_DATA,PDF,event);
		//weightHel = 1.0;
		weightRW = weightHel * weightBPT * weightKK;
		//
		if (e % 10000 == 0) {
			cout << "Weights I am using are:" << endl;
			cout << "\t" << "BPT: " << weightBPT << endl;
			cout << "\t" << "Helicity: " << weightHel << endl;
			cout << "\t" << "phi vs. phi: " << weightKK << endl;
		}

		// ONLY PHI PHI SPECIFIC PART
		Mathematics::getBs2PhiPhiAngularFunctions(f[0], f[1], f[2], f[3], f[4], f[5], 
		f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13], f[14], ct1, ct2, phi);
		// The method sums the angular factor f_i divided by the sum_i(A_i*f_i)
		// for each accepted event. I'm not sure if dividing by Evaluate is exactly the
		// same here, particularly if you look at untagged events.
		evalPDFraw = PDF_DATA->Evaluate( event );
		// Now need to calculate the normalisation when integrated over the 3 angles
		vector<string> dontIntegrate = PDF_DATA->GetDoNotIntegrateList();
		dontIntegrate.push_back("time");
		evalPDFnorm = PDF_DATA->EvaluateTimeOnly(event);
		//evalPDFnorm = PDF->Integral(event,boundary); // Line relevant when not time dependent
	
		// Code to make the weights
		// Calculate val based on method of weight calculation
		val = evalPDFraw/evalPDFnorm/weightRW;
		//}
		/*	
		else if(AccLoc==2) {
			double passval = rnd->Uniform(0.,1.);
			if(CleverInterpolate(histo,ct1,ct2,phi)>passval) {
			//if(histo->GetBinContent(histo->FindBin(ct1,ct2,phi))>passval) {
				if( (theFunction->GetWeightsWereUsed()) && (weightval > 0.0)) {
					val = evalPDFraw/evalPDFnorm*weightval;
				}
				else {
					val = evalPDFraw/evalPDFnorm;
				}
		
			}
		}
		else if(AccLoc==3) {
			double passval = rnd->Uniform(0.,1.);
			double y[3];
			y[0]=ct1;y[1]=ct2;y[2]=phi;
			Bs2PhiPhi_AccCorr FIT = Bs2PhiPhi_AccCorr();
			if(FIT.MDF(y)>passval) {
				if( (theFunction->GetWeightsWereUsed()) && (weightval > 0.0))
					val = evalPDFraw/evalPDFnorm*weightval;
				else {
					val = evalPDFraw/evalPDFnorm;
				}
		
			}
		}
		*/
		for (int i = 0; i < numAngularTerms; i++)
		{
			Sum[i] = Sum[i] + f[i]/val;
			xi[i] += f[i]/val;
			for (int k = 0; k < numAngularTerms; k++)
			{
			Sum_sq[i][k] += f[i]/val*f[k]/val;
			}
		}
	}
	delete rnd;
	double AvNonZero = (xi[0]+xi[1]+xi[2]+xi[6]+xi[7])/5./numEvents;
	// Print out the results and assign the weights we will use in PDF
	for (int i = 0; i < numAngularTerms; i++)
	{
		for (int k = 0; k < numAngularTerms; k++)
		{
			cov[i][k] = 1./numEvents/numEvents * ( Sum_sq[i][k] - Sum[i]*Sum[k]/numEvents);
			cout << cov[i][k]/AvNonZero/AvNonZero << "\t";
		}
		cout << endl;
	}
	cout << "Weights:" << endl;
	for (int i = 0; i < numAngularTerms; i++)
	{
			cout << xi[i]/numEvents/AvNonZero << " \\pm " << sqrt(cov[i][i])/AvNonZero << endl;
	}	
	cout << "AvNonZero = " << AvNonZero << endl;	
	w1 = xi[0]/numEvents/AvNonZero;
	w2 = xi[1]/numEvents/AvNonZero;
	w3 = xi[2]/numEvents/AvNonZero;
	w4 = xi[3]/numEvents/AvNonZero;
	w5 = xi[4]/numEvents/AvNonZero;
	w6 = xi[5]/numEvents/AvNonZero;
	
	w7 = xi[6]/numEvents/AvNonZero;
	w8 = xi[7]/numEvents/AvNonZero;
	w9 = xi[8]/numEvents/AvNonZero;
	w10 = xi[9]/numEvents/AvNonZero;
	w11 = xi[10]/numEvents/AvNonZero;
	w12 = xi[11]/numEvents/AvNonZero;
	w13 = xi[12]/numEvents/AvNonZero;
	w14 = xi[13]/numEvents/AvNonZero;
	w15 = xi[14]/numEvents/AvNonZero;

	// WRITE WEIGHTS TO FILE;
	ofstream file;
	file.open(fname.str().c_str());
	file << w1 << endl;
	file << w2 << endl;
	file << w3 << endl;
	file << w4 << endl;
	file << w5 << endl;
	file << w6 << endl;
	file << w7 << endl;
	file << w8 << endl;
	file << w9 << endl;
	file << w10 << endl;
	file << w11 << endl;
	file << w12 << endl;
	file << w13 << endl;
	file << w14 << endl;
	file << w15 << endl;
	file.close();

	return;	

}

vector<double> AccCorr_v3::ReturnWeights() {
	vector<double> retVect;
	retVect.push_back(w1);
	retVect.push_back(w2);
	retVect.push_back(w3);
	retVect.push_back(w4);
	retVect.push_back(w5);
	retVect.push_back(w6);

	retVect.push_back(w7);
	retVect.push_back(w8);
	retVect.push_back(w9);
	retVect.push_back(w10);
	retVect.push_back(w11);
	retVect.push_back(w12);
	retVect.push_back(w13);
	retVect.push_back(w14);
	retVect.push_back(w15);

	return retVect;
}

void AccCorr_v3::WeightsFromFile(string fname) {
	ifstream file;
	file.open(fname.c_str());
	string line;
	int n=0;
	cout << "Requested weights from file: " << fname.c_str() << endl;
	cout << "Weights found:" << endl;

	while(getline(file, line))    //read stream line by line
	{
		stringstream in(line);
		cout << "in = " << in.str().c_str() << endl;
		if(n==0) in >> w1;
		if(n==1) in >> w2;
		if(n==2) in >> w3;
		if(n==3) in >> w4;
		if(n==4) in >> w5;
		if(n==5) in >> w6;
		if(n==6) in >> w7;

		if(n==7) in >> w8;
		if(n==8) in >> w9;
		if(n==9) in >> w10;
		if(n==10) in >> w11;
		if(n==11) in >> w12;
		if(n==12) in >> w13;
		if(n==13) in >> w14;
		if(n==14) in >> w15;

		n++;
	}
	file.close();
}

double AccCorr_v3::CleverInterpolate(TH3F* hist, double x, double y, double z) {
	double ret=0.;

	if( (x<hist->GetBinCenter(1)) || (x>hist->GetBinCenter(hist->GetXaxis()->GetNbins())) ||
			(y<hist->GetBinCenter(1)) || (y>hist->GetBinCenter(hist->GetYaxis()->GetNbins())) ||
			(z<hist->GetBinCenter(1)) || (z>hist->GetBinCenter(hist->GetZaxis()->GetNbins())) ) {
	
		ret = hist->GetBinContent(hist->FindBin(x,y,z));
	}
	else {
		ret = hist->Interpolate(x,y,z);
	}

	return ret;
}

double AccCorr_v3::histoArea(TH3F* hist) {
	double ret=0.0;

	for(int nx=1;nx<hist->GetXaxis()->GetNbins()+1;nx++) {
		for(int ny=1;ny<hist->GetYaxis()->GetNbins()+1;ny++) {
			for(int nz=1;nz<hist->GetZaxis()->GetNbins()+1;nz++) {
				ret+=hist->GetBinContent(nx,ny,nz);
			}
		}
	}
	return ret;
}

double AccCorr_v3::histoArea(TH2F* hist) {
	double ret=0.0;

	for(int nx=1;nx<hist->GetXaxis()->GetNbins()+1;nx++) {
		for(int ny=1;ny<hist->GetYaxis()->GetNbins()+1;ny++) {
			ret+=hist->GetBinContent(nx,ny);
		}
	}
	return ret;
}

double AccCorr_v3::getHelicityWeight(PDFConfigurator* config,IPDF* PDF, IPDF* PDF_GEN,DataPoint* event){

	double ret=0.0;
/*
	vector<pair<string, string> >* XMLOverrideList = new vector<pair<string,string> >;

	XMLConfigReader* xmlFileACC = new XMLConfigReader( config->getConfigurationValue( "DATAXML_PARAM" ), XMLOverrideList);
	XMLConfigReader* xmlFileACC_GEN = new XMLConfigReader( config->getConfigurationValue( "ACCEPXMLname" ), XMLOverrideList);
	
	// Get the PDF and parameters from the data XML
	PDFWithData * pdfAndData = xmlFileACC->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( xmlFileACC->GetFitParameters() );
	IPDF * PDF = pdfAndData->GetPDF();	
	// Get the PDF and parameters from the MC XML
	PDFWithData * pdfAndData_GEN = xmlFileACC_GEN->GetPDFsAndData()[0];
	pdfAndData_GEN->SetPhysicsParameters( xmlFileACC_GEN->GetFitParameters() );
	IPDF * PDF_GEN = pdfAndData_GEN->GetPDF();	

	double num = PDF->Evaluate( event )/PDF->EvaluateTimeOnly(event);
	double den = PDF_GEN->Evaluate( event )/PDF_GEN->EvaluateTimeOnly(event);
*/
	double num = PDF->Evaluate( event );
	double den = PDF_GEN->Evaluate( event );
	ret = num/den;
	return ret;
}

double AccCorr_v3::getBWeight(string type, double BPT, double chi, double ann, double bdt){
	double ret=1.0;
	
	stringstream ss_type;
	ss_type << "/afs/cern.ch/work/s/sbenson/Rapidfit_Current/trunk/pdfs/configdata/BPTrw_" << type << ".root";
	// Extra weights from review
	stringstream ss_sel;
	//ss_sel << "/afs/cern.ch/work/s/sbenson/Rapidfit_Current/trunk/pdfs/configdata/Sel_Weights.root";
	ss_sel << "/afs/cern.ch/work/s/sbenson/Rapidfit_Current/trunk/pdfs/configdata/BDTrw_" << type << ".root";
	double bdt_min,bdt_max;
	if(type=="2011"){
		bdt_min=0.14;
		bdt_max=0.5;
	} else {
		bdt_min=0.06;
		bdt_max=0.4;
	}

	//TFile f_B(ss_type.str().c_str());
	//TH1F* h_B = (TH1F*)f_B.Get("BPTratio");
	//if(BPT>1000.0&&BPT<40000) ret = h_B->GetBinContent( h_B->FindBin(BPT) );
	//else ret=1.0;
	// Multiply by extra weights
	TFile f_sel(ss_sel.str().c_str());
	//TH1F* h_chi = (TH1F*)f_sel.Get("r_maxKchi2");
	//TH1F* h_ann = (TH1F*)f_sel.Get("r_minANN");
	TH1F* h_bdt = (TH1F*)f_sel.Get("BDTratio");
	//if(chi>0.8&&chi<3.0) ret *= h_chi->GetBinContent( h_chi->FindBin(chi) );
	//if(ann>0.0&&ann<1.0) ret *= h_ann->GetBinContent( h_ann->FindBin(ann) );
	if(bdt>bdt_min&&bdt<bdt_max) ret *= h_bdt->GetBinContent( h_bdt->FindBin(bdt) );

	return ret;
}

TH2F AccCorr_v3::makeRWHistKK(string it_num, string type, PDFConfigurator* config, IDataSet* dset_num, IPDF* data_pdf, IDataSet* dset_den, IPDF* pdf_mc){
	
	double binninglow[]={0.0,2300.0,50000.0};
	double binninghigh[]={0.0,4600.0,50000.0};

	int numEvents = dset_den->GetDataNumber();	
	int numEvents_num = dset_num->GetDataNumber();	
	//
	double kk1pt,kk2pt,bpt,bpt_weight,hel_weight;
	double kk1pt_d,kk2pt_d,sw;
	float chi,ann,bdt;

	TH2F h_2D_KKpt_d("Kpt2D_d","Kpt2D_d",2,binninglow,2,binninghigh);
	TH2F h_2D_KKpt("Kpt2D","Kpt2D",2,binninglow,2,binninghigh);
	TH2F h_2D_KKpt_ratio("Kpt2D_ratio","Kpt2D_ratio",2,binninglow,2,binninghigh);

	// Fill for data
	cout << "Filling data histogram" << endl;
	for(int i=0;i<numEvents_num;i++){
		// Get the observables for the event
		DataPoint * event_d = dset_num->GetDataPoint(i);
		kk1pt_d = event_d->GetObservable("B_s0_LOKI_phi1_PT")->GetValue();
		kk2pt_d = event_d->GetObservable("B_s0_LOKI_phi2_PT")->GetValue();
		sw = event_d->GetObservable("Nsig_sw_PB_LbLoKi")->GetValue();
		h_2D_KKpt_d.Fill(kk1pt_d,kk2pt_d,sw);
	}

	// Fill for MC
	cout << "Filling MC histogram" << endl;
	for(int i=0;i<numEvents;i++){
		// Get the observables for the event
		DataPoint * event = dset_den->GetDataPoint(i);
		kk1pt = event->GetObservable("B_s0_LOKI_phi1_PT")->GetValue();
		kk2pt = event->GetObservable("B_s0_LOKI_phi2_PT")->GetValue();
		bpt = event->GetObservable("B_s0_LOKI_B_s0_PT")->GetValue();
		chi = event->GetObservable("maxTRACK_CHI2NDOF_MVA")->GetValue();
		ann = event->GetObservable("minProbK_MVA")->GetValue();
		bdt = event->GetObservable("bdt")->GetValue();
		//
		hel_weight = getHelicityWeight(config,data_pdf,pdf_mc,event);
		//hel_weight = 1.0;
		bpt_weight = getBWeight(type,bpt,chi,ann,bdt);
		//
		h_2D_KKpt.Fill(kk1pt,kk2pt,hel_weight*bpt_weight);
		//h_2D_KKpt.Fill(kk1pt,kk2pt);
	}

	// Make and return ratio
	h_2D_KKpt_ratio.Divide(&h_2D_KKpt_d,&h_2D_KKpt,histoArea(&h_2D_KKpt),histoArea(&h_2D_KKpt_d));

	stringstream ss_out;
	ss_out << "phi_vs_phi_Iter_" << it_num << "_" << type << ".root";
	TFile f_out(ss_out.str().c_str(),"RECREATE");
	h_2D_KKpt_ratio.Write();
	f_out.Close();

	//cout << "Rndm content" << h_2D_KKpt_ratio.GetBinContent(3)  << endl;
	return h_2D_KKpt_ratio;
}
