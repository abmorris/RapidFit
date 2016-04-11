/** @class Bs2PhiPhi_v4 Bs2PhiPhi_v4.cpp
 *
 *  RapidFit PDF for Bs2PhiPhi_v4
 *
 *  @author Sean Benson
 *  @date 5-10-2011
 */

#include "Bs2PhiPhi_v4.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "PDFConfigurator.h"
#define DEBUGFLAG true
#include "AccCorr_v3.h"
#include "AccCorr_v3.h"

PDF_CREATOR( Bs2PhiPhi_v4 )

//Constructor
Bs2PhiPhi_v4::Bs2PhiPhi_v4(PDFConfigurator* config) : 
     	normalisationCacheValid(false)
    	, evaluationCacheValid(false)
    	,gamma_s(0.0),gamma_l(0.0),gamma_h(0.0),deltaM(0.0),Phi_s(0.0),delta_1(0.0),delta_2(0.0)
    	,Azero_sq(0.0),Aperp_sq(0.0),Apara_sq(0.0),ctheta_1(0.0),ctheta_2(0.0),cach_v1(0.0),cach_v2(0.0)
    	,cach_Azero(0.0),cach_Apara(0.0),cach_Aperp(0.0),cach_sinPhis(0.0),cach_cosPhis(0.0),cach_cosDeltaSubtr(0.0)
    	,phi(0.0),time(0.0),tlo(0.0),thi(0.0),_useAccepWeights(false),weightsVect()
	,intexpcos(0.0),intexpsin(0.0),IntExpGammas(0.0),IntExpGammal(0.0),IntExpGammah(0.0),intcos(0.0),intsin(0.0),w1n(0.0),w2n(0.0),returnValueInit(0.0)
	, deltaGamma(0.0),lambda(0.0)
    	/*
	,tag(0),_mistag(0.0)
	, _mistagP0(0.0), _mistagP1(0.0), _mistagSetPoint(0.0)
    	, _mistagDeltaP0(0.0), _mistagDeltaP1(0.0), _mistagDeltaSetPoint(0.0)
	*/
	,_SS(0.0),_DD(0.0),_CC(0.0)
	,expH(),expL(),expCos(),expSin()
	, resolutionScale(0.0)
	, resolutionOffset(0.0)
	, perEvTime(0.0)
	, resolution(0.0)
	, resolution1(0.0)
	,tlow(0.0),thigh(0.0)
	,Csp(0.0)
	// Physics parameters
	, gamma_sName ( config->getName("gamma") )
	, deltaGammaName ( config->getName("deltaGamma") )
	, Ass_sqName     ( config->getName("Ass_sq") )
	, As_sqName     ( config->getName("As_sq") )
	, deltaSSName     ( config->getName("deltaSS") )
	, deltaSName     ( config->getName("deltaS") )
  	, deltaMName      ( config->getName("deltaM") )
    
	, Phi_sName       ( config->getName("Phi_s") )
    	, Azero_sqName    ( config->getName("Azero_sq") )
    	, Aperp_sqName    ( config->getName("Aperp_sq") )
    	, Apara_sqName    ( config->getName("Apara_sq") )
    	, delta_1Name     ( config->getName("delta_1") )
    	, delta_2Name     ( config->getName("delta_2") )

    	, lambdaName     ( config->getName("lambda") )
    	
	, CspName     ( config->getName("Csp") )

	// Detector parameters
	, resScaleName			( config->getName("timeResolutionScale") )
	, resOffsetName			( config->getName("timeResolutionOffset") )
	, res1Name			( config->getName("timeResolution1") )
	// per event time resolution
	, perEvTimeName			( config->getName("B_s0_LOKI_DTF_CTAUERR") )
    	// Observables
    	, timeName   ( config->getName("time") )
    	, ctheta_1Name( config->getName("ctheta_1") )
    	, ctheta_2Name( config->getName("ctheta_2") )
    	, phiName    ( config->getName("phi") )
    	/*
	, tagName    ( config->getName("tag_calib_Dec13") )
    	, mistagName ( config->getName("mistag_calib_Dec13") )
    	//, tagName    ( config->getName("tag") )
    	//, mistagName ( config->getName("mistag") )
    	, mistagP0Name ( config->getName("mistagP0") )
    	, mistagP1Name ( config->getName("mistagP1") )
    	, mistagSetPointName( config->getName("mistagSetPoint") )
    	, mistagDeltaP0Name ( config->getName("mistagDeltaP0") )
    	, mistagDeltaP1Name ( config->getName("mistagDeltaP1") )
    	, mistagDeltaSetPointName( config->getName("mistagDeltaSetPoint") )
	*/
	// Acceptances
    	,accep1()
    	,timeAcc(NULL)
	// sWave
    	,Ass_sq(0.0),As_sq(0.0),deltaS(0.0),deltaSS(0.0),_useSWAVE(false)
    	,As_sq_in(0.0),Ass_sq_in(0.0)
	// Options
	, _useTimeAcceptance(false)
    	,_useAccepInEval(false)
    	,_usePerEvTime(false)
	// Components
	, _usePlotComponents(false)
	, performingComponentProjection(false)
{
	//this->UseGSLNumericalIntegration( true );
	//this->SetCopyConstructorSafe( false );
	//
	// For the mistag calibration
	mCalib = new MistagCalib3fb(config);
	//
	componentIndex = 0;
	_usePerEvTime = config->isTrue( "perEvTime" ) ;
	_usePlotComponents = config->isTrue( "PlotComponents" ) ;
	//this->SetNumericalNormalisation(true);
	// Do we want to use sWave
	_useSWAVE = config->isTrue( "UseSWave" );
	// Do we want angular acceptance
	_useAccepInEval = config->isTrue( "UseAccepInEval" );
	_useAccepWeights = config->isTrue( "UseAccepWeights" );
	if ( useAccepWeights() ) {
		accep1 = new AccCorr_v3();
		accep1->AccCorr_v3Init(config);
		//accep1->MakeWeights(config);
		weightsVect = accep1->ReturnWeights();
		//if( useAccepInEval() ) {
		//	// Test Petes class
		//	angAccPETE = new AngularAcceptance( config->getConfigurationValue( "AccLocation" ), true ) ;			
		//}
	}
	else {
		for (int k1=0;k1<15;k1++){
			double fillVal=0.0;
			if(k1==0||k1==1||k1==2||k1==6||k1==7) fillVal=1.0;
			weightsVect.push_back(fillVal);
		}
	}
	//*****************************************************************************
	// Configure to use time acceptance machinery 
	_useTimeAcceptance = config->isTrue( "UseTimeAcceptance" ) ;
	
	if( useTimeAcceptance() ) {
		if( config->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0157 ) ;
			cout << "Bs2JpsiPhi_SignalAlt_MO_v4:: Constructing timeAcc: Upper time acceptance beta=0.0157 [0 < t < 14] " << endl ;
		}
		else if( config->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , config->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2PhiPhi_v4:: Constructing timeAcc: using file: " << config->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 10. ) ;
		cout << "Bs2JpsiPhi_SignalAlt_MO_v4:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
	

    MakePrototypes();
    cout << "Making PhiPhi" << endl;
}

//Make the data point and parameter set
void Bs2PhiPhi_v4::MakePrototypes()
{
    	//Make the DataPoint prototype
    	allObservables.push_back( timeName );
    	allObservables.push_back( ctheta_1Name );
    	allObservables.push_back( ctheta_2Name );
    	allObservables.push_back( phiName );
    	//allObservables.push_back( tagName );
    	//allObservables.push_back( mistagName );
	mCalib->addObservables(allObservables);
	if( usePerEvTime() ){
		allObservables.push_back( perEvTimeName );		
	}

    	//Make the parameter set
    	vector<string> parameterNames;
    	parameterNames.push_back( gamma_sName );
    	parameterNames.push_back( deltaGammaName );
    	parameterNames.push_back( Aperp_sqName );
    	parameterNames.push_back( Azero_sqName );
    	if( useSWAVE() ) {
		parameterNames.push_back( Ass_sqName );
		parameterNames.push_back( As_sqName );
		parameterNames.push_back( deltaSSName );
		parameterNames.push_back( deltaSName );
    		parameterNames.push_back( CspName );
	}
    	parameterNames.push_back( delta_1Name );
    	parameterNames.push_back( delta_2Name );
    	parameterNames.push_back( deltaMName );
    	parameterNames.push_back( Phi_sName );
    	parameterNames.push_back( lambdaName );
    	/*
	parameterNames.push_back( mistagP0Name );
    	parameterNames.push_back( mistagP1Name );
	parameterNames.push_back( mistagSetPointName );
    	parameterNames.push_back( mistagDeltaP0Name );
    	parameterNames.push_back( mistagDeltaP1Name );
	parameterNames.push_back( mistagDeltaSetPointName );
	*/
	parameterNames.push_back( resScaleName );
	parameterNames.push_back( resOffsetName );
	parameterNames.push_back( res1Name );

	mCalib->addParameters(parameterNames);

	allParameters = *( new ParameterSet(parameterNames) );
	cout << "Finished making prototypes" << endl;
}

//................................................
//	Copy Constructor
Bs2PhiPhi_v4::Bs2PhiPhi_v4( const Bs2PhiPhi_v4& input )	: BasePDF( (BasePDF) input ),
gamma_sName(input.gamma_sName),
Phi_sName( input.Phi_sName ),Azero_sqName(input.Azero_sqName), Apara_sqName(input.Apara_sqName), Aperp_sqName(input.Aperp_sqName), lambdaName(input.lambdaName),
deltaGammaName(input.deltaGammaName),
//mistagP0Name(input.mistagP0Name),mistagP1Name(input.mistagP1Name),mistagSetPointName(input.mistagSetPointName),
//mistagDeltaP0Name(input.mistagDeltaP0Name),mistagDeltaP1Name(input.mistagDeltaP1Name),mistagDeltaSetPointName(input.mistagDeltaSetPointName),
delta_1Name(input.delta_1Name), delta_2Name(input.delta_2Name), deltaMName(input.deltaMName),
ctheta_1Name(input.ctheta_1Name), ctheta_2Name(input.ctheta_2Name), phiName(input.phiName), 
//mistagName(input.mistagName), tagName(input.tagName),
timeName(input.timeName),
Aperp_sq(input.Aperp_sq), Apara_sq(input.Apara_sq), Azero_sq(input.Azero_sq),gamma_s(input.gamma_s),gamma_l(input.gamma_l),gamma_h(input.gamma_h),deltaM(input.deltaM),Phi_s(input.Phi_s),ctheta_1(input.ctheta_1),ctheta_2(input.ctheta_2),cach_v1(input.cach_v1),cach_v2(input.cach_v2),
delta_1(input.delta_1), delta_2(input.delta_2), 
lambda(input.lambda),
//_mistagP0(input._mistagP0), _mistagP1(input._mistagP1), _mistagSetPoint(input._mistagSetPoint),
//_mistagDeltaP0(input._mistagDeltaP0), _mistagDeltaP1(input._mistagDeltaP1), _mistagDeltaSetPoint(input._mistagDeltaSetPoint),
	deltaGamma(input.deltaGamma)
//,_mistag(input._mistag),tag(input.tag)
,phi(input.phi),cach_Azero(input.cach_Azero),cach_Apara(input.cach_Apara),cach_Aperp(input.cach_Aperp),cach_sinPhis(input.cach_sinPhis),cach_cosPhis(input.cach_cosPhis),cach_cosDeltaSubtr(input.cach_cosDeltaSubtr),
tlo(input.tlo), thi(input.thi),_useTimeAcceptance(input._useTimeAcceptance),_useAccepWeights(input._useAccepWeights),weightsVect(input.weightsVect),w1n(input.w1n),
normalisationCacheValid(input.normalisationCacheValid),evaluationCacheValid(input.evaluationCacheValid),time(input.time),timeAcc(NULL)
// Amplitude Integral members
,w2n(input.w2n),returnValueInit(input.returnValueInit)
,intexpcos(input.intexpcos),intexpsin(input.intexpsin),IntExpGammas(input.IntExpGammas),IntExpGammal(input.IntExpGammal),IntExpGammah(input.IntExpGammah),intcos(input.intcos),intsin(input.intsin)
// sWave
,Ass_sqName(input.Ass_sqName),Ass_sq(input.Ass_sq),As_sqName(input.As_sqName),As_sq(input.As_sq),_useSWAVE(input._useSWAVE)
,deltaSSName(input.deltaSSName),deltaSName(input.deltaSName),deltaSS(input.deltaSS),deltaS(input.deltaS)
,As_sq_in(input.As_sq_in),Ass_sq_in(input.Ass_sq_in)
,_useAccepInEval(input._useAccepInEval)
,accep1(input.accep1)
,angAccPETE(input.angAccPETE)
,resScaleName(input.resScaleName),res1Name(input.res1Name),resOffsetName(input.resOffsetName),perEvTimeName(input.perEvTimeName)
,resolutionScale(input.resolutionScale),resolution(input.resolution),resolution1(input.resolution1)
,resolutionOffset(input.resolutionOffset),perEvTime(input.perEvTime)
,_usePerEvTime(input._usePerEvTime)
,expL(input.expL),expH(input.expH),expSin(input.expSin),expCos(input.expCos)
,tlow(input.tlow),thigh(input.thigh)
,Csp(input.Csp),CspName(input.CspName)
,_SS(input._SS),_DD(input._DD),_CC(input._CC)
, _usePlotComponents(input._usePlotComponents)
, componentIndex(input.componentIndex)
, performingComponentProjection( input.performingComponentProjection )
// For the mistag calibration
//,mCalib(input.mCalib)
{
	timeAcc = new SlicedAcceptance( *(input.timeAcc) );
	mCalib = new MistagCalib3fb(this->GetConfigurator());
}

//Destructor
Bs2PhiPhi_v4::~Bs2PhiPhi_v4()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2PhiPhi_v4::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
    	normalisationCacheValid = false;
    	evaluationCacheValid = false;
    	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
    
    	gamma_s      			= allParameters.GetPhysicsParameter( gamma_sName )->GetValue();
    	deltaGamma      		= allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
    	if( useSWAVE() ) {
		Ass_sq_in      	= allParameters.GetPhysicsParameter( Ass_sqName )->GetValue();
		As_sq_in      	= allParameters.GetPhysicsParameter( As_sqName )->GetValue();
		deltaSS      	= allParameters.GetPhysicsParameter( deltaSSName )->GetValue();
		deltaS      	= allParameters.GetPhysicsParameter( deltaSName )->GetValue();
		Csp      	= allParameters.GetPhysicsParameter( CspName )->GetValue();
	}
    	deltaM       			= allParameters.GetPhysicsParameter( deltaMName )->GetValue();
    	Phi_s        			= allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
    	lambda        			= allParameters.GetPhysicsParameter( lambdaName )->GetValue();
    	Azero_sq     			= allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
    	Aperp_sq     			= allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
    	delta_1 			= allParameters.GetPhysicsParameter( delta_1Name )->GetValue();
    	delta_2 			= allParameters.GetPhysicsParameter( delta_2Name )->GetValue();
    	/*
	_mistagP0 			= allParameters.GetPhysicsParameter( mistagP0Name )->GetValue();
    	_mistagP1 			= allParameters.GetPhysicsParameter( mistagP1Name )->GetValue();
    	_mistagSetPoint			= allParameters.GetPhysicsParameter( mistagSetPointName )->GetValue();
    	_mistagDeltaP0 			= allParameters.GetPhysicsParameter( mistagDeltaP0Name )->GetValue();
    	_mistagDeltaP1 			= allParameters.GetPhysicsParameter( mistagDeltaP1Name )->GetValue();
    	_mistagDeltaSetPoint		= allParameters.GetPhysicsParameter( mistagDeltaSetPointName )->GetValue();
	*/
	mCalib->setParameters(allParameters);
	// Detector parameters
	resolutionScale			= allParameters.GetPhysicsParameter( resScaleName )->GetValue();
	resolutionOffset		= allParameters.GetPhysicsParameter( resOffsetName )->GetValue();
	resolution1         		= allParameters.GetPhysicsParameter( res1Name )->GetValue();
    	
	gamma_l = gamma_s + deltaGamma/2.;
	gamma_h = gamma_s - deltaGamma/2.;

	//Apara_sq = 1.0 - Azero_sq - Aperp_sq - As_sq_in - Ass_sq_in;
	Apara_sq = 1.0 - Azero_sq - Aperp_sq;

    	if( useSWAVE() ) {
    		As_sq=As_sq_in/(1.-As_sq_in-Ass_sq_in);
    		Ass_sq=Ass_sq_in/(1.-As_sq_in-Ass_sq_in);
	}
    	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bs2PhiPhi_v4::GetDoNotIntegrateList()
{
    	vector<string> list;
    	//list.push_back(mistagName);
    	mCalib->addObservables(list);
	if( usePerEvTime() ) list.push_back(perEvTimeName) ;
	return list;
}

//Calculate the function value
double Bs2PhiPhi_v4::Evaluate(DataPoint * measurement)
{
	this->prepareCDS();	
	double returnValue=0.0;
	
    	double evalres1=0.0;
	// Observables (the stuff your experiment measures)
    	ctheta_1 = measurement->GetObservable( ctheta_1Name )->GetValue();
    	ctheta_2 = measurement->GetObservable( ctheta_2Name )->GetValue();
    	phi     = measurement->GetObservable( phiName )->GetValue();
    	//
	//tag = (int)measurement->GetObservable( tagName )->GetValue(); //-1, 0 or +1
    	//_mistag = measurement->GetObservable( mistagName )->GetValue();    
    	mCalib->setObservables(measurement);
	//
	time = measurement->GetObservable( timeName )->GetValue();    
    	if(usePerEvTime()) perEvTime = measurement->GetObservable( perEvTimeName )->GetValue();    
	
	// Min zero resolution
	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = this->EvaluateBase();
	}
	else if(usePerEvTime()){
			resolution = resolutionOffset + perEvTime * resolutionScale;
			evalres1 = this->EvaluateBase();
			returnValue = evalres1 ;				
			//cout << "using resolution (evaluate): " << resolution << "\t";
	}
	else {		
			resolution = resolution1 * resolutionScale ;
			evalres1 = this->EvaluateBase();
			returnValue = evalres1 ;				
	}
	
	// TIME SLICES **************************************************
	if( useTimeAcceptance() ) returnValue = returnValue * timeAcc->getValue(time);
	// ***************************************************************
	/*
	if( std::isnan(evalres) ) {
		cout << "PDF evals to NAN" << endl;
		cout << "littleD = " << littleD << endl;
		cout << "Apara_sq = " << Apara_sq << endl;
		throw 10 ;
	}
	if( evalres <= 0. ) {
		cout << "PDF evals to -ve #" << endl;
		throw 10 ;
	}
	*/
	// USE ACCEPTANCE IN EVAL **************************************************
	// Angular part of PDF
	vector<double> obsvect;
	obsvect.push_back(ctheta_1);
	obsvect.push_back(ctheta_2);
	obsvect.push_back(phi);
	if ( useAccepInEval() ) {
		//return returnValue*accep1->AccEval(obsvect); // Method for my class

		// Method for Pete's class
		Observable* thetaK_obs = measurement->GetObservable( ctheta_1Name );
		Observable* thetaL_obs = measurement->GetObservable( ctheta_2Name );
		Observable* hphi_obs = measurement->GetObservable( phiName );
		//return returnValue*angAccPETE->getValue( thetaK_obs, thetaL_obs, hphi_obs );  // Histogram is generated in PDF basis!
		return 1.05*returnValue*CT1A()*CT2A()*PHIA();
	}
	// ***************************************************************
	//cout << evalres << endl;
	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (time>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	//cout << "D1 = " << D1() << endl;
	//cout << "D2 = " << D2() << endl;
	//mCalib->Print();
	return returnValue;
}

vector<string> Bs2PhiPhi_v4::PDFComponents()
{
	vector<string> this_component_list;
	if( _usePlotComponents ) {
		this_component_list.push_back( "CP-Even" );
		this_component_list.push_back( "CP-Odd" );
		if( allParameters.GetPhysicsParameter(As_sqName)->GetValue() > 1E-10 )
		{
			this_component_list.push_back( "As" );
		}
		this_component_list.push_back( "0" );
	}
	return this_component_list;
}
double Bs2PhiPhi_v4::EvaluateComponent( DataPoint* input, ComponentRef* Component )
{
	performingComponentProjection = true;
	componentIndex = Component->getComponentNumber();
	if( componentIndex == -1 )
	{
		string ComponentName = Component->getComponentName();
		if( ComponentName.compare( "CP-Even" ) == 0 )
		{
			Component->setComponentNumber( 1 );
			componentIndex = 1;
		}
		else if( ComponentName.compare( "CP-Odd" ) == 0 )
		{
			Component->setComponentNumber( 2 );
			componentIndex = 2;
		}
		else if( ComponentName.compare( "As" ) == 0 )
		{
			Component->setComponentNumber( 3 );
			componentIndex = 3;
		}
		else
		{
			Component->setComponentNumber( 0 );
			componentIndex = 0;
		}
	}
	double return_value = this->Evaluate( input );
	componentIndex = 0;
	
	performingComponentProjection = false;
	return return_value;
}

double Bs2PhiPhi_v4::EvaluateBase()  const {
	
	preCalculateTimeFactors();
	double evalres;

	switch( componentIndex )
	{
		case 1:		//	CP-Even		CP-Odd=0 && S-Wave=0
			evalres = HangleFactorA0A0() * timeA0A0()
    	          	+ HangleFactorAPAP() * timeAPAP()
    	          	+ HangleFactorReA0AP() * timeA0AP();
			break;
		case 2:		//	CP-Odd		CP-Even=0 && S-Wave=0
    	          	evalres = HangleFactorATAT() * timeATAT();
			break;
		case 3:		//	S-Wave		CP-Even=0 && CP-Odd=0
			if ( useSWAVE() ) {
				evalres = f7() * timeASSASS()
				+f8() * timeASAS();
			}
			break;
		default:	//	Everything
			//PELC - This turned out to be an important debugging tool
			//switch it on to see the values of PDF being returend.  If ANY go negative, it means there is a sign wrong in one or more of the terms
			//You need to enable in the .h file as well
			//histOfPdfValues->Fill(xsec) ;
			//histCounter++ ;
			//if( histCounter > 10000 ) {
			//	histOfPdfValues->Draw() ;
			//	c0->Update() ;
			//	c0->SaveAs( "histOfPdfValues-from-Evaluate.eps" ) ;
			//	histCounter = 0 ;
			//}
			evalres = HangleFactorA0A0() * timeA0A0()
    			          + HangleFactorAPAP() * timeAPAP()
    			          + HangleFactorATAT() * timeATAT()
    			          + HangleFactorImAPAT() * timeAPAT()
    			          + HangleFactorReA0AP() * timeA0AP()
    			          + HangleFactorImA0AT() * timeA0AT(); 
				if ( useSWAVE() ) {
					evalres += f7() * timeASSASS()
					+f8() * timeASAS()
					+f8a() * timeASAS()
					+f9() * timeASASS()
					+f10() * timeA0ASS()
					+f11() * timeAPASS()
					+f12() * timeATASS()
					+f13() * timeA0AS()
					+f14() * timeAPAS()
					+f15() * timeATAS();
				}
			break;
	}
			/* ORIGINAL
	evalres = HangleFactorA0A0() * timeA0A0()
    	          + HangleFactorAPAP() * timeAPAP()
    	          + HangleFactorATAT() * timeATAT()
    	          + HangleFactorImAPAT() * timeAPAT()
    	          + HangleFactorReA0AP() * timeA0AP()
    	          + HangleFactorImA0AT() * timeA0AT(); 
		//cout << "time factor A0 Apara = " << timeA0AP() << endl;
		//cout << "time factor A0sq = " << timeA0A0() << endl;
		//cout << "time factor Aparasq = " << timeAPAP() << endl;
		//cout << "time factor Aperpsq = " << timeATAT() << endl;
		//cout << "time factor A0 Apara = " << timeA0AP() << endl;
		//cout << "time factor AS ASS = " << timeASASS() << endl;
		//cout << "time factor A0 ASS = " << timeA0ASS() << endl;
		//cout << "Evalres = " << evalres << endl;	
	if ( useSWAVE() ) {
		//cout << "time factor A0 Apara = " << timeA0AP() << endl;
		//cout << "time factor AS ASS = " << timeASASS() << endl;
		//cout << "time factor A0 ASS = " << timeA0ASS() << endl;
		evalres += f7() * timeASSASS()
			+f8() * timeASAS()
			+f8a() * timeASAS()
			+f9() * timeASASS()
			+f10() * timeA0ASS()
			+f11() * timeAPASS()
			+f12() * timeATASS()
			+f13() * timeA0AS()
			+f14() * timeAPAS()
			+f15() * timeATAS();
	}
			*/
	return evalres;
}

//Method for weights generation
double Bs2PhiPhi_v4::EvaluateTimeOnly(DataPoint * measurement)
{
	this->prepareCDS();	
	double returnValue=0.0;
    	double evalres1=0.0;
	// Observables (the stuff your experiment measures)
    	ctheta_1 = measurement->GetObservable( ctheta_1Name )->GetValue();
    	ctheta_2 = measurement->GetObservable( ctheta_2Name )->GetValue();
    	phi     = measurement->GetObservable( phiName )->GetValue();
    	//
	//tag = (int)measurement->GetObservable( tagName )->GetValue(); //-1, 0 or +1
    	//_mistag = measurement->GetObservable( mistagName )->GetValue();    
    	mCalib->setObservables(measurement);
	//
	time = measurement->GetObservable( timeName )->GetValue();    
    	if(usePerEvTime()) perEvTime = measurement->GetObservable( perEvTimeName )->GetValue();    
	// Min zero resolution
	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = this->EvaluateBaseTimeOnly();
	}
	else if(usePerEvTime()){
			resolution = resolutionOffset + perEvTime * resolutionScale;
			evalres1 = this->EvaluateBaseTimeOnly();
			returnValue = evalres1 ;				
	}
	else {		
			resolution = resolution1 * resolutionScale ;
			evalres1 = this->EvaluateBaseTimeOnly();
			returnValue = evalres1 ;				
	}
	
	// TIME SLICES **************************************************
	if( useTimeAcceptance() ) returnValue = returnValue * timeAcc->getValue(time);
	// ***************************************************************
	
	// USE ACCEPTANCE IN EVAL **************************************************
	// Angular part of PDF
	vector<double> obsvect;
	obsvect.push_back(ctheta_1);
	obsvect.push_back(ctheta_2);
	obsvect.push_back(phi);
	//if ( useAccepInEval() ) return returnValue*accep1->AccEval(obsvect);
	// ***************************************************************
	//cout << evalres << endl;
	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (time>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	return returnValue;
}

double Bs2PhiPhi_v4::EvaluateBaseTimeOnly()  const {
	
	preCalculateTimeFactors();
	double evalres;

	evalres = timeA0A0();
    	          + timeAPAP();
    	          + timeATAT();
	
	if ( useSWAVE() ) {
		//cout << "time factor A0 Apara = " << timeA0AP() << endl;
		//cout << "time factor AS ASS = " << timeASASS() << endl;
		//cout << "time factor A0 ASS = " << timeA0ASS() << endl;
		evalres += timeASSASS()
			+ 2.*timeASAS();
	}

	return evalres;
}

void Bs2PhiPhi_v4::AssignTAI(PhaseSpaceBoundary * boundary){
    	setupTimeAmplitudeIntegrals(boundary);
	//cout << Int_timeA0A0() << "\t" << endl;
	returnValueInit =  weightsVect[0]*Int_timeA0A0() + weightsVect[1]*Int_timeAPAP() + weightsVect[2]*Int_timeATAT()
		+ weightsVect[3]*Int_timeAPAT() + weightsVect[4]*Int_timeA0AP() + weightsVect[5]*Int_timeA0AT();

	//cout << "Int_timeA0A0() = " << Int_timeA0A0() << endl;
	//cout << "Int_timeAPAP() = " << Int_timeA0A0() << endl;
	//cout << "Int_timeATAT() = " << Int_timeA0A0() << endl;
	//cout << "IntExpGammal = " << IntExpGammal << endl;
	//cout << "IntExpGammah = " << IntExpGammah << endl;
	//cout << "intexpsin = " << intexpsin << endl;
	//cout << "intexpcos = " << intexpcos << endl;
	if( useSWAVE() ) {
		//cout << "Int. time factor A0 Apara = " << Int_timeA0AP() << endl;
		//cout << "Int. time factor AS ASS = " << Int_timeASASS() << endl;
		//cout << "Int. time factor A0 ASS = " << Int_timeA0ASS() << endl;
		//cout << "deltaSS = " << deltaSS << endl;
		//cout << "deltaS = " << deltaS << endl;
		//cout << "Ass_sq = " << Ass_sq << endl;
		//cout << "As_sq = " << As_sq << endl;
		returnValueInit += weightsVect[6]*Int_timeASSASS()
			+ 2.*(weightsVect[7] + 0.5*2.0*weightsVect[9]*Csp*Csp)*Int_timeASAS() // v2
			//+ 2.*(weightsVect[7] + 2.0*weightsVect[9]*Csp*Csp)*Int_timeASAS() // v1
			//+ 2.*weightsVect[7]*Int_timeASAS() // Original
			+ weightsVect[8]*Int_timeASASS()
			+ weightsVect[9]*Int_timeA0ASS()
			+ weightsVect[10]*Int_timeAPASS()
			+ weightsVect[11]*Int_timeATASS()
			+ weightsVect[12]*Int_timeA0AS() 
			+ weightsVect[13]*Int_timeAPAS() 
			+ weightsVect[14]*Int_timeATAS();
	}
}

double Bs2PhiPhi_v4::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	
	double returnValue = 0.0;
	this->prepareCDS();	
	double normres1=0.0;

	//
	//tag = (int)measurement->GetObservable( tagName )->GetValue(); //-1, 0 or +1
    	//_mistag = measurement->GetObservable( mistagName )->GetValue();    
    	mCalib->setObservables(measurement);
	//
	if(usePerEvTime()) perEvTime = measurement->GetObservable( perEvTimeName )->GetValue();    
  	/*
	if(tag==0){
		_mistag=0.5;
	}
	else{
		_mistag=0.368;
	} 
	*/
	// Min zero resolution
	//return -1.;
	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = NormalisationBase(boundary);
	}
	else if(usePerEvTime()){
			resolution = resolutionOffset + perEvTime * resolutionScale;
			normres1 = this->NormalisationBase(boundary);
			returnValue = normres1 ;				
			//cout << "using resolution (normalisation): " << resolution << "\t";
	}
	else {		
			resolution = resolution1 * resolutionScale ;
			normres1 = NormalisationBase(boundary);
			returnValue = normres1 ;				
	}
	//cout << "Norm. returns " << returnValue << endl;
	// Conditions to throw exception
	bool c1 = std::isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;
	if( DEBUGFLAG && (c1 || c2 ) ) {
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	return returnValue;

}

double Bs2PhiPhi_v4::NormalisationBase(PhaseSpaceBoundary * boundary)
{
    	
	double returnValue=0.;
	
    	// Get time boundaries into member variables
    	IConstraint * timeBound = boundary->GetConstraint( timeName );
    	if ( timeBound->GetUnit() == "NameNotFoundError" ) {
		cerr << "Bound on time not provided" << endl;
		return 0;
    	}
    	else {
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
    	}
    
	double tlo_boundary = tlo ;
    	double thi_boundary = thi ;
	
	// TIME ACCEPTANCE SLICES *******************************************************************************************	
	if( useTimeAcceptance() ) {		    // Set to true because seleting false makes a single slice for 0 --> 14. 
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
		{
			tlo = tlo_boundary > timeAcc->getSlice((int)islice)->tlow() ? tlo_boundary : timeAcc->getSlice((int)islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice((int)islice)->thigh() ? thi_boundary : timeAcc->getSlice((int)islice)->thigh() ;			
			//tlo = timeAcc->getSlice((int)islice)->tlow();
			//thi = timeAcc->getSlice((int)islice)->thigh();			
			if(thi<tlo) continue;
			//if(thi<tlo) cout << "thi = "<< thi << " , tlo = " << tlo << "\n";
			AssignTAI(boundary);
			if( thi > tlo ) returnValue+= returnValueInit * timeAcc->getSlice((int)islice)->height() ;
		}
	}	
	else {
		AssignTAI(boundary);
		returnValue = returnValueInit ;
	}
	// *****************************************************************************************************************
	return returnValue;
}

void Bs2PhiPhi_v4::setupTimeAmplitudeIntegrals( PhaseSpaceBoundary * boundary )
{
   	//double tlow = 0.;
    	//double thigh = 0.;
    	IConstraint * timeBound = boundary->GetConstraint("time");
    	if ( timeBound->GetUnit() == "NameNotFoundError" )
    	{
        	cerr << "Bound on time not provided" << endl;
        	return;
    	}
    	else if ( useTimeAcceptance() )
    	{
        	tlow = tlo;
        	thigh = thi;
    	}
    	else 
    	{
        	tlow = timeBound->GetMinimum();
        	thigh = timeBound->GetMaximum();
    	}

	preCalculateTimeIntegrals();

   	return;
}

void Bs2PhiPhi_v4::preCalculateTimeFactors( ) const
{
	expL = Mathematics::Exp( time, gamma_l, resolution ) ;
	expH = Mathematics::Exp( time, gamma_h, resolution ) ;
	expSin = Mathematics::ExpSin( time, gamma_s, deltaM, resolution ) ;
	expCos = Mathematics::ExpCos( time, gamma_s, deltaM, resolution ) ;
	return ;
}

void Bs2PhiPhi_v4::preCalculateTimeIntegrals( ) const
{
	IntExpGammal = Mathematics::ExpInt( tlow, thigh, gamma_l, resolution )  ;
	IntExpGammah = Mathematics::ExpInt( tlow, thigh, gamma_h, resolution )  ;
	intexpsin = Mathematics::ExpSinInt( tlow, thigh, gamma_s, deltaM, resolution ) ; 
	intexpcos = Mathematics::ExpCosInt( tlow, thigh, gamma_s, deltaM, resolution ) ; 
	return ;
}

//....................................................
// New to prepare all of the coeefficients needed in the time dependen terms
void Bs2PhiPhi_v4::prepareCDS()
{
	double lambda_sq = lambda*lambda;
	double inv_lambda = 1./(1.0 + lambda_sq);

	double F1 = 2.0*lambda *inv_lambda;
	double F2 = (1.0 - lambda_sq) *inv_lambda;

	_SS = sin(Phi_s) * F1;
	_DD = cos(Phi_s) * F1;
	_CC = F2;

}
