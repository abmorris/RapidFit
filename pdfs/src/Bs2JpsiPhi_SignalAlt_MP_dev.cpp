// $Id: Bs2JpsiPhi_SignalAlt_MP_dev.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MP_dev Bs2JpsiPhi_SignalAlt_MP_dev.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-01-28
 */

#include "Bs2JpsiPhi_SignalAlt_MP_dev.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

#define DEBUGFLAG true

#define DOUBLE_TOLERANCE 1E-6

//......................................
//Constructor(s)

//....... 
// old one
Bs2JpsiPhi_SignalAlt_MP_dev::Bs2JpsiPhi_SignalAlt_MP_dev() : 
	  Bs2JpsiPhi_SignalAlt_BaseClass()
	, normalisationCacheValid(false)
{
	MakePrototypes();	
	std::cout << "Constructing PDF: Bs2JpsiPhi_SignalAlt_MP_dev " << std::endl ;
}


//...........
// New with configurator
Bs2JpsiPhi_SignalAlt_MP_dev::Bs2JpsiPhi_SignalAlt_MP_dev( PDFConfigurator config) : 
Bs2JpsiPhi_SignalAlt_BaseClass(config)
, normalisationCacheValid(false)
{
	MakePrototypes();	
	std::cout << "Constructing PDF: Bs2JpsiPhi_SignalAlt_MP_dev " << std::endl ;
}


//.......................................
//Make the data point and parameter set
void Bs2JpsiPhi_SignalAlt_MP_dev::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName.first );
	allObservables.push_back( cosThetaName.first );
	allObservables.push_back( phiName.first );
	allObservables.push_back( cosPsiName.first );
	allObservables.push_back( tagName.first );
	allObservables.push_back( timeAcceptanceCategoryName.first );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName.first );
	parameterNames.push_back( deltaGammaName.first );
	parameterNames.push_back( Aperp_sqName.first );
	parameterNames.push_back( Azero_sqName.first );
	parameterNames.push_back( As_sqName.first );
	parameterNames.push_back( delta_paraName.first );
	parameterNames.push_back( delta_perpName.first );
	parameterNames.push_back( delta_zeroName.first );
	parameterNames.push_back( delta_sName.first );
	parameterNames.push_back( deltaMName.first );

	if( useCosAndSin ) {
		parameterNames.push_back( cosphisName.first );
		parameterNames.push_back( sinphisName.first );
	}
	else{
		parameterNames.push_back( Phi_sName.first );
	}
	
	parameterNames.push_back( mistagName.first );
	parameterNames.push_back( mistagP1Name.first );
	parameterNames.push_back( mistagP0Name.first );
	parameterNames.push_back( mistagSetPointName.first );
	parameterNames.push_back( res1FractionName.first );
	parameterNames.push_back( res1Name.first );
	parameterNames.push_back( res2Name.first );
	parameterNames.push_back( timeOffsetName.first );
	parameterNames.push_back( angAccI1Name.first );
	parameterNames.push_back( angAccI2Name.first );
	parameterNames.push_back( angAccI3Name.first );
	parameterNames.push_back( angAccI4Name.first );
	parameterNames.push_back( angAccI5Name.first );
	parameterNames.push_back( angAccI6Name.first );
	parameterNames.push_back( angAccI7Name.first );
	parameterNames.push_back( angAccI8Name.first );
	parameterNames.push_back( angAccI9Name.first );
	parameterNames.push_back( angAccI10Name.first );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_MP_dev::~Bs2JpsiPhi_SignalAlt_MP_dev()
{
}

//........................................................
//Set the physics parameters into member variables
//Indicate that the cache is no longer valid

bool Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( &gammaName )->GetValue();
    dgam      = allParameters.GetPhysicsParameter( &deltaGammaName )->GetValue();
	
	Azero_sq = allParameters.GetPhysicsParameter( &Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;	}	
	Aperp_sq = allParameters.GetPhysicsParameter( &Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;	}	
	As_sq = allParameters.GetPhysicsParameter( &As_sqName )->GetValue();
	if( (!allowNegativeAsSq&&(As_sq < 0.)) || (As_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: As_sq <0 or >1 but left as is" <<  endl ;	}	

	Apara_sq = (1. - Azero_sq - Aperp_sq  - As_sq) ;
	if( Apara_sq < 0. ) {
		cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}	
		
	delta_zero = allParameters.GetPhysicsParameter( &delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( &delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( &delta_perpName )->GetValue();
	delta_s	   = allParameters.GetPhysicsParameter( &delta_sName )->GetValue();
	delta1 = delta_perp -  delta_para ;    
	delta2 = delta_perp -  delta_zero ;

	_mistag			= allParameters.GetPhysicsParameter( &mistagName )->GetValue();
	_mistagP1		= allParameters.GetPhysicsParameter( &mistagP1Name )->GetValue();
	_mistagP0		= allParameters.GetPhysicsParameter( &mistagP0Name )->GetValue();
	_mistagSetPoint = allParameters.GetPhysicsParameter( &mistagSetPointName )->GetValue();

	delta_ms  = allParameters.GetPhysicsParameter( &deltaMName )->GetValue();

	if(useCosAndSin){
		_cosphis = allParameters.GetPhysicsParameter( &cosphisName )->GetValue();
		_sinphis = allParameters.GetPhysicsParameter( &sinphisName )->GetValue();
	}
	else{
		phi_s     = allParameters.GetPhysicsParameter( &Phi_sName )->GetValue();
		_cosphis = cos(phi_s) ;
		_sinphis = sin(phi_s) ;
	}
	
	// Resolution parameters
	resolution1Fraction = allParameters.GetPhysicsParameter( &res1FractionName )->GetValue();
	resolution1         = allParameters.GetPhysicsParameter( &res1Name )->GetValue();
	resolution2         = allParameters.GetPhysicsParameter( &res2Name )->GetValue();
	timeOffset          = allParameters.GetPhysicsParameter( &timeOffsetName )->GetValue();
	
	// Angular acceptance factors
	angAccI1 = allParameters.GetPhysicsParameter( &angAccI1Name )->GetValue();
	angAccI2 = allParameters.GetPhysicsParameter( &angAccI2Name )->GetValue();
	angAccI3 = allParameters.GetPhysicsParameter( &angAccI3Name )->GetValue();
	angAccI4 = allParameters.GetPhysicsParameter( &angAccI4Name )->GetValue();
	angAccI5 = allParameters.GetPhysicsParameter( &angAccI5Name )->GetValue();
	angAccI6 = allParameters.GetPhysicsParameter( &angAccI6Name )->GetValue();
	angAccI7 = allParameters.GetPhysicsParameter( &angAccI7Name )->GetValue();
	angAccI8 = allParameters.GetPhysicsParameter( &angAccI8Name )->GetValue();
	angAccI9 = allParameters.GetPhysicsParameter( &angAccI9Name )->GetValue();
	angAccI10 = allParameters.GetPhysicsParameter( &angAccI10Name )->GetValue();
	
	// Do a test to ensure user is not using upper time acceptance wrongly
	if( ((fabs(resolution1-0.0)>DOUBLE_TOLERANCE) || (fabs(resolution2-0.0)>DOUBLE_TOLERANCE) || (fabs(mistag()-0.5)>DOUBLE_TOLERANCE) || (fabs(phi_s-0.0)>DOUBLE_TOLERANCE)) && useUpperTimeAcceptance() )
	{
		cout << " You appear to be trying to use the upper time acceptance but are using either resolution or are doing a tagged fit" << endl ;
		cout << " This is not possible at present" << endl ;
		cout << " Resolution1 : " << resolution1 << endl ;
		cout << " Resolution2 : " << resolution2 << endl ;
		cout << " Mistag : " << mistag() << endl ;
		cout << " Phi_s : " << phi_s <<  endl ;
		throw(10);
	}
	
	return result;
}

//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_SignalAlt_MP_dev::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}

//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_SignalAlt_MP_dev::Evaluate(DataPoint * measurement) 
{
	// Get observables into member variables
	t = measurement->GetObservable( &timeName )->GetValue() - timeOffset ;
	ctheta_tr = measurement->GetObservable( &cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( &phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( &cosPsiName )->GetValue();	
	tag = (int)measurement->GetObservable( &tagName )->GetValue();
	timeAcceptanceCategory = (int)measurement->GetObservable( &timeAcceptanceCategoryName )->GetValue();
	
	double val1=0, val2=0 ;
	double returnValue=0 ;
	
	if(resolution1Fraction >= 0.9999 ) {
		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;
		returnValue = this->diffXsec( );
	}
	else {
		// Set the member variable for time resolution to the first value and calculate
		resolution = resolution1 ;
		val1 = this->diffXsec( );
		// Set the member variable for time resolution to the second value and calculate
		resolution = resolution2 ;
		val2 = this->diffXsec( );
		
		returnValue = resolution1Fraction*val1 + (1. - resolution1Fraction)*val2 ;				
	}
	
	//conditions to throw exception
	bool c1 = isnan(returnValue) ;
	bool c2 = ((resolution1>0.)||(resolution2>0.)) && (returnValue <= 0.) ;
	bool c3 = ((fabs(resolution1-0.)<DOUBLE_TOLERANCE)&&((fabs(resolution2-0.)<DOUBLE_TOLERANCE))) && (returnValue <= 0.) && (t>0.) ;
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		cout << endl ;
		cout << " Bs2JpsiPhi_SignalAlt_MP_dev::evaluate() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		cout << " For event with:  " << endl ;
		cout << "   time      " << t << endl ;
		cout << "   ctheta_tr " << ctheta_tr << endl ;
		cout << "   ctheta_1 " << ctheta_1 << endl ;
		cout << "   phi_tr " << phi_tr << endl ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
	if( useLowerTimeAcceptance() ) return returnValue * timeAcceptance.acceptance(t);
	else return returnValue ;
	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_SignalAlt_MP_dev::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary) 
{
		
	// Get observables into member variables
	t = measurement->GetObservable( &timeName )->GetValue() - timeOffset;
	ctheta_tr = measurement->GetObservable( &cosThetaName )->GetValue();
	phi_tr      = measurement->GetObservable( &phiName )->GetValue();
	ctheta_1   = measurement->GetObservable( &cosPsiName )->GetValue();	
	timeAcceptanceCategory = (int)measurement->GetObservable( &timeAcceptanceCategoryName )->GetValue();
	
	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" ) {
		cerr << "Bound on time not provided" << endl;
		exit(1);
	}
	else {
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}
	
	// Recalculate cached values if Physics parameters have changed
	// Must do this for each of the two resolutions.
	if( ! normalisationCacheValid )  {
		for( tag = -1; tag <= 1; ++tag ) {
			if(resolution1Fraction >= 0.9999 ){
				resolution =  resolution1 ;
				normalisationCacheValueRes1[tag+1] = this->diffXsecCompositeNorm1( );
			}
			else {
				resolution =  resolution1 ;
				normalisationCacheValueRes1[tag+1] = this->diffXsecCompositeNorm1( );
				resolution =  resolution2 ;
				normalisationCacheValueRes2[tag+1] = this->diffXsecCompositeNorm1( );
			}
		}
		normalisationCacheValid = true ;
	}	
	
	// calculate return value according to tag 
	tag = (int)measurement->GetObservable( &tagName )->GetValue();
	double returnValue  ;
	if(resolution1Fraction >= 0.9999 )
	{
		returnValue = normalisationCacheValueRes1[tag+1] ;
	}
	else
	{
		returnValue = resolution1Fraction*normalisationCacheValueRes1[tag+1] + (1. - resolution1Fraction)*normalisationCacheValueRes2[tag+1] ;
	}
	
	//conditions to throw exception
	bool c1 = isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;	
	if( DEBUGFLAG && (c1 || c2 ) ) {
		cout << endl ;
		cout << " Bs2JpsiPhi_SignalAlt_MP_dev::Normaisation() returns <=0 or nan :" << returnValue << endl ;
		cout << "   gamma " << gamma() << endl ;
		cout << "   gl    " << gamma_l() << endl ;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl ;
		cout << "   AS^2    " << AS()*AS() << endl ;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
		cout << "   delta_ms       " << delta_ms << endl ;
		cout << "   mistag         " << mistag() << endl ;
		cout << "   mistagP1       " << _mistagP1 << endl ;
		cout << "   mistagP0       " << _mistagP0 << endl ;
		cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
		
	return returnValue ;
}

