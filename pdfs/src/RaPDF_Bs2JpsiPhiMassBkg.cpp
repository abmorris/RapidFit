// $Id: RaPDF_Bs2JpsiPhiMassBkg.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class RaPDF_Bs2JpsiPhiMassBkg RaPDF_Bs2JpsiPhiMassBkg.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi mass background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "RaPDF_Bs2JpsiPhiMassBkg.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
RaPDF_Bs2JpsiPhiMassBkg::RaPDF_Bs2JpsiPhiMassBkg() : 
	// Physics parameters
	  alphaM_prName	( "alphaM_pr" )
        // Observables
        , recoMassName  ( "mass")
{
	MakePrototypes();
}

//Make the data point and parameter set
void RaPDF_Bs2JpsiPhiMassBkg::MakePrototypes()
{
        allObservables.push_back( recoMassName );

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( alphaM_prName );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
RaPDF_Bs2JpsiPhiMassBkg::~RaPDF_Bs2JpsiPhiMassBkg()
{
}


//Calculate the function value
double RaPDF_Bs2JpsiPhiMassBkg::Evaluate(DataPoint * measurement)
{
  	double alphaM_pr = allParameters.GetPhysicsParameter( alphaM_prName )->GetValue();
	
	// Get the observable
        double mass = measurement->GetObservable( recoMassName )->GetValue();
	
	double val = exp( -alphaM_pr * mass);
  	return val;
}


double RaPDF_Bs2JpsiPhiMassBkg::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	double mhigh, mlow ;
	
	IConstraint * massBound = boundary->GetConstraint("mass");
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on mass not provided in RaPDF_Bs2JpsiPhiMassBkg" << endl;
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}
	
	double alphaM_pr = allParameters.GetPhysicsParameter( alphaM_prName )->GetValue();

	double integral = (1.0/alphaM_pr)* (exp(-alphaM_pr*mlow) - exp(-alphaM_pr*mhigh)) ;
	
	return integral;
	
}
