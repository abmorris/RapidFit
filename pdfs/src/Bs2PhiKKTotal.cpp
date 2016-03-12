/** @class Bs2PhiKKTotal Bs2PhiKKTotal.cpp
 *
 *  RapidFit PDF for Bs2PhiKKTotal
 *
 *  @author Adam Morris
 *  @date Feb 2016
 */
// Self
#include "Bs2PhiKKTotal.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// ROOT Libraries
#include "TComplex.h"
// RapidFit
#include "PDFConfigurator.h"
#include "DPHelpers.hh"
#define DEBUGFLAG true
PDF_CREATOR( Bs2PhiKKTotal )
/*****************************************************************************/
// Constructor
Bs2PhiKKTotal::Bs2PhiKKTotal(PDFConfigurator* config) :
  // Dependent variables
    mKK     (0.0)
  , phi     (0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Dependent variable names
  , mKKName     ( config->getName("mKK"     ) )
  , phiName     ( config->getName("phi"     ) )
  , ctheta_1Name( config->getName("ctheta_1") )
  , ctheta_2Name( config->getName("ctheta_2") )
  // mKK boundaries
  , mKKmin( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMinimum() )
  , mKKmax( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMaximum() )
  // Options
  , init(true)
  , compIndex(0)
{
  // Set physics parameters to zero for now
  ANonRes   = 0;
  ASsq      = 0;
  deltaS    = 0;
  for(unsigned short i = 0; i < 3; i++)
  {
    APsq[i]   = 0;
    ADsq[i]   = 0;
    deltaP[i] = 0;
    deltaD[i] = 0;
  }
  ANonResName = config->getName("ANonRes2");
  // Magnitude-squared of helicity amplitudes
  ASsqName    = config->getName("ASzero2" );
  APsqName[0] = config->getName("APminus2");
  APsqName[1] = config->getName("APzero2" );
  APsqName[2] = config->getName("APplus2" );
  ADsqName[0] = config->getName("ADminus2");
  ADsqName[1] = config->getName("ADzero2" );
  ADsqName[2] = config->getName("ADplus2" );
  // Phases of helicity amplitudes
  deltaSName    = config->getName("deltaSzero" );
  deltaPName[0] = config->getName("deltaPminus");
  deltaPName[1] = config->getName("deltaPzero" );
  deltaPName[2] = config->getName("deltaPplus" );
  deltaDName[0] = config->getName("deltaDminus");
  deltaDName[1] = config->getName("deltaDzero" );
  deltaDName[2] = config->getName("deltaDplus" );
  MakePrototypes(); // Should only ever go in the constructor. Never put this in the copy constructor!!
  acc = new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile"));
  Initialise();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiKKTotal::Bs2PhiKKTotal(const Bs2PhiKKTotal& copy) : 
    BasePDF( (BasePDF) copy)
  // Dependent variables
  , mKK     (copy.mKK     )
  , phi     (copy.phi     )
  , ctheta_1(copy.ctheta_1)
  , ctheta_2(copy.ctheta_2)
  // Dependent variable names
  , mKKName     (copy.mKKName     )
  , phiName     (copy.phiName     )
  , ctheta_1Name(copy.ctheta_1Name)
  , ctheta_2Name(copy.ctheta_2Name)
  // mKK boundaries
  , mKKmin(copy.mKKmin)
  , mKKmax(copy.mKKmax)
  // Options
  , compIndex(0)
{
  ANonRes     = copy.ANonRes    ;
  ASsq        = copy.ASsq       ;
  deltaS      = copy.deltaS     ;
  ANonResName = copy.ANonResName;
  ASsqName    = copy.ASsqName   ;
  deltaSName  = copy.deltaSName ;
  for(unsigned short i = 0; i < 3; i++)
  {
    APsq[i]       = copy.APsq[i]      ;
    ADsq[i]       = copy.ADsq[i]      ;
    deltaP[i]     = copy.deltaP[i]    ;
    deltaD[i]     = copy.deltaD[i]    ;
    APsqName[i]   = copy.APsqName[i]  ;
    ADsqName[i]   = copy.ADsqName[i]  ;
    deltaPName[i] = copy.deltaPName[i];
    deltaDName[i] = copy.deltaDName[i];
  }
  acc = new LegendreMomentShape(*copy.acc);
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKTotal::~Bs2PhiKKTotal()
{
  delete Swave;
  delete Pwave;
  delete Dwave;
  delete NonRes;
//  delete acc;
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKTotal::Initialise()
{
  // Typical values of barrier factor radius are 3 and 1.7 inverse GeV
  double RBs = 1.7;
  double RKK = 1.7; // TODO: Get these from the config
  double mphi = Bs2PhiKKComponent::mphi;
  // Initialise the signal components
  Swave  = new Bs2PhiKKComponent(0, 980,100    ,"FT",RBs,RKK);
  Pwave  = new Bs2PhiKKComponent(1,mphi,  4.266,"BW",RBs,RKK);
  Dwave  = new Bs2PhiKKComponent(2,1525, 73    ,"BW",RBs,RKK);
  NonRes = new Bs2PhiKKComponent(0,0   ,  0    ,"NR",RBs,RKK);
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( true );
  this->TurnCachingOff();
  componentlist.push_back("Swave");
  componentlist.push_back("Pwave-odd");
  componentlist.push_back("Pwave-even");
  componentlist.push_back("Dwave-odd");
  componentlist.push_back("Dwave-even");
  componentlist.push_back("nonresonant");
  componentlist.push_back("interference");
  SetComponentAmplitudes();
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKTotal::MakePrototypes()
{
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName      );
  allObservables.push_back(phiName      );
  allObservables.push_back(ctheta_1Name );
  allObservables.push_back(ctheta_2Name );
  // Make the parameter set
  vector<string> parameterNames;
  parameterNames.push_back(ANonResName);
  parameterNames.push_back(ASsqName   );
  // Separate loops for P-wave and D-wave for readability in fit output
  for(unsigned short i = 0; i < 3; i++) parameterNames.push_back(APsqName[i]  );
  for(unsigned short i = 0; i < 3; i++) parameterNames.push_back(ADsqName[i]  );
  parameterNames.push_back(deltaSName );
  for(unsigned short i = 0; i < 3; i++) parameterNames.push_back(deltaPName[i]);
  for(unsigned short i = 0; i < 3; i++) parameterNames.push_back(deltaDName[i]);
  allParameters = *( new ParameterSet(parameterNames) );
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKTotal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  // Retrieve the physics parameters
  // Non-resonant
  ANonRes   = allParameters.GetPhysicsParameter(ANonResName)->GetValue();
  // S-wave
  ASsq      = allParameters.GetPhysicsParameter(ASsqName   )->GetValue();
  deltaS    = allParameters.GetPhysicsParameter(deltaSName )->GetValue();
  // P- and D-wave
  for(unsigned short i = 0; i < 3; i++)
  {
    APsq[i]   = allParameters.GetPhysicsParameter(APsqName[i]  )->GetValue();
    ADsq[i]   = allParameters.GetPhysicsParameter(ADsqName[i]  )->GetValue();
    deltaP[i] = allParameters.GetPhysicsParameter(deltaPName[i])->GetValue();
    deltaD[i] = allParameters.GetPhysicsParameter(deltaDName[i])->GetValue();
  }
  // Set the amplitudes
  SetComponentAmplitudes();
  return isOK;
}
/*****************************************************************************/
void Bs2PhiKKTotal::SetComponentAmplitudes()
{
  Swave->SetHelicityAmplitudes(0, sqrt(ASsq), deltaS);
  for(unsigned short i = 0; i < 3; i++)
  {
    Pwave->SetHelicityAmplitudes(i, sqrt(APsq[i]), deltaP[i]);
    Dwave->SetHelicityAmplitudes(i, sqrt(ADsq[i]), deltaD[i]);
  }
  NonRes->SetHelicityAmplitudes(0, sqrt(ANonRes), 0);
//  Swave->Print();
//  Pwave->Print();
//  Dwave->Print();
//  NonRes->Print();
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKTotal::PDFComponents()
{
  return componentlist;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKTotal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
  string compName = component->getComponentName();
  compIndex = component->getComponentNumber();
  if( compIndex == -1)
  {
    component->setComponentNumber(0);
    compIndex = 0;
    for(int i = 0; i < componentlist.size(); i++)
    {
      if(componentlist[i] == compName)
      {
        component->setComponentNumber(i+1);
        compIndex = i+1;
      }
    }
  }
  return Evaluate(measurement);
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKTotal::Evaluate(DataPoint* measurement)
{
  // Get values from the datapoint
  mKK       = measurement->GetObservable(mKKName      )->GetValue();
  phi       = measurement->GetObservable(phiName      )->GetValue();
  ctheta_1  = measurement->GetObservable(ctheta_1Name )->GetValue();
  ctheta_2  = measurement->GetObservable(ctheta_2Name )->GetValue();
  // Check if the datapoint makes sense
  if(phi < -TMath::Pi() || phi > TMath::Pi()
       || ctheta_1 < -1 || ctheta_1 > 1
       || ctheta_2 < -1 || ctheta_2 > 1
       || mKK<2*Bs2PhiKKComponent::mK)
  {
    cout << "Received unphysical datapoint" << endl;
    measurement->Print();
  }
  // Evaluate the PDF at this point
  double evalres = 0;
  double Gamma = 0;
  switch(compIndex)
  {
    case 1:
      // f0(980)
      Gamma =  Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 2:
      // CP-odd phi(1020)
      Gamma =  Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2();
      break;
    case 3:
      // CP-even phi(1020)
      Gamma =  Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2, "even").Rho2();
      break;
    case 4:
      // CP-odd f2'(1525)
      Gamma =  Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2();
      break;
    case 5:
      // CP-even f2'(1525)
      Gamma =  Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2, "even").Rho2();
      break;
    case 6:
      // non-resonant
      Gamma = NonRes->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 7:
      // interference
      Gamma =    TotalAmplitude().Rho2()
            -  Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            - NonRes->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    default:
      Gamma = TotalAmplitude().Rho2();
      break;
  }
  // Phase space part
  const double mK   = Bs2PhiKKComponent::mK;
  const double mBs  = Bs2PhiKKComponent::mBs;
  const double mPhi = Bs2PhiKKComponent::mphi;
  double pR = DPHelpers::daughterMomentum(mKK, mK, mK);
  double pB = DPHelpers::daughterMomentum(mBs,mKK,mPhi);
  double pRpB = pR*pB;
  double phasespace = TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
  evalres = Gamma * Acceptance(mKK, phi, ctheta_1, ctheta_2) * phasespace;
  return evalres>0 && compIndex!=4 ? evalres : 1e-37;
}
/*****************************************************************************/
TComplex Bs2PhiKKTotal::TotalAmplitude()
{
  TComplex amplitude =  Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                     +  Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                     +  Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                     + NonRes->Amplitude(mKK, phi, ctheta_1, ctheta_2);
  return amplitude;
}
/*****************************************************************************/
// Get the angular acceptance
double Bs2PhiKKTotal::Acceptance(double _mKK, double _phi, double _ctheta_1, double _ctheta_2)
{
 return acc->Evaluate(_mKK, _phi, _ctheta_1, _ctheta_2);
}
/*****************************************************************************/
// Insert analitical integral here, if one exists
double Bs2PhiKKTotal::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  double norm = -1;
  return norm;
}
/*****************************************************************************/
