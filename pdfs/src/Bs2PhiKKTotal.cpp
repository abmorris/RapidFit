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
    mKK(0.0)
  , phi(0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Dependent variable names
  , mKKName       ( config->getName("mKK"     ) )
  , ctheta_1Name  ( config->getName("ctheta_1") )
  , ctheta_2Name  ( config->getName("ctheta_2") )
  , phiName       ( config->getName("phi"     ) )
  // Pre-calculated angular parts
  , ReFSzero(0.0)
  , ReFPminus(0.0)
  , ImFPminus(0.0)
  , ReFPzero(0.0)
  , ReFPplus(0.0)
  , ImFPplus(0.0)
  , ReFDminus(0.0)
  , ImFDminus(0.0)
  , ReFDzero(0.0)
  , ReFDplus(0.0)
  , ImFDplus(0.0)
  // Pre-calculated angular parts names
  , ReFSzeroName ( config->getName("ReFSzero" ) )
  , ReFPminusName( config->getName("ReFPminus") )
  , ImFPminusName( config->getName("ImFPminus") )
  , ReFPzeroName ( config->getName("ReFPzero" ) )
  , ReFPplusName ( config->getName("ReFPplus" ) )
  , ImFPplusName ( config->getName("ImFPplus" ) )
  , ReFDminusName( config->getName("ReFDminus") )
  , ImFDminusName( config->getName("ImFDminus") )
  , ReFDzeroName ( config->getName("ReFDzero" ) )
  , ReFDplusName ( config->getName("ReFDplus" ) )
  , ImFDplusName ( config->getName("ImFDplus" ) )
  // mKK boundaries
  , mKKmin( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMinimum() )
  , mKKmax( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMaximum() )
  // Options
  , init(true)
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
  , mKK(copy.mKK          )
  , phi(copy.phi          )
  , ctheta_1(copy.ctheta_1)
  , ctheta_2(copy.ctheta_2)
  // Dependent variable names
  , mKKName       (copy.mKKName     )
  , ctheta_1Name  (copy.ctheta_1Name)
  , ctheta_2Name  (copy.ctheta_2Name)
  , phiName       (copy.phiName     )
  // Pre-calculated angular parts
  , ReFSzero (copy.ReFSzero )
  , ReFPminus(copy.ReFPminus)
  , ImFPminus(copy.ImFPminus)
  , ReFPzero (copy.ReFPzero )
  , ReFPplus (copy.ReFPplus )
  , ImFPplus (copy.ImFPplus )
  , ReFDminus(copy.ReFDminus)
  , ImFDminus(copy.ImFDminus)
  , ReFDzero (copy.ReFDzero )
  , ReFDplus (copy.ReFDplus )
  , ImFDplus (copy.ImFDplus )
  // Pre-calculated angular part names
  , ReFSzeroName (copy.ReFSzeroName )
  , ReFPminusName(copy.ReFPminusName)
  , ImFPminusName(copy.ImFPminusName)
  , ReFPzeroName (copy.ReFPzeroName )
  , ReFPplusName (copy.ReFPplusName )
  , ImFPplusName (copy.ImFPplusName )
  , ReFDminusName(copy.ReFDminusName)
  , ImFDminusName(copy.ImFDminusName)
  , ReFDzeroName (copy.ReFDzeroName )
  , ReFDplusName (copy.ReFDplusName )
  , ImFDplusName (copy.ImFDplusName )
  // mKK boundaries
  , mKKmin(copy.mKKmin)
  , mKKmax(copy.mKKmax)
{
  ANonRes = copy.ANonRes;
  ASsq = copy.ASsq;
  deltaS = copy.deltaS;
  ANonResName = copy.ANonResName;
  ASsqName = copy.ASsqName;
  deltaSName = copy.deltaSName;
  for(unsigned short i = 0; i < 3; i++)
  {
    APsq[i] = copy.APsq[i];
    ADsq[i] = copy.ADsq[i];
    deltaP[i] = copy.deltaP[i];
    deltaD[i] = copy.deltaD[i];
    APsqName[i] = copy.APsqName[i];
    ADsqName[i] = copy.ADsqName[i];
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
  componentlist.push_back("Pwave");
  componentlist.push_back("Dwave");
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
  allObservables.push_back(ReFSzeroName );
  allObservables.push_back(ReFPminusName);
  allObservables.push_back(ImFPminusName);
  allObservables.push_back(ReFPzeroName );
  allObservables.push_back(ReFPplusName );
  allObservables.push_back(ImFPplusName );
  allObservables.push_back(ReFDminusName);
  allObservables.push_back(ImFDminusName);
  allObservables.push_back(ReFDzeroName );
  allObservables.push_back(ReFDplusName );
  allObservables.push_back(ImFDplusName );
  // Make the parameter set
  vector<string> parameterNames;
  parameterNames.push_back(ANonResName);
  parameterNames.push_back(ASsqName   );
  parameterNames.push_back(deltaSName );
  // Separate loops for P-wave and D-wave for readability in fit output
  for(unsigned short i = 0; i < 3; i++)
  {
    parameterNames.push_back(APsqName[i]  );
    parameterNames.push_back(deltaPName[i]);
  }
  for(unsigned short i = 0; i < 3; i++)
  {
    parameterNames.push_back(ADsqName[i]  );
    parameterNames.push_back(deltaDName[i]);
  }
  allParameters = *( new ParameterSet(parameterNames) );
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKTotal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  double sumq = 0; // Sum in quadriture of amplitudes
  // Retrieve the physics parameters
  // Non-resonant
  ANonRes   = allParameters.GetPhysicsParameter(ANonResName)->GetValue();
  // S-wave
  ASsq      = allParameters.GetPhysicsParameter(ASsqName   )->GetValue();
  deltaS    = allParameters.GetPhysicsParameter(deltaSName )->GetValue();
  sumq     += ANonRes + ASsq;
  // P- and D-wave
  for(unsigned short i = 0; i < 3; i++)
  {
    APsq[i]   = allParameters.GetPhysicsParameter(APsqName[i]  )->GetValue();
    ADsq[i]   = allParameters.GetPhysicsParameter(ADsqName[i]  )->GetValue();
    deltaP[i] = allParameters.GetPhysicsParameter(deltaPName[i])->GetValue();
    deltaD[i] = allParameters.GetPhysicsParameter(deltaDName[i])->GetValue();
    sumq     += APsq[i] + ADsq[i];
  }
  /* Doesn't do what I thought it would do
  // Normalise the amplitudes
  if(abs(sumq-1.0) > 0.00001)
  {
    ANonRes /= sumq;
    ASsq /= sumq;
    allParameters.GetPhysicsParameter(ANonResName)->SetValue(ANonRes);
    allParameters.GetPhysicsParameter(ASsqName   )->SetValue(ASsq);
    for(unsigned short i = 0; i < 3; i++)
    {
      APsq[i] /= sumq;
      ADsq[i] /= sumq;
      allParameters.GetPhysicsParameter(APsqName[i])->SetValue(APsq[i]);
      allParameters.GetPhysicsParameter(ADsqName[i])->SetValue(ADsq[i]);
    }
  }
  */
  // Set the amplitudes
  SetComponentAmplitudes();
  return isOK;
}
/*****************************************************************************/
void Bs2PhiKKTotal::SetComponentAmplitudes()
{
  Swave->SetHelicityAmplitudes(0,sqrt(ASsq), deltaS);
  for(unsigned short i = 0; i < 3; i++)
  {
    Pwave->SetHelicityAmplitudes(i,sqrt(APsq[i]), deltaP[i]);
    Dwave->SetHelicityAmplitudes(i,sqrt(ADsq[i]), deltaD[i]);
  }
  NonRes->SetHelicityAmplitudes(0,sqrt(ANonRes), 0);
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKTotal::PDFComponents()
{
  return componentlist;
}
/*****************************************************************************/
void Bs2PhiKKTotal::ReadDataPoint(DataPoint* measurement)
{
  // Get values from the datapoint
  mKK       = measurement->GetObservable(mKKName      )->GetValue();
  ctheta_1  = measurement->GetObservable(ctheta_1Name )->GetValue();
  ctheta_2  = measurement->GetObservable(ctheta_2Name )->GetValue();
  phi       = measurement->GetObservable(phiName      )->GetValue();
  ReFSzero  = measurement->GetObservable(ReFSzeroName )->GetValue();
  ReFPminus = measurement->GetObservable(ReFPminusName)->GetValue();
  ImFPminus = measurement->GetObservable(ImFPminusName)->GetValue();
  ReFPzero  = measurement->GetObservable(ReFPzeroName )->GetValue();
  ReFPplus  = measurement->GetObservable(ReFPplusName )->GetValue();
  ImFPplus  = measurement->GetObservable(ImFPplusName )->GetValue();
  ReFDminus = measurement->GetObservable(ReFDminusName)->GetValue();
  ImFDminus = measurement->GetObservable(ImFDminusName)->GetValue();
  ReFDzero  = measurement->GetObservable(ReFDzeroName )->GetValue();
  ReFDplus  = measurement->GetObservable(ReFDplusName )->GetValue();
  ImFDplus  = measurement->GetObservable(ImFDplusName )->GetValue();
  // Check if the datapoint makes sense
  if(phi < -TMath::Pi() || phi > TMath::Pi()
       || ctheta_1 < -1 || ctheta_1 > 1
       || ctheta_2 < -1 || ctheta_2 > 1
       || mKK<2*Bs2PhiKKComponent::mK)
  {
    cout << "Received unphysical datapoint:" << endl;
    cout << "m(KK)      :\t" << mKK << endl;
    cout << "phi        :\t" << phi << endl;
    cout << "cos(theta1):\t" << ctheta_1 << endl;
    cout << "cos(theta2):\t" << ctheta_2 << endl;
  }
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKTotal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
  string compName = component->getComponentName();
  compIndex = component->getComponentNumber();
  if( compIndex == -1)
  {
    if(componentlist[0] == compName)
    {
      // f0(980)
      component->setComponentNumber(1);
      compIndex = 1;
    }
    else if(componentlist[1] == compName)
    {
      // phi(1020)
      component->setComponentNumber(2);
      compIndex = 2;
    }
    else if(componentlist[2] == compName)
    {
      // f2'(1525)
      component->setComponentNumber(3);
      compIndex = 3;
    }
    else if(componentlist[3] == compName)
    {
      // non-resonant
      component->setComponentNumber(4);
      compIndex = 4;
    }  
    else if(componentlist[4] == compName)
    {
      // interference
      component->setComponentNumber(5);
      compIndex = 5;
    }
    else
    {
      // total
      component->setComponentNumber(0);
      compIndex = 0;
    }
  }
  return Evaluate(measurement);
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKTotal::Evaluate(DataPoint* measurement)
{
  ReadDataPoint(measurement);
  // Evaluate the PDF at this point
  double evalres = 0;
  double Gamma = 0;
  switch(compIndex)
  {
    case 1:
      // f0(980)
      Gamma = Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 2:
      // phi(1020)
      Gamma = Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 3:
      // f2'(1525)
      Gamma = Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 4:
      // non-resonant
      Gamma = NonRes->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    case 5:
      // interference
      Gamma =    TotalAmplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2()
            - NonRes->Amplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
      break;
    default:
      Gamma = TotalAmplitude(mKK, phi, ctheta_1, ctheta_2).Rho2();
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
TComplex Bs2PhiKKTotal::TotalAmplitude(double _mKK, double _phi, double _ctheta_1, double _ctheta_2)
{
//  TComplex amplitude =  Swave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
//                     +  Pwave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
//                     +  Dwave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
//                     + NonRes->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2);
  TComplex amplitude =  Swave->Amplitude(_mKK, TComplex(0,0), TComplex(ReFSzero,0), TComplex(0,0))
                     +  Pwave->Amplitude(_mKK, TComplex(ReFPminus,ImFPminus), TComplex(ReFPzero,0), TComplex(ReFPplus,ImFPplus))
                     +  Dwave->Amplitude(_mKK, TComplex(ReFDminus,ImFDminus), TComplex(ReFDzero,0), TComplex(ReFDplus,ImFDplus))
                     + NonRes->Amplitude(_mKK, TComplex(0,0), TComplex(ReFSzero,0), TComplex(0,0));
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
