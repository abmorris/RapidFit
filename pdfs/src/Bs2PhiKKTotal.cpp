/** @class Bs2PhiKKTotal Bs2PhiKKTotal.cpp
 *
 *  RapidFit PDF for Bs2PhiKKTotal
 *
 *  @author Adam Morris
 *  @date Feb 2016
 */
// Self
#include "Bs2PhiKKTotal.h"
#include "Bs2PhiKKNonResonant.h"
// Std Libraries
#include <iostream>
// ROOT Libraries
#include "TComplex.h"
// RapidFit
#include "PDFConfigurator.h"
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
  // mKK boundaries
  , mKKmin( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMinimum() )
  , mKKmax( config->GetPhaseSpaceBoundary()->GetConstraint("mKK")->GetMaximum() )
  // Options
  , init(true)
  , plotComponents( config->isTrue( "PlotComponents" ) )
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
  , mKK(copy.mKK)
  , phi(copy.phi)
  , ctheta_1(copy.ctheta_1)
  , ctheta_2(copy.ctheta_2)
  // Dependent variable names
  , mKKName       (copy.mKKName)
  , ctheta_1Name  (copy.ctheta_1Name)
  , ctheta_2Name  (copy.ctheta_2Name)
  , phiName       (copy.phiName)
  // mKK boundaries
  , mKKmin(copy.mKKmin)
  , mKKmax(copy.mKKmax)
  // Options
  , plotComponents(copy.plotComponents)
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
  Swave = new Bs2PhiKKComponent(0, 980,100    ,"FT",RBs,RKK);
  Pwave = new Bs2PhiKKComponent(1,mphi,  4.266,"BW",RBs,RKK);
  Dwave = new Bs2PhiKKComponent(2,1525, 73    ,"BW",RBs,RKK);
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( true );
  this->TurnCachingOff();
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKTotal::MakePrototypes()
{
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
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
  // Set the amplitudes
  Swave->SetHelicityAmplitudes(0,sqrt(ASsq), deltaS);
  for(unsigned short i = 0; i < 3; i++)
  {
    Pwave->SetHelicityAmplitudes(i,sqrt(APsq[i]), deltaP[i]);
    Dwave->SetHelicityAmplitudes(i,sqrt(ADsq[i]), deltaD[i]);
  }
  return isOK;
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKTotal::PDFComponents()
{
  // The ordering matters for the EvaluateComponent function
  if(plotComponents)
  {
    componentlist.push_back("f0(980)");
    componentlist.push_back("phi(1020)");
    componentlist.push_back("f2'(1525)");
    componentlist.push_back("non-resonant");
    componentlist.push_back("interference");
  }
  return componentlist;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKTotal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
  // Get values from the datapoint
  mKK      = measurement->GetObservable(mKKName     )->GetValue();
  ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
  phi      = measurement->GetObservable(phiName     )->GetValue();
  // Get the real or imaginary part of the amplitude of the component
  double amplitude = 0;
  if(componentlist[0].compare(component->getComponentName()) == 0)
  {
    // f0(980)
    amplitude = Swave->Amplitude(mKK,ctheta_1,ctheta_2,phi).Re();
  }
  if(componentlist[1].compare(component->getComponentName()) == 0)
  {
    // phi(1020)
    amplitude = Pwave->Amplitude(mKK,ctheta_1,ctheta_2,phi).Re();
  }
  if(componentlist[2].compare(component->getComponentName()) == 0)
  {
    // f2'(1525)
    amplitude = Dwave->Amplitude(mKK,ctheta_1,ctheta_2,phi).Re();
  }
  if(componentlist[3].compare(component->getComponentName()) == 0)
  {
    // non-resonant
    amplitude = sqrt(ANonRes*Bs2PhiKKNonResonant::Evaluate(mKK));
  }  
  if(componentlist[4].compare(component->getComponentName()) == 0)
  {
    // interference
    TComplex total = Swave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                   + Pwave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                   + Dwave->Amplitude(mKK, phi, ctheta_1, ctheta_2)
                   + TComplex(sqrt(ANonRes*Bs2PhiKKNonResonant::Evaluate(mKK)),0);
    amplitude = total.Im();
  }
  return amplitude * amplitude * Acceptance(mKK, phi, ctheta_1, ctheta_2) * 1e12;
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKTotal::Evaluate(DataPoint* measurement)
{
  // Get values from the datapoint
  mKK      = measurement->GetObservable(mKKName     )->GetValue();
  ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
  phi      = measurement->GetObservable(phiName     )->GetValue();
  // Check if the datapoint makes sense
  if(phi < -TMath::Pi() || phi > TMath::Pi()
       || ctheta_1 < -1 || ctheta_1 > 1
       || ctheta_2 < -1 || ctheta_2 > 1
       || mKK<2*Bs2PhiKKComponent::mK)
  {
    cout << "Bs2PhiKKComponent::Amplitude() received unphysical datapoint:" << endl;
    cout << "m(KK)      :\t" << mKK << endl;
    cout << "phi        :\t" << phi << endl;
    cout << "cos(theta1):\t" << ctheta_1 << endl;
    cout << "cos(theta2):\t" << ctheta_2 << endl;
  }
  // Evaluate the PDF at this point
  double evalres;
  // Only do convolution around the phi.. or don't do it at all
  evalres = /* TMath::Abs(mKK-Bs2PhiKKComponent::mphi)<20 ? Convolution() : */ EvaluateBase(mKK, phi, ctheta_1, ctheta_2);
  return evalres>0 ? evalres : 1e-9;
}
/*****************************************************************************/
// Base function for evaluation
double Bs2PhiKKTotal::EvaluateBase(double _mKK, double _phi, double _ctheta_1, double _ctheta_2)
{
  // Sum-square of component amplitudes 
  TComplex amplitude = Swave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
                     + Pwave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
                     + Dwave->Amplitude(_mKK, _phi, _ctheta_1, _ctheta_2)
                     + TComplex(sqrt(ANonRes*Bs2PhiKKNonResonant::Evaluate(_mKK)),0);
  double Gamma = amplitude.Rho2()*1e12; // Sensible order of magnitude
  return  (Gamma) * Acceptance(_mKK, _phi, _ctheta_1, _ctheta_2);
}
/*****************************************************************************/
// Do Gaussian convolution
double Bs2PhiKKTotal::Convolution()
{
  double MassResolution = 0.747;
  unsigned int nsteps = 100;
  // Numerical integration method (super slow!)
  // Range of convolution integral
  double xlo = -5 * MassResolution;
  double xhi = +5 * MassResolution;
  // Step size
  double stepSize = (xhi-xlo) / (double)nsteps;
  // Result
  double sum=0.;
  // Shift away from mean so it's not calculated twice
  double running_xlo=xlo-0.5*stepSize;
  double running_xhi=xhi+0.5*stepSize;
  // Do the integral
  for(unsigned int i=0; i < nsteps/2; i++)
  {
    running_xlo += stepSize;
    running_xhi -= stepSize;
    sum += EvaluateBase(mKK - running_xlo, phi, ctheta_1, ctheta_2)
         * TMath::Gaus( running_xlo, 0, MassResolution );
    sum += EvaluateBase(mKK - running_xhi, phi, ctheta_1, ctheta_2)
         * TMath::Gaus( running_xhi, 0, MassResolution );
  }
  return sum * stepSize;
}
/*****************************************************************************/
// Get the angular acceptance
double Bs2PhiKKTotal::Acceptance(double _mKK, double _phi, double _ctheta_1, double _ctheta_2)
{
 return acc->Evaluate(_mKK, _phi, _ctheta_1, _ctheta_2)/0.04;
}
/*****************************************************************************/
// Insert analitical integral here, if one exists
double Bs2PhiKKTotal::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  double norm = -1;
  /*
  norm += ANonRes;
  norm += ASsq;
  for(unsigned short i = 0; i < 3; i++)
  {
    norm+=APsq[i];
    norm+=ADsq[i];
  }
  */
  return norm;
}
/*****************************************************************************/
