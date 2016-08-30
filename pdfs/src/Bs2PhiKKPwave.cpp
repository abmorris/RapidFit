/** @class Bs2PhiKKPwave Bs2PhiKKPwave.cpp
 *
 *  RapidFit PDF for Bs2PhiKKPwave
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
// Self
#include "Bs2PhiKKPwave.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// ROOT Libraries
#include "TComplex.h"
// RapidFit
#include "PDFConfigurator.h"
#include "DPHelpers.hh"

PDF_CREATOR( Bs2PhiKKPwave )
/*****************************************************************************/
// Constructor
Bs2PhiKKPwave::Bs2PhiKKPwave(PDFConfigurator* config) :
  // Dependent variables
    mKK     (0.0)
  , phi     (0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Dependent variable names
  , mKKName     (config->getName("mKK"     ))
  , phiName     (config->getName("phi"     ))
  , ctheta_1Name(config->getName("ctheta_1"))
  , ctheta_2Name(config->getName("ctheta_2"))
  // Physics parameters
  , Aperpsq(0.0)
  , Azerosq(0.0)
  , Aparasq(0.0)
  , deltaperp(0.0)
  , deltazero(0.0)
  , deltapara(0.0)
  // Magnitude-squared of helicity amplitudes
  , AperpsqName(config->getName("APperp2"))
  , AzerosqName(config->getName("APzero2"))
  // Phases of helicity amplitudes
  , deltaperpName(config->getName("deltaPperp"))
  , deltazeroName(config->getName("deltaPzero"))
  , deltaparaName(config->getName("deltaPpara"))
  // Options
  , compName("0")
{
  MakePrototypes();
  bool drawAll = config->isTrue("DrawAll");
  if(config->isTrue("DrawComponents"))
  {
    componentlist.push_back("0");
    if(drawAll || config->isTrue("DrawPwave"))
    {
      componentlist.push_back("Pwave-even");
      componentlist.push_back("Pwave-odd");
    }
    componentlist.push_back("0");
  }
  Initialise();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiKKPwave::Bs2PhiKKPwave(const Bs2PhiKKPwave& copy) :
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
  // Physics parameters
  , Aperpsq(copy.Aperpsq)
  , Azerosq(copy.Azerosq)
  , Aparasq(copy.Aparasq)
  , deltaperp(copy.deltaperp)
  , deltazero(copy.deltazero)
  , deltapara(copy.deltapara)
  // Magnitude-squared of helicity amplitudes
  , AperpsqName(copy.AperpsqName)
  , AzerosqName(copy.AzerosqName)
  // Phases of helicity amplitudes
  , deltaperpName(copy.deltaperpName)
  , deltazeroName(copy.deltazeroName)
  , deltaparaName(copy.deltaparaName)
  , compName("0")
  , componentlist(copy.componentlist)
{
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKPwave::~Bs2PhiKKPwave()
{
  delete Pwave;
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKPwave::Initialise()
{
  // Typical values of barrier factor radius are 3 and 1.7 inverse GeV
  double RBs = 1.7;
  double RKK = 1.7; // TODO: Get these from the config
  double mphi = Bs2PhiKKComponent::mphi;
  // Initialise the signal components
  Pwave  = new Bs2PhiKKComponent(1,mphi,  4.266,"BW",RBs,RKK);
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( false );
  this->TurnCachingOff();
  SetComponentAmplitudes();
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKPwave::MakePrototypes()
{
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
  // Make the parameter set
  vector<string> parameterNames;
  parameterNames.push_back(AperpsqName);
  parameterNames.push_back(AzerosqName);
  parameterNames.push_back(deltaperpName);
  parameterNames.push_back(deltazeroName);
  parameterNames.push_back(deltaparaName);
  // Separate loops for P-wave and D-wave for readability in fit output
  allParameters = *( new ParameterSet(parameterNames) );
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKPwave::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  // Retrieve the physics parameters
  Aperpsq = allParameters.GetPhysicsParameter(AperpsqName)->GetValue();
  Azerosq = allParameters.GetPhysicsParameter(AzerosqName)->GetValue();
  Aparasq = 1.0 - Aperpsq - Azerosq;
  deltaperp = allParameters.GetPhysicsParameter(deltaperpName)->GetValue();
  deltazero = allParameters.GetPhysicsParameter(deltazeroName)->GetValue();
  deltapara = allParameters.GetPhysicsParameter(deltaparaName)->GetValue();
  SetComponentAmplitudes();
  return isOK;
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKPwave::PDFComponents()
{
  return componentlist;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKPwave::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
  compName = component->getComponentName();
  if( component->getComponentNumber() == -1 )
  {
    // Make the index the same as in the array, for the sake of draw order, or something.
    for(unsigned int i = 0; i < componentlist.size(); i++)
    {
      if(compName == componentlist[i])
      {
        component->setComponentNumber(i);
      }
    }
  }
  return Evaluate(measurement);
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKPwave::Evaluate(DataPoint* measurement)
{
  // Get values from the datapoint
  mKK       = measurement->GetObservable(mKKName     )->GetValue();
  phi       = measurement->GetObservable(phiName     )->GetValue();
  ctheta_1  = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2  = measurement->GetObservable(ctheta_2Name)->GetValue();
  // Check if the datapoint makes sense
  if(phi < -TMath::Pi() || phi > TMath::Pi()
       || ctheta_1 < -1 || ctheta_1 > 1
       || ctheta_2 < -1 || ctheta_2 > 1
  )
  {
    cout << "Received unphysical datapoint" << endl;
    measurement->Print();
  }
  // Evaluate the PDF at this point
  double Gamma = 0;
  if(compName=="Pwave-odd")
  {
    // CP-odd phi(1020)
    Gamma =  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
          +  Pwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2();
  }
  else if(compName=="Pwave-even")
  {
    // CP-even phi(1020)
    Gamma =  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
          +  Pwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2, "even").Rho2();
  }
  else
  {
    Gamma =  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
          +  Pwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2).Rho2();
  }
  return Gamma * PhaseSpace(mKK);
}
/*****************************************************************************/
void Bs2PhiKKPwave::SetComponentAmplitudes()
{
  TComplex Aperp, Azero, Apara;
  Aperp = TComplex(sqrt(Aperpsq),deltaperp,true);
  Azero = TComplex(sqrt(Azerosq),deltazero,true);
  Apara = TComplex(sqrt(Aparasq),deltapara,true);
  TComplex Aplus  = (Apara + Aperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
  TComplex Aminus = (Apara - Aperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
  Pwave->SetHelicityAmplitudes(0, Aminus.Rho(), Aminus.Theta());
  Pwave->SetHelicityAmplitudes(1, Azero.Rho() , Azero .Theta());
  Pwave->SetHelicityAmplitudes(2, Aplus.Rho() , Aplus .Theta());
}
/*****************************************************************************/
double Bs2PhiKKPwave::PhaseSpace(double _mKK)
{
  _mKK/=1000;
  const double mK   = Bs2PhiKKComponent::mK/1000;
  const double mBs  = Bs2PhiKKComponent::mBs/1000;
  const double mPhi = Bs2PhiKKComponent::mphi/1000;
  double pR = DPHelpers::daughterMomentum(_mKK, mK, mK);
  double pB = DPHelpers::daughterMomentum(mBs,_mKK,mPhi);
  double pRpB = pR*pB;
  return TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
}
/*****************************************************************************/
double Bs2PhiKKPwave::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  return 1;
}

