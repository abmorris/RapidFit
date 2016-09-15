/** @class Bs2PhiKKSignal Bs2PhiKKSignal.cpp
 *
 *  RapidFit PDF for Bs2PhiKKSignal
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
// Self
#include "Bs2PhiKKSignal.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// ROOT Libraries
#include "TComplex.h"
// RapidFit
#include "PDFConfigurator.h"
#include "DPHelpers.hh"

PDF_CREATOR( Bs2PhiKKSignal )
/*****************************************************************************/
// Constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(PDFConfigurator* config) :
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
  // Variables for physics amplitudes
  // Scaling factors for each component
  , SwaveFrac(0.0)
  , PwaveFrac(0.0)
  , DwaveFrac(0.0)
  // Magnitude-squared of amplitudes
  , ASzerosq(0.0)
  , APperpsq(0.0)
  , APzerosq(0.0)
  , APparasq(0.0)
  , ADperpsq(0.0)
  , ADzerosq(0.0)
  , ADparasq(0.0)
  // Phases of amplitudes
  , deltaSzero(0.0)
  , deltaPperp(0.0)
  , deltaPzero(0.0)
  , deltaPpara(0.0)
  , deltaDperp(0.0)
  , deltaDzero(0.0)
  , deltaDpara(0.0)
  // Fit parameter names for physics amplitudes
  // Scaling factors for each component
  , SwaveFracName(config->getName("SwaveFrac"))
  , PwaveFracName(config->getName("PwaveFrac"))
  , DwaveFracName(config->getName("DwaveFrac"))
  // Magnitude-squared of amplitudes
  , APperpsqName(config->getName("APperp2"))
  , APzerosqName(config->getName("APzero2"))
  , ADperpsqName(config->getName("ADperp2"))
  , ADzerosqName(config->getName("ADzero2"))
  // Phases of amplitudes
  , deltaSzeroName(config->getName("deltaSzero"))
  , deltaPperpName(config->getName("deltaPperp"))
  , deltaPzeroName(config->getName("deltaPzero"))
  , deltaPparaName(config->getName("deltaPpara"))
  , deltaDperpName(config->getName("deltaDperp"))
  , deltaDzeroName(config->getName("deltaDzero"))
  , deltaDparaName(config->getName("deltaDpara"))
  // Options
  , compName("0")
{
  MakePrototypes();
  bool drawAll = config->isTrue("DrawAll");
  if(config->isTrue("DrawComponents"))
  {
    componentlist.push_back("0");
    if(drawAll || config->isTrue("DrawSwave"))
    {
      componentlist.push_back("Swave");
    }
    if(drawAll || config->isTrue("DrawPwave"))
    {
      componentlist.push_back("Pwave-even");
      componentlist.push_back("Pwave-odd");
    }
    if(drawAll || config->isTrue("DrawDwave"))
    {
      componentlist.push_back("Dwave-even");
      componentlist.push_back("Dwave-odd");
    }
    if(drawAll || config->isTrue("DrawInterference"))
    {
      componentlist.push_back("interference");
    }
  }
  Initialise();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(const Bs2PhiKKSignal& copy) :
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
  // Variables for physics amplitudes
  // Scaling factors for each component
  , SwaveFrac(copy.SwaveFrac)
  , PwaveFrac(copy.PwaveFrac)
  , DwaveFrac(copy.DwaveFrac)
  // Magnitude-squared of amplitudes
  , ASzerosq(copy.ASzerosq)
  , APperpsq(copy.APperpsq)
  , APzerosq(copy.APzerosq)
  , APparasq(copy.APparasq)
  , ADperpsq(copy.ADperpsq)
  , ADzerosq(copy.ADzerosq)
  , ADparasq(copy.ADparasq)
  // Phases of amplitudes
//  , deltaSzero(copy.deltaSzero)
  , deltaPperp(copy.deltaPperp)
  , deltaPzero(copy.deltaPzero)
  , deltaPpara(copy.deltaPpara)
  , deltaDperp(copy.deltaDperp)
  , deltaDzero(copy.deltaDzero)
  , deltaDpara(copy.deltaDpara)
  // Fit parameter names for physics amplitudes
  // Scaling factors for each component
  , SwaveFracName(copy.SwaveFracName)
  , PwaveFracName(copy.PwaveFracName)
  , DwaveFracName(copy.DwaveFracName)
  // Magnitude-squared of amplitudes
  , APperpsqName(copy.APperpsqName)
  , APzerosqName(copy.APzerosqName)
  , ADperpsqName(copy.ADperpsqName)
  , ADzerosqName(copy.ADzerosqName)
  // Phases of amplitudes
  , deltaSzeroName(copy.deltaSzeroName)
  , deltaPperpName(copy.deltaPperpName)
  , deltaPzeroName(copy.deltaPzeroName)
  , deltaPparaName(copy.deltaPparaName)
  , deltaDperpName(copy.deltaDperpName)
  , deltaDzeroName(copy.deltaDzeroName)
  , deltaDparaName(copy.deltaDparaName)
  // Drawing components
  , compName("0")
  , componentlist(copy.componentlist)
{
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKSignal::~Bs2PhiKKSignal()
{
  delete Swave;
  delete Pwave;
  delete Dwave;
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKSignal::Initialise()
{
  // Typical values of barrier factor radius are 3 and 1.7 inverse GeV
  double RBs = 1.7;
  double RKK = 1.7; // TODO: Get these from the config
  double mphi = Bs2PhiKKComponent::mphi;
  // Initialise the signal components
  Swave  = new Bs2PhiKKComponent(0, 939.9,0    ,"FT",RBs,RKK);
  Pwave  = new Bs2PhiKKComponent(1,mphi,  4.266,"BW",RBs,RKK);
  Dwave  = new Bs2PhiKKComponent(2,1525, 73    ,"BW",RBs,RKK);
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( true );
  this->TurnCachingOff();
  SetComponentAmplitudes();
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKSignal::MakePrototypes()
{
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
  // Make the parameter set
  vector<string> parameterNames;
  // Scaling factors for each component
  parameterNames.push_back(SwaveFracName);
  parameterNames.push_back(PwaveFracName);
  parameterNames.push_back(DwaveFracName);
  // Magnitude-squared of amplitudes
  parameterNames.push_back(APperpsqName);
  parameterNames.push_back(APzerosqName);
  parameterNames.push_back(ADperpsqName);
  parameterNames.push_back(ADzerosqName);
  // Phases of amplitudes
  parameterNames.push_back(deltaSzeroName);
  parameterNames.push_back(deltaPperpName);
  parameterNames.push_back(deltaPzeroName);
  parameterNames.push_back(deltaPparaName);
  parameterNames.push_back(deltaDperpName);
  parameterNames.push_back(deltaDzeroName);
  parameterNames.push_back(deltaDparaName);
  allParameters = *( new ParameterSet(parameterNames) );
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKSignal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  // Retrieve the physics parameters
  // Scaling factors for each component
  SwaveFrac  = allParameters.GetPhysicsParameter(SwaveFracName)->GetValue();
  PwaveFrac  = allParameters.GetPhysicsParameter(PwaveFracName)->GetValue();
  DwaveFrac  = allParameters.GetPhysicsParameter(DwaveFracName)->GetValue();
  // Magnitude-squared of amplitudes
  ASzerosq   = 1.0;
  APperpsq   = allParameters.GetPhysicsParameter(APperpsqName)->GetValue();
  APzerosq   = allParameters.GetPhysicsParameter(APzerosqName)->GetValue();
  APparasq   = 1.0 - APperpsq - APzerosq;
  ADperpsq   = allParameters.GetPhysicsParameter(ADperpsqName)->GetValue();
  ADzerosq   = allParameters.GetPhysicsParameter(ADzerosqName)->GetValue();
  ADparasq   = 1.0 - ADperpsq - ADzerosq;
  // Phases of amplitudes
  deltaSzero = allParameters.GetPhysicsParameter(deltaSzeroName)->GetValue();
  deltaPperp = allParameters.GetPhysicsParameter(deltaPperpName)->GetValue();
  deltaPzero = allParameters.GetPhysicsParameter(deltaPzeroName)->GetValue();
  deltaPpara = allParameters.GetPhysicsParameter(deltaPparaName)->GetValue();
  deltaDperp = allParameters.GetPhysicsParameter(deltaDperpName)->GetValue();
  deltaDzero = allParameters.GetPhysicsParameter(deltaDzeroName)->GetValue();
  deltaDpara = allParameters.GetPhysicsParameter(deltaDparaName)->GetValue();
  SetComponentAmplitudes();
  return isOK;
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKSignal::PDFComponents()
{
  return componentlist;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKSignal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
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
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
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
//       || mKK < 2*Bs2PhiKKComponent::mK
  )
  {
    cout << "Received unphysical datapoint" << endl;
    measurement->Print();
    Swave->Print();
    Pwave->Print();
    Dwave->Print();
  }
  // Evaluate the PDF at this point
  double Gamma = 0;
  if(compName=="Swave")
  {
    Gamma =  Swave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
          +  Swave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2).Rho2()
          ;
  }
  else if(compName=="Pwave-odd")
  {
    // CP-odd phi(1020)
    Gamma =  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
          +  Pwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
          ;
  }
  else if(compName=="Pwave-even")
  {
    // CP-even phi(1020)
    Gamma =  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
          +  Pwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
          ;
  }
  else if(compName=="Dwave-odd")
  {
    // CP-odd f2'(1525)
    Gamma =  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
          +  Dwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
          ;
  }
  else if(compName=="Dwave-even")
  {
    // CP-even f2'(1525)
    Gamma =  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
          +  Dwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
          ;
  }
  else if(compName=="interference")
  {
    // interference
    Gamma =    TotalAmplitude(false).Rho2()
          -  Swave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
          -  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
          -  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
          // and then the conjugate terms
          +    TotalAmplitude(true).Rho2()
          -  Swave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
          -  Pwave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
          -  Dwave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
          ;
  }
  else
  {
    Gamma = TotalAmplitude(false).Rho2()
          + TotalAmplitude(true ).Rho2()
          ;
  }
  return Gamma * PhaseSpace(mKK)/0.0105327;
}
/*****************************************************************************/
TComplex Bs2PhiKKSignal::TotalAmplitude(bool conjHelAmp)
{
  TComplex amplitude =  Swave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     +  Pwave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     +  Dwave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     ;
  return amplitude;
}
/*****************************************************************************/
void Bs2PhiKKSignal::SetComponentAmplitudes()
{
  TComplex              ASzero(0,0),
           APperp(0,0), APzero(0,0), APpara(0,0),
           ADperp(0,0), ADzero(0,0), ADpara(0,0);
  // S-wave complex amplitude
  ASzero = TComplex(sqrt(ASzerosq),deltaSzero,true)*sqrt(SwaveFrac);
  // P-wave complex amplitudes
  APperp = TComplex(sqrt(APperpsq),deltaPperp,true)*sqrt(PwaveFrac);
  APzero = TComplex(sqrt(APzerosq),deltaPzero,true)*sqrt(PwaveFrac);
  APpara = TComplex(sqrt(APparasq),deltaPpara,true)*sqrt(PwaveFrac);
  TComplex APplus  = (APpara + APperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
  TComplex APminus = (APpara - APperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
  // D-wave complex amplitudes
  ADperp = TComplex(sqrt(ADperpsq),deltaDperp,true)*sqrt(DwaveFrac);
  ADzero = TComplex(sqrt(ADzerosq),deltaDzero,true)*sqrt(DwaveFrac);
  ADpara = TComplex(sqrt(ADparasq),deltaDpara,true)*sqrt(DwaveFrac);
  TComplex ADplus  = (ADpara + ADperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
  TComplex ADminus = (ADpara - ADperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
  // Set the S-wave helicity amplitude
  Swave->SetHelicityAmplitudes(0, ASzero .Rho(), ASzero .Theta());
  // Set the P-wave helicity amplitudes
  Pwave->SetHelicityAmplitudes(0, APminus.Rho(), APminus.Theta());
  Pwave->SetHelicityAmplitudes(1, APzero .Rho(), APzero .Theta());
  Pwave->SetHelicityAmplitudes(2, APplus .Rho(), APplus .Theta());
  // Set the D-wave helicity amplitudes
  Dwave->SetHelicityAmplitudes(0, ADminus.Rho(), ADminus.Theta());
  Dwave->SetHelicityAmplitudes(1, ADzero .Rho(), ADzero .Theta());
  Dwave->SetHelicityAmplitudes(2, ADplus .Rho(), ADplus .Theta());
}
/*****************************************************************************/
double Bs2PhiKKSignal::PhaseSpace(double _mKK)
{
//  _mKK/=1000.;
  const double mK   = Bs2PhiKKComponent::mK;//1000.;
  const double mBs  = Bs2PhiKKComponent::mBs;//1000.;
  const double mPhi = Bs2PhiKKComponent::mphi;//1000.;
  double pR = DPHelpers::daughterMomentum(_mKK, mK, mK);
  double pB = DPHelpers::daughterMomentum(mBs,_mKK,mPhi);
  double pRpB = pR*pB;
  return TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  return -1;
}

