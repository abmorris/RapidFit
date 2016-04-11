/** @class Bs2PhiKKTotal Bs2PhiKKTotal.cpp
 *
 *  RapidFit PDF for Bs2PhiKKTotal
 *
 *  @author Adam Morris
 *  @date Mar 2016
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
  , compName("0")
  , useTimeIntPwavePDF(config->isTrue("UseTimeIntegratedPwavePDF"))
  , useTimeIntDwavePDF(config->isTrue("UseTimeIntegratedDwavePDF"))
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
  if(useTimeIntPwavePDF && useTimeIntDwavePDF)
  {
    cout << "WARNING: Cannot currently plot P-wave and D-wave. Will only plot P-wave PDF." << endl;
  }
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
  , useTimeIntPwavePDF(copy.useTimeIntPwavePDF)
  , useTimeIntDwavePDF(copy.useTimeIntDwavePDF)
  , compName("0")
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
  delete acc;
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
  componentlist.push_back("Pwave-even");
  componentlist.push_back("Pwave-odd");
  componentlist.push_back("Dwave-even");
  componentlist.push_back("Dwave-odd");
  componentlist.push_back("nonresonant");
  componentlist.push_back("interference");
//  componentlist.push_back("acceptance"); // Don't expect to use this very much 
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
  // Phases
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
  compName = component->getComponentName();
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
  double Gamma = 0;
  if(useTimeIntPwavePDF || useTimeIntDwavePDF)
  {
    // Trig identities
    double ctheta_1_sq = ctheta_1*ctheta_1;          // cos²(θ1)
    double stheta_1_sq = 1 - ctheta_1_sq;            // sin²(θ1)
    double s2theta_1 = 2*ctheta_1*sqrt(stheta_1_sq); // sin(2θ1)
    double ctheta_2_sq = ctheta_2*ctheta_2;          // cos²(θ2)
    double stheta_2_sq = 1 - ctheta_2_sq;            // sin²(θ2)
    double s2theta_2 = 2*ctheta_2*sqrt(stheta_2_sq); // sin(2θ2)
    TComplex Aminus, Azero, Aplus;
    if(useTimeIntPwavePDF)
    {
      Aminus = TComplex(sqrt(APsq[0]),deltaP[0],true);
      Azero  = TComplex(sqrt(APsq[1]),deltaP[1],true);
      Aplus  = TComplex(sqrt(APsq[2]),deltaP[2],true);
    }
    else
    {
      Aminus = TComplex(sqrt(ADsq[0]),deltaD[0],true);
      Azero  = TComplex(sqrt(ADsq[1]),deltaD[1],true);
      Aplus  = TComplex(sqrt(ADsq[2]),deltaD[2],true);
    }
    TComplex Apara = (Aplus + Aminus)/sqrt(2.); // A‖ = (A+ + A−)/sqrt(2)
    TComplex Aperp = (Aplus - Aminus)/sqrt(2.); // A⊥ = (A+ − A−)/sqrt(2)
    // Coefficients 
    double K[6] = {0};
    K[0] = Azero.Rho2(); // |A0|²
    K[1] = Apara.Rho2(); // |A‖|²
    K[2] = Aperp.Rho2(); // |A⊥|²
    // K[3] is zero
    K[4] = Azero.Rho() * Apara.Rho() * cos((Apara/Azero).Theta()); // |A0||A‖|cos(arg(A‖/A0)) or Re(A‖A0*)
    // K[5] is zero
    // Angular functions
    double F[6] = {0};
    if(useTimeIntPwavePDF)
    {
      F[0] = 4*ctheta_1_sq*ctheta_2_sq;                // 4cos²(θ1)cos²(θ2)
      F[1] = stheta_1_sq*stheta_2_sq*(1+cos(2*phi));   // sin²(θ1)sin²(θ2)(1+cos(2Φ))
      F[2] = stheta_1_sq*stheta_2_sq*(1-cos(2*phi));   // sin²(θ1)sin²(θ2)(1−cos(2Φ))
      // K[3] is zero so don't calculate F[3]
      F[4] = sqrt(2)*s2theta_1*s2theta_2*cos(phi);     // √2sin(2θ1)sin(2θ2)cos(Φ)
      // K[5] is zero so don't calculate F[5]
    }
    else
    {
      F[0] = 4*ctheta_1_sq*pow(3*ctheta_2_sq-1,2);                    // 4cos²(θ1)(3cos²(θ2)-1)²
      F[1] = 3*stheta_1_sq*stheta_2_sq*(1+cos(2*phi));                // 3sin²(θ1)sin²(θ2)(1+cos(2Φ))
      F[2] = 3*stheta_1_sq*stheta_2_sq*(1-cos(2*phi));                // 3sin²(θ1)sin²(θ2)(1−cos(2Φ))
      // K[3] is zero so don't calculate F[3]
      F[4] = 2*sqrt(6)*s2theta_1*s2theta_2*(3*ctheta_2_sq-1)*cos(phi); // 2√6sin(2θ1)sin(2θ2)(3*ctheta_2_sq-1)cos(Φ)
      // K[5] is zero so don't calculate F[5]
    }
    if(compName == "Pwave-odd" || compName == "Dwave-odd")
    {
      Gamma += K[2]*F[2];
    }
    else if (compName == "Pwave-even" || compName == "Dwave-even")
    {
      Gamma += K[0]*F[0] + K[1]*F[1] + K[4]*F[4];
    }
    else
    {
      for(int i = 0; i < 6; i++)
      {
        Gamma += K[i] * F[i];
      }
    }
  }
  else
  {
    // Get the square of the amplitude for the chosen component
    if(compName=="Swave")
    {
      // f0(980)
      Gamma =  Swave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            +  Swave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2).Rho2();
    }
    else if(compName=="Pwave-odd")
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
    else if(compName=="Dwave-odd")
    {
      // CP-odd f2'(1525)
      Gamma =  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2()
            +  Dwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2,  "odd").Rho2();
    }
    else if(compName=="Dwave-even")
    {
      // CP-even f2'(1525)
      Gamma =  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2, "even").Rho2()
            +  Dwave->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2, "even").Rho2();
    }
    else if(compName=="nonresonant")
    {
      // non-resonant
      Gamma = NonRes->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            + NonRes->Amplitude(true,  mKK, phi, ctheta_1, ctheta_2).Rho2();
    }
    else if(compName=="interference")
    {
      // interference
      Gamma =    TotalAmplitude(false).Rho2()
            -  Swave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Pwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Dwave->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            - NonRes->Amplitude(false, mKK, phi, ctheta_1, ctheta_2).Rho2()
            // and then the conjugate terms
            +    TotalAmplitude(true).Rho2()
            -  Swave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Pwave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
            -  Dwave->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2()
            - NonRes->Amplitude(true, mKK, phi, ctheta_1, ctheta_2).Rho2();
    }
    else if(compName=="acceptance")
      Gamma = 2*Pwave->Amplitude(false, 1030, 0, 0, 0).Rho2()*PhaseSpace(1015)/PhaseSpace(mKK); // Just a visual check
    else
    {
      Gamma = TotalAmplitude(true).Rho2() + TotalAmplitude(false).Rho2();
    }
  }
  Gamma/=2.0;
//  return Gamma * PhaseSpace(mKK);
  return Gamma * Acceptance(mKK, phi, ctheta_1, ctheta_2) * PhaseSpace(mKK);
}
/*****************************************************************************/
TComplex Bs2PhiKKTotal::TotalAmplitude(bool conjHelAmp)
{
  TComplex amplitude =  Swave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     +  Pwave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     +  Dwave->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2)
                     + NonRes->Amplitude(conjHelAmp, mKK, phi, ctheta_1, ctheta_2);
  return amplitude;
}
/*****************************************************************************/
// Get the angular acceptance
double Bs2PhiKKTotal::Acceptance(double _mKK, double _phi, double _ctheta_1, double _ctheta_2)
{
 return acc->Evaluate(_mKK, _phi, _ctheta_1, _ctheta_2);
}
/*****************************************************************************/
// Two-body phase space probability function
double Bs2PhiKKTotal::PhaseSpace(double _mKK)
{
  const double mK   = Bs2PhiKKComponent::mK;
  const double mBs  = Bs2PhiKKComponent::mBs;
  const double mPhi = Bs2PhiKKComponent::mphi;
  double pR = DPHelpers::daughterMomentum(_mKK, mK, mK);
  double pB = DPHelpers::daughterMomentum(mBs,_mKK,mPhi);
  double pRpB = pR*pB;
  return TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
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
