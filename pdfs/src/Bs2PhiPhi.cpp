/** @class Bs2PhiPhi Bs2PhiPhi.cpp
 *
 *  RapidFit PDF for Bs2PhiPhi
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
// Self
#include "Bs2PhiPhi.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// ROOT Libraries
#include "TComplex.h"
// RapidFit
#include "PDFConfigurator.h"

PDF_CREATOR( Bs2PhiPhi )
/*****************************************************************************/
// Constructor
Bs2PhiPhi::Bs2PhiPhi(PDFConfigurator* config) :
  // Dependent variables
    phi     (0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Dependent variable names
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
  , numerical(config->isTrue("Numerical"))
{
  MakePrototypes();
  bool drawAll = config->isTrue("DrawAll");
  if(drawAll || config->isTrue("DrawPwaveComponents"))
  {
    componentlist.push_back("Pwave-even");
    componentlist.push_back("Pwave-odd");
  }
  Initialise();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiPhi::Bs2PhiPhi(const Bs2PhiPhi& copy) :
    BasePDF( (BasePDF) copy)
  // Dependent variables
  , phi     (copy.phi     )
  , ctheta_1(copy.ctheta_1)
  , ctheta_2(copy.ctheta_2)
  // Dependent variable names
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
  , numerical(copy.numerical)
{
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiPhi::~Bs2PhiPhi()
{
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiPhi::Initialise()
{
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( numerical );
  this->TurnCachingOff();
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiPhi::MakePrototypes()
{
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
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
bool Bs2PhiPhi::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  // Retrieve the physics parameters
  Aperpsq = allParameters.GetPhysicsParameter(AperpsqName)->GetValue();
  Azerosq = allParameters.GetPhysicsParameter(AzerosqName)->GetValue();
  Aparasq = 1.0 - Aperpsq - Azerosq;
  deltaperp = allParameters.GetPhysicsParameter(deltaperpName)->GetValue();
  deltazero = allParameters.GetPhysicsParameter(deltazeroName)->GetValue();
  deltapara = allParameters.GetPhysicsParameter(deltaparaName)->GetValue();
  return isOK;
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiPhi::PDFComponents()
{
  return componentlist;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiPhi::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
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
double Bs2PhiPhi::Evaluate(DataPoint* measurement)
{
  // Get values from the datapoint
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
  // Trig identities
  double ctheta_1_sq = ctheta_1*ctheta_1;           // cos²(θ1)
  double stheta_1_sq = 1 - ctheta_1_sq;             // sin²(θ1)
  double s2theta_1   = 2*ctheta_1*sqrt(stheta_1_sq);// sin(2θ1)
  double ctheta_2_sq = ctheta_2*ctheta_2;           // cos²(θ2)
  double stheta_2_sq = 1 - ctheta_2_sq;             // sin²(θ2)
  double s2theta_2   = 2*ctheta_2*sqrt(stheta_2_sq);// sin(2θ2)
  // Coefficients
  TComplex Aperp, Azero, Apara;
  Aperp = TComplex(sqrt(Aperpsq),deltaperp,true);
  Azero = TComplex(sqrt(Azerosq),deltazero,true);
  Apara = TComplex(sqrt(Aparasq),deltapara,true);
  double K[6] = {0};
  K[0] = Azero.Rho2();// |A0|²
  K[1] = Apara.Rho2();// |A‖|²
  K[2] = Aperp.Rho2();// |A⊥|²
  // K[3] is zero
  K[4] = Azero.Rho() * Apara.Rho() * cos((Apara/Azero).Theta());// |A0||A‖|cos(arg(A‖/A0)) or Re(A‖A0*)
  // K[5] is zero
  // Angular functions
  double F[6] = {0};
  F[0] = 4*ctheta_1_sq*ctheta_2_sq;             // 4cos²(θ1)cos²(θ2)
  F[1] = stheta_1_sq*stheta_2_sq*(1+cos(2*phi));// sin²(θ1)sin²(θ2)(1+cos(2Φ))
  F[2] = stheta_1_sq*stheta_2_sq*(1-cos(2*phi));// sin²(θ1)sin²(θ2)(1−cos(2Φ))
  // K[3] is zero so don't calculate F[3]
  F[4] = sqrt(2)*s2theta_1*s2theta_2*cos(phi);  // √2sin(2θ1)sin(2θ2)cos(Φ)
  // K[5] is zero so don't calculate F[5]
  if(compName == "Pwave-odd")
  {
    Gamma = K[2]*F[2];
  }
  else if(compName == "Pwave-even")
  {
    Gamma = K[0]*F[0] + K[1]*F[1] + K[4]*F[4];
  }
  else
  {
    for(int i = 0; i < 6; i++)
    {
      Gamma += K[i] * F[i];
    }
  }
  Gamma *= 9.0 / (32.0*TMath::Pi());
  return Gamma;
}
/*****************************************************************************/
double Bs2PhiPhi::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  return 1;
}

