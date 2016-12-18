/** @class Bs2PhiKKComponent Bs2PhiKKComponent.cpp
 *
 *  RapidFit PDF for Bs2PhiKKComponent
 *
 *  @author Adam Morris
 *  @date Mar 2016
 */
// Self 
#include "Bs2PhiKKComponent.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
// ROOT Libraries
#include "TMath.h"
// RapidFit Dalitz Plot Libraries
#include "DPBWResonanceShape.hh"
#include "DPFlatteShape.hh"
#include "DPNonresonant.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPHelpers.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"

using std::cout;
using std::cerr;
using std::endl;

double Bs2PhiKKComponent::mBs  = 5.36677;
double Bs2PhiKKComponent::mK   = 0.493677;
double Bs2PhiKKComponent::mpi  = 0.139570;

// Constructor
Bs2PhiKKComponent::Bs2PhiKKComponent(PDFConfigurator* config, string _phiname, string _KKname, int _JKK, string _lineshape) :
    phiname(_phiname)
  , KKname(_KKname)
  , JKK(_JKK)
  , lineshape(_lineshape)
{
  string RBs_str = config->getConfigurationValue("RBs");
  double RBs = std::atof(RBs_str.c_str());
  string RKK_str = config->getConfigurationValue("RKK");
  double RKK = std::atof(RKK_str.c_str());
  fraction = PhysPar(config,KKname+"_fraction");
  // Decide and store the possible helicities
  int lambda_max = TMath::Min(1, _JKK); // Maximum helicity
  for(int lambda = -lambda_max; lambda <= lambda_max; lambda++)
    helicities.push_back(lambda);
  long unsigned int n = helicities.size();
  // Phi resonance parameters
  phipars.push_back(PhysPar(config,phiname+"_mass"));
  phipars.push_back(PhysPar(config,phiname+"_width"));
  // KK resonance parameters
  // Breit Wigner
  if(lineshape=="BW")
  {
    KKLineShape = new DPBWResonanceShape(0, 0, JKK, mK, mK, RKK);
    KKpars.push_back(PhysPar(config,KKname+"_mass"));
    KKpars.push_back(PhysPar(config,KKname+"_width"));
  }
  // Flatte
  else if(lineshape=="FT")
  {
    KKLineShape = new DPFlatteShape(0, 0, mpi, mpi, 0, mK, mK);
    KKpars.push_back(PhysPar(config,KKname+"_mass"));
    KKpars.push_back(PhysPar(config,KKname+"_gpipi"));
    KKpars.push_back(PhysPar(config,KKname+"_Rg"));
  }
  else
  {
    if(lineshape!="NR")
    {
      lineshape = "NR";
      cerr << "Bs2PhiKKComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant." << endl;
    }
    KKLineShape = new DPNonresonant();
  }
  // Helicity amplitude observables
  switch(n)
  {
    case 1:
      magsqs.push_back(PhysPar(config,KKname+"Azerosq"));
      phases.push_back(PhysPar(config,KKname+"deltazero"));
      break;
    case 3:
      magsqs.push_back(PhysPar(config,KKname+"Aperpsq"));
      magsqs.push_back(PhysPar(config,KKname+"Azerosq"));
      phases.push_back(PhysPar(config,KKname+"deltaperp"));
      phases.push_back(PhysPar(config,KKname+"deltazero"));
      phases.push_back(PhysPar(config,KKname+"deltapara"));
      break;
    default:
      throw std::range_error("Bs2PhiKKComponent can't handle this many helicities");
      break;
  }
  // Make the helicity amplitude vector
  while(Ahel.size() < n)
    Ahel.push_back(TComplex(sqrt(1. / (double)n), 0, true));
  Initialise();
}
void Bs2PhiKKComponent::Initialise()
{
  // Build the barrier factor and Wigner function objects
  Bsbarrier = new DPBarrierL0(RBs);
  wignerPhi = new DPWignerFunctionJ1();
  switch (JKK)
  {
    case 0:
      KKbarrier = new DPBarrierL0(RKK);
      wignerKK = new DPWignerFunctionJ0();
      break;
    case 1:
      KKbarrier = new DPBarrierL1(RKK);
      wignerKK = new DPWignerFunctionJ1();
      break;
    case 2:
      KKbarrier = new DPBarrierL2(RKK);
      wignerKK = new DPWignerFunctionJ2();
      break;
    default:
      cerr << "Bs2PhiKKComponent WARNING: Do not know which barrier factor to use." << endl;
      KKbarrier = new DPBarrierL0(RKK);
      wignerKK = new DPWignerFunctionGeneral(JKK); // This shouldn't happen
      break;
  }
}
Bs2PhiKKComponent::Bs2PhiKKComponent(const Bs2PhiKKComponent& other) :
  // Floatable parameters
    fraction(other.fraction)
  , Ahel(other.Ahel)
  , helicities(other.helicities)
  , magsqs(other.magsqs)
  , phases(other.phases)
  , phipars(other.phipars)
  , KKpars(other.KKpars)
  // Fixed parameters
  , phiname(other.phiname)
  , KKname(other.KKname)
  , JKK(other.JKK)
  , lineshape(other.lineshape)
{
  if(other.KKLineShape == nullptr)
  {
    throw std::runtime_error("Somehow the KKLineShape was uninitialised or deleted.");
    return;
  }
  *KKLineShape = *other.KKLineShape;
  Initialise();
}
Bs2PhiKKComponent::~Bs2PhiKKComponent()
{
  delete KKLineShape;
  delete Bsbarrier;
  delete KKbarrier;
  delete wignerPhi;
  delete wignerKK;
}
// Get the corresponding helicity amplitude for a given value of helicity, instead of using array indices
TComplex Bs2PhiKKComponent::A(int lambda)
{
  if(abs(lambda) > helicities.back()) return TComplex(0, 0); //safety
  int i = lambda + helicities.back();
  return Ahel[i];
}
// Angular part of the amplitude
TComplex Bs2PhiKKComponent::F(int lambda, double Phi, double ctheta_1, double ctheta_2)
{
  if(lineshape == "NR")
    return TComplex(1, 0);
  TComplex exparg = TComplex::I()*Phi;
  exparg*=lambda;
  return wignerPhi->function(ctheta_1, lambda, 0) * wignerKK->function(ctheta_2, lambda, 0) * TComplex::Exp(exparg);
}
// Orbital and barrier factor
double Bs2PhiKKComponent::OFBF(double mKK)
{
  if(lineshape=="NR")
    return 1;
  // Orbital factor
  // Masses
  double Mphi = phipars[0].value;
  double Mres = KKpars[0].value;
  double m_min  = mK + mK;
  double m_max  = mBs - Mphi;
  double m0_eff = m_min + (m_max - m_min) * (1 + TMath::TanH((Mres - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
  // Momenta
  double pBs  = DPHelpers::daughterMomentum(mBs,  Mphi, mKK   );
  double pKK  = DPHelpers::daughterMomentum(mKK,  mK,   mK    );
  double pBs0 = DPHelpers::daughterMomentum(mBs,  Mphi, m0_eff);
  double pKK0 = DPHelpers::daughterMomentum(Mres, mK,   mK    );
  double orbitalFactor = //TMath::Power(pBs/mBs,   0)* // == 1 so don't bother
                         TMath::Power(pKK/mKK, JKK);
  // Barrier factors
  double barrierFactor = Bsbarrier->barrier(pBs0, pBs)*
                         KKbarrier->barrier(pKK0, pKK);
  return orbitalFactor * barrierFactor;
}
// The full amplitude
TComplex Bs2PhiKKComponent::Amplitude(double mKK, double phi, double ctheta_1, double ctheta_2)
{
  return Amplitude(mKK, phi, ctheta_1, ctheta_2, "");
}
// The full amplitude with an option
TComplex Bs2PhiKKComponent::Amplitude(double mKK, double phi, double ctheta_1, double ctheta_2, string option)
{
  // Mass-dependent part
  TComplex massPart = KKLineShape->massShape(mKK);
  // Angular part
  TComplex angularPart(0, 0);
  if(option.find("odd") != string::npos || option.find("even") != string::npos)
  {
    TComplex Aperp(sqrt(magsqs[0].value),phases[0].value,true);
    TComplex Apara(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value,true);
    // Temporary helicity amplitudes
    TComplex* HelAmp;
    // CP-odd component
    if(option.find("odd") != string::npos)
      HelAmp = new TComplex[3]{-Aperp/sqrt(2.), TComplex(0, 0), Aperp/sqrt(2.)};
    // CP-even component
    else
      HelAmp = new TComplex[3]{Apara/sqrt(2.), A(0), Apara/sqrt(2.)};
    for(int lambda : helicities)
      angularPart += HelAmp[lambda + helicities.back()] * F(lambda, phi, ctheta_1, ctheta_2);
    delete[] HelAmp;
  }
  else // assume full amplitude, don't thow an error
    for(int lambda : helicities)
      angularPart += A(lambda) * F(lambda, phi, ctheta_1, ctheta_2);
  // Result
  return sqrt(fraction.value) * massPart * angularPart * OFBF(mKK);
}
void Bs2PhiKKComponent::SetPhysicsParameters(ParameterSet& fitpars)
{
  // Update everything from the parameter set
  for(auto set: {{fraction},magsqs,phases,phipars,KKpars})
    for(auto par: set)
      par.Update(fitpars);
  // Update the helicity amplitudes
  switch(helicities.size())
  {
    case 1:
      Ahel[0] = TComplex(sqrt(magsqs[0].value),phases[0].value,true);
      break;
    case 3:
      { // new scope here because we need to declare temporary variables
        TComplex Aperp(sqrt(magsqs[0].value),phases[0].value,true);
        TComplex Azero(sqrt(magsqs[1].value),phases[1].value,true);
        TComplex Apara(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value,true);
        Ahel[0] = (Apara - Aperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
        Ahel[1] = Azero;
        Ahel[2] = (Apara + Aperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
      }
      break;
    default:
      throw std::range_error("Bs2PhiKKComponent can't handle this many helicities");
      break;
  }
  // Update the resonance line shape
  vector<double> respars;
  if(lineshape == "BW")
  {
    respars.push_back(KKpars[0].value); // mass
    respars.push_back(KKpars[1].value); // width
  }
  else if(lineshape == "FT")
  {
    respars.push_back(KKpars[0].value); // mass
    respars.push_back(KKpars[1].value); // gpipi
    respars.push_back(KKpars[1].value*KKpars[2].value); // gKK = gpipi*Rg
  }
  KKLineShape->setParameters(&respars[0]);
}
vector<ObservableRef> Bs2PhiKKComponent::GetPhysicsParameters()
{
  vector<ObservableRef> parameters;
  for(auto set: {{fraction},magsqs,phases,phipars,KKpars})
    for(auto par: set)
      parameters.push_back(par.name);
  return parameters;
}

