// Self
#include "Bs2PhiKKComponent.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
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

double Bs2PhiKKComponent::mBs  = 5.36677;
double Bs2PhiKKComponent::mK   = 0.493677;
double Bs2PhiKKComponent::mpi  = 0.139570;

// Constructor
Bs2PhiKKComponent::Bs2PhiKKComponent(PDFConfigurator* config, std::string _phiname, std::string _KKname, int _JKK, std::string _lineshape) :
    phiname(_phiname)
  , KKname(_KKname)
  , JKK(_JKK)
  , lineshape(_lineshape)
{
  // Barrier factors
  std::string RBs_str = config->getConfigurationValue("RBs");
  std::string RKK_str = config->getConfigurationValue("RKK");
  RBs = std::atof(RBs_str.c_str());
  RKK = std::atof(RKK_str.c_str());
  // Component scale factor
  fraction = PhysPar(config,KKname+"_fraction");
  // KK resonance parameters
  // Breit Wigner
  if(lineshape=="BW")
  {
    KKpars.push_back(PhysPar(config,KKname+"_mass"));
    KKpars.push_back(PhysPar(config,KKname+"_width"));
  }
  // Flatte
  else if(lineshape=="FT")
  {
    KKpars.push_back(PhysPar(config,KKname+"_mass")); 
    KKpars.push_back(PhysPar(config,KKname+"_gpipi"));
    KKpars.push_back(PhysPar(config,KKname+"_Rg"));
  }
  else if(lineshape!="NR")
  {
    lineshape = "NR";
    std::cerr << "Bs2PhiKKComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant." << std::endl;
  }
  // Helicity amplitude observables
  int lambda_max = std::min(1, _JKK); // Maximum helicity
  if(lineshape!="NR")
    for(int lambda = -lambda_max; lambda <= lambda_max; lambda++)
      helicities.push_back(lambda);
  long unsigned int n = helicities.size();
  switch(n)
  {
    case 0:
      break;
    case 1:
      phases.push_back(PhysPar(config,KKname+"_deltazero"));
      break;
    case 3:
      magsqs.push_back(PhysPar(config,KKname+"_Aperpsq"));
      magsqs.push_back(PhysPar(config,KKname+"_Azerosq"));
      phases.push_back(PhysPar(config,KKname+"_deltaperp"));
      phases.push_back(PhysPar(config,KKname+"_deltazero"));
      phases.push_back(PhysPar(config,KKname+"_deltapara"));
      break;
    default:
      throw std::range_error("Bs2PhiKKComponent can't handle this many helicities");
      break;
  }
  // Make the helicity amplitude vector
  while(Ahel.size() < n)
    Ahel.push_back(std::polar<double>(sqrt(1. / (double)n), 0));
  Initialise();
}
void Bs2PhiKKComponent::Initialise()
{
  // Breit Wigner
  if(lineshape=="BW")
    KKLineShape = std::unique_ptr<DPBWResonanceShape>(new DPBWResonanceShape(KKpars[0].value, KKpars[2].value, JKK, mK, mK, RKK));
  // Flatte
  else if(lineshape=="FT")
    KKLineShape = std::unique_ptr<DPFlatteShape>(new DPFlatteShape(KKpars[0].value, KKpars[1].value, mpi, mpi, KKpars[1].value*KKpars[2].value, mK, mK));
  else
    KKLineShape = std::unique_ptr<DPNonresonant>(new DPNonresonant());
  // Build the barrier factor and Wigner function objects
  Bsbarrier = std::unique_ptr<DPBarrierL0>(new DPBarrierL0(RBs));
  wignerPhi = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
  switch (JKK)
  {
    case 0:
      KKbarrier = std::unique_ptr<DPBarrierL0>(new DPBarrierL0(RKK));
      wignerKK  = std::unique_ptr<DPWignerFunctionJ0>(new DPWignerFunctionJ0());
      break;
    case 1:
      KKbarrier = std::unique_ptr<DPBarrierL1>(new DPBarrierL1(RKK));
      wignerKK  = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
      break;
    case 2:
      KKbarrier = std::unique_ptr<DPBarrierL2>(new DPBarrierL2(RKK));
      wignerKK  = std::unique_ptr<DPWignerFunctionJ2>(new DPWignerFunctionJ2());
      break;
    default:
      std::cerr << "Bs2PhiKKComponent WARNING: Do not know which barrier factor to use." << std::endl;
      KKbarrier = std::unique_ptr<DPBarrierL0>(new DPBarrierL0(RKK));
      wignerKK  = std::unique_ptr<DPWignerFunctionGeneral>(new DPWignerFunctionGeneral(JKK)); // This shouldn't happen
      break;
  }
  UpdateParameters();
}
Bs2PhiKKComponent::Bs2PhiKKComponent(const Bs2PhiKKComponent& other) :
  // Floatable parameters
    fraction(other.fraction)
  , Ahel(other.Ahel)
  , helicities(other.helicities)
  , magsqs(other.magsqs)
  , phases(other.phases)
  , mphi(other.mphi)
  , KKpars(other.KKpars)
  // Fixed parameters
  , phiname(other.phiname)
  , KKname(other.KKname)
  , JKK(other.JKK)
  , lineshape(other.lineshape)
{
  Initialise();
}
Bs2PhiKKComponent::~Bs2PhiKKComponent()
{
}
// Get the corresponding helicity amplitude for a given value of helicity, instead of using array indices
std::complex<double> Bs2PhiKKComponent::A(int lambda)
{
  if(abs(lambda) > helicities.back()) return std::complex<double>(0, 0); //safety
  int i = lambda + helicities.back();
  return Ahel[i];
}
// Angular part of the amplitude
std::complex<double> Bs2PhiKKComponent::F(int lambda, double Phi, double ctheta_1, double ctheta_2)
{
  return wignerPhi->function(ctheta_1, lambda, 0) * wignerKK->function(ctheta_2, lambda, 0) * std::polar<double>(1, lambda*Phi);
}
// Orbital and barrier factor
double Bs2PhiKKComponent::OFBF(double mKK)
{
  if(mKK < 2*mK) return 0;
  if(lineshape=="NR")
    return 1;
  // Orbital factor
  // Masses
  double Mres = KKpars[0].value;
  double m_min  = mK + mK;
  double m_max  = mBs - mphi;
  double m0_eff = m_min + (m_max - m_min) * (1 + std::tanh((Mres - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
  // Momenta
  double pBs  = DPHelpers::daughterMomentum(mBs,  mphi, mKK   );
  double pKK  = DPHelpers::daughterMomentum(mKK,  mK,   mK    );
  double pBs0 = DPHelpers::daughterMomentum(mBs,  mphi, m0_eff);
  double pKK0 = DPHelpers::daughterMomentum(Mres, mK,   mK    );
  double orbitalFactor = //std::pow(pBs/mBs,   0)* // == 1 so don't bother
                         std::pow(pKK/mKK, JKK);
  // Barrier factors
  double barrierFactor = Bsbarrier->barrier(pBs0, pBs)*
                         KKbarrier->barrier(pKK0, pKK);
  return orbitalFactor * barrierFactor;
}
// The full amplitude
std::complex<double> Bs2PhiKKComponent::Amplitude(double mKK, double phi, double ctheta_1, double ctheta_2)
{
  return Amplitude(mKK, phi, ctheta_1, ctheta_2, "");
}
// The full amplitude with an option
std::complex<double> Bs2PhiKKComponent::Amplitude(double mKK, double phi, double ctheta_1, double ctheta_2, std::string option)
{
  // Mass-dependent part
  std::complex<double> massPart = KKLineShape->massShape(mKK);
  // Angular part
  std::complex<double> angularPart(0, 0);
  if(option.find("odd") != std::string::npos || option.find("even") != std::string::npos)
  {
    std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
    std::complex<double> Apara = std::polar(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value);
    // Temporary helicity amplitudes
    std::vector<std::complex<double>> HelAmp;
    // CP-odd component
    if(option.find("odd") != std::string::npos)
      HelAmp = {-Aperp/sqrt(2.), std::complex<double>(0, 0), Aperp/sqrt(2.)};
    // CP-even component
    else
      HelAmp = {Apara/sqrt(2.), A(0), Apara/sqrt(2.)};
    for(int lambda : helicities)
      angularPart += HelAmp[lambda + helicities.back()] * F(lambda, phi, ctheta_1, ctheta_2);
  }
  else // assume full amplitude, don't thow an error
    for(int lambda : helicities)
      angularPart += A(lambda) * F(lambda, phi, ctheta_1, ctheta_2);
  if(helicities.size()==0)
    angularPart = std::complex<double>(1,0);
  // Result
  double frac = fraction.value;
  double ofbf = OFBF(mKK);
  return frac * massPart * angularPart * ofbf;
}
void Bs2PhiKKComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
  // Update everything from the parameter set
  mphi = fitpars->GetPhysicsParameter(phiname+"_mass")->GetValue();
  fraction.Update(fitpars);
  for(auto& par: magsqs) par.Update(fitpars);
  for(auto& par: phases) par.Update(fitpars);
  for(auto& par: KKpars) par.Update(fitpars);
  UpdateParameters();
}
void Bs2PhiKKComponent::UpdateParameters()
{
  // Update the helicity amplitudes
  switch(helicities.size())
  {
    case 0:
      break;
    case 1:
      Ahel[0] = std::polar<double>(1.0,phases[0].value);
      break;
    case 3:
      { // new scope here because we need to declare temporary variables
        std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
        std::complex<double> Azero = std::polar(sqrt(magsqs[1].value),phases[1].value);
        std::complex<double> Apara = std::polar(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value);
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
  for(auto set: {{fraction},magsqs,phases,KKpars})
    for(auto par: set)
      parameters.push_back(par.name);
  return parameters;
}

