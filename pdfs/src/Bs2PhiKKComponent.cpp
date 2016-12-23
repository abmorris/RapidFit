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
  RBs = std::atof(RBs_str.c_str());
  string RKK_str = config->getConfigurationValue("RKK");
  RKK = std::atof(RKK_str.c_str());
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
  else
  {
    if(lineshape!="NR")
    {
      lineshape = "NR";
      cerr << "Bs2PhiKKComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant." << endl;
    }
  }
  // Helicity amplitude observables
  int lambda_max = TMath::Min(1, _JKK); // Maximum helicity
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
    Ahel.push_back(TComplex(sqrt(1. / (double)n), 0, true));
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
  wignerPhi =  std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
  switch (JKK)
  {
    case 0:
      KKbarrier = std::unique_ptr<DPBarrierL0>(new DPBarrierL0(RKK));
      wignerKK = std::unique_ptr<DPWignerFunctionJ0>(new DPWignerFunctionJ0());
      break;
    case 1:
      KKbarrier = std::unique_ptr<DPBarrierL1>(new DPBarrierL1(RKK));
      wignerKK = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
      break;
    case 2:
      KKbarrier = std::unique_ptr<DPBarrierL2>(new DPBarrierL2(RKK));
      wignerKK = std::unique_ptr<DPWignerFunctionJ2>(new DPWignerFunctionJ2());
      break;
    default:
      cerr << "Bs2PhiKKComponent WARNING: Do not know which barrier factor to use." << endl;
      KKbarrier = std::unique_ptr<DPBarrierL0>(new DPBarrierL0(RKK));
      wignerKK = std::unique_ptr<DPWignerFunctionGeneral>(new DPWignerFunctionGeneral(JKK)); // This shouldn't happen
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
TComplex Bs2PhiKKComponent::A(int lambda)
{
  if(abs(lambda) > helicities.back()) return TComplex(0, 0); //safety
  int i = lambda + helicities.back();
  return Ahel[i];
}
// Angular part of the amplitude
TComplex Bs2PhiKKComponent::F(int lambda, double Phi, double ctheta_1, double ctheta_2)
{
  TComplex exparg = TComplex::I()*Phi;
  exparg*=lambda;
  return wignerPhi->function(ctheta_1, lambda, 0) * wignerKK->function(ctheta_2, lambda, 0) * TComplex::Exp(exparg);
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
  double m0_eff = m_min + (m_max - m_min) * (1 + TMath::TanH((Mres - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
  // Momenta
  double pBs  = DPHelpers::daughterMomentum(mBs,  mphi, mKK   );
  double pKK  = DPHelpers::daughterMomentum(mKK,  mK,   mK    );
  double pBs0 = DPHelpers::daughterMomentum(mBs,  mphi, m0_eff);
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
  if(helicities.size()==0)
    angularPart = TComplex(1,0);
  // Result
  double frac = fraction.value;
  double ofbf = OFBF(mKK);
  /*
  if(std::isnan(frac))
    cerr << this->GetName() << " fraction is not a number." << endl;
  if(std::isnan(massPart.Re()) || std::isnan(massPart.Im()))
    cerr << this->GetName() << " mass-dependent part evaluates to not a number." << endl;
  if(std::isnan(angularPart.Re()) || std::isnan(angularPart.Im()))
    cerr << this->GetName() << " angular part evaluates to not a number." << endl;
  if(std::isnan(ofbf))
    cerr << this->GetName() << " orbital & barrier factors evaluate to not a number." << endl;
  */
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
      Ahel[0] = TComplex(1.0,phases[0].value,true);
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
  for(auto set: {{fraction},magsqs,phases,KKpars})
    for(auto par: set)
      parameters.push_back(par.name);
  return parameters;
}

