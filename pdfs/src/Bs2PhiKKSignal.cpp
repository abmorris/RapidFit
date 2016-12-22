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
#include <stdexcept>
// ROOT Libraries
#include "TComplex.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
// RapidFit
#include "StringProcessing.h"
#include "PDFConfigurator.h"
#include "DPHelpers.hh"

PDF_CREATOR( Bs2PhiKKSignal )
/*****************************************************************************/
// Constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(PDFConfigurator* config) :
  // Dependent variables
    mKK     (0.)
  , phi     (0.)
  , ctheta_1(0.)
  , ctheta_2(0.)
  // Dependent variable names
  , mKKName     (config->getName("mKK"     ))
  , phiName     (config->getName("phi"     ))
  , ctheta_1Name(config->getName("ctheta_1"))
  , ctheta_2Name(config->getName("ctheta_2"))
  , acceptance_moments((string)config->getConfigurationValue("CoefficientsFile") != "")
  , acceptance_histogram((string)config->getConfigurationValue("HistogramFile") != "")
{
  string phiname = config->getConfigurationValue("phiname");
  phimass = PhysPar(config,phiname+"_mass");
  dGsGs = PhysPar(config,"dGsGs");
  vector<string> reslist = StringProcessing::SplitString( config->getConfigurationValue("resonances"), ' ' );
  printf("========================\n");
  printf("Building resonance model\n");
  printf("------------------------\n");
  for(auto name: reslist)
  {
    if(name=="") continue;
    Bs2PhiKKComponent comp = ParseComponent(config,phiname,name);
    components.push_back(comp);
    componentnames.push_back(comp.GetName());
  }
  printf("========================\n");
  if(componentnames.size() > 1)
    componentnames.push_back("interference");
  if(acceptance_moments) acc_m = unique_ptr<LegendreMomentShape>(new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile")));
  else if(acceptance_histogram)
  {
    TFile* histfile = TFile::Open(config->getConfigurationValue("HistogramFile").c_str());
    acc_h = shared_ptr<NDHist_Adaptive>(new NDHist_Adaptive(histfile));
    acc_h->LoadFromTree((TTree*)histfile->Get("AccTree"));
    acc_h->SetDimScales({1e-3,0.1,1.,1.}); // TODO: read this from the file
  }
  printf("Calling Initialiser\n");
  Initialise();
  printf("Making Prototypes\n");
  MakePrototypes();
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
  // Width splitting
  , dGsGs(copy.dGsGs)
  // Phi mass
  , phimass(copy.phimass)
  // PDF components
  , components(copy.components)
  // Plotting components
  , componentnames(copy.componentnames)
  // Options
  , acceptance_moments(copy.acceptance_moments)
  , acceptance_histogram(copy.acceptance_histogram)
{
  if(acceptance_moments) acc_m = unique_ptr<LegendreMomentShape>(new LegendreMomentShape(*copy.acc_m));
  else if(acceptance_histogram) acc_h = copy.acc_h;
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKSignal::~Bs2PhiKKSignal()
{
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKSignal::Initialise()
{
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( true );
  this->TurnCachingOff();
}
/*****************************************************************************/
// Build a component object from a passed option
Bs2PhiKKComponent Bs2PhiKKSignal::ParseComponent(PDFConfigurator* config, string phiname, string option)
{
  // Syntax: <resonance name>(<spin>,<lineshape>)
  // - the list is space-delimited: no extra spaces, please!
  // - resonance name must be alphanumeric
  // - spin must be a single digit 0, 1 or 2, or you have to implement the code for higher spins
  // - lineshape name must be 2 capital letters
  // - see Bs2PhiKKComponent.cpp for implemented lineshapes (BW, FT, NR... but more can be added)
  // Example: "phi1020(1,BW)"
  string KKname = option.substr(0,option.find('('));
  int JKK = atoi(option.substr(option.find('(')+1,1).c_str());
  string lineshape = option.substr(option.find(',')+1,2);
  printf("Bs -> %s %s \t (spin-%i %s shape)\n", phiname.c_str(), KKname.c_str(), JKK, lineshape.c_str());
  return Bs2PhiKKComponent(config, phiname, KKname, JKK, lineshape);
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKSignal::MakePrototypes()
{
  // Make the DataPoint prototype
  printf("Prototype datapoint:\t%s, %s, %s, %s\n",mKKName.Name().c_str(),phiName.Name().c_str(),ctheta_1Name.Name().c_str(),ctheta_2Name.Name().c_str());
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
  // Make the parameter set
  vector<string> parameterNames;
  // Resonance parameters
  parameterNames.push_back(dGsGs.name);
  parameterNames.push_back(phimass.name);
  for(auto comp: components)
    for(string par: comp.GetPhysicsParameters())
      if(par!=phimass.name.Name()) parameterNames.push_back(par);
  printf("Full parameter set\n");
  for(string par: parameterNames)
    printf("%s\n",par.c_str());
  allParameters = *( new ParameterSet(parameterNames) );
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKSignal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  dGsGs.Update(&allParameters);
  phimass.Update(&allParameters);
  for(auto& comp: components)
    comp.SetPhysicsParameters(&allParameters);
  return isOK;
}
/*****************************************************************************/
// List of components
vector<string> Bs2PhiKKSignal::PDFComponents()
{
  return componentnames;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKSignal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
  string compName = component->getComponentName();
  if(compName=="0") // should quicken things up slightly?
    return Evaluate(measurement);
  ReadDataPoint(measurement);
  if(mKK < 2 * Bs2PhiKKComponent::mK) return 0;
  // Evaluate the PDF at this point
  double evalRes = -1; // Important that this initalises to negative value
  if(compName=="interference")
  {
    evalRes = TotalDecayRate();
    for(auto comp : components)
      evalRes -= ComponentDecayRate(comp);
  }
  else
  {
    for(auto comp: components)
      if(compName.find(comp.GetName()) != string::npos)
        evalRes = ComponentDecayRate(comp, compName);
    // If none of the component names matched, assume total PDF
    if(evalRes<0) return Evaluate(measurement);
  }
  return evalRes;
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
  ReadDataPoint(measurement);
  if(mKK < 2 * Bs2PhiKKComponent::mK) return 0;
  return TotalDecayRate();
}
/*****************************************************************************/
void Bs2PhiKKSignal::ReadDataPoint(DataPoint* measurement)
{
  // Get values from the datapoint
  mKK       = measurement->GetObservable(mKKName     )->GetValue();
  phi       = measurement->GetObservable(phiName     )->GetValue();
  ctheta_1  = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2  = measurement->GetObservable(ctheta_2Name)->GetValue();
  CalculateAcceptance(); // Evaluate this here so we can check for negative acceptance
  // Check if the datapoint makes sense
  if(phi < -TMath::Pi() || phi > TMath::Pi()
       || ctheta_1 < -1 || ctheta_1 > 1
       || ctheta_2 < -1 || ctheta_2 > 1
       || acceptance < 0
  )
  {
    measurement->Print();
  }
  phi+=TMath::Pi();
}
/*****************************************************************************/
double Bs2PhiKKSignal::TotalDecayRate()
{
  map<bool,TComplex> TotalAmp;
  for(auto anti : {false,true})
  {
    TotalAmp[anti] = TComplex(0,0);
    for(auto comp : components)
      TotalAmp[anti] += anti ? comp.Amplitude(mKK, -phi, -ctheta_1, -ctheta_2) : comp.Amplitude(mKK, phi, ctheta_1, ctheta_2);
  }
  return TimeIntegratedDecayRate(TotalAmp[false],TotalAmp[true]);
}
double Bs2PhiKKSignal::ComponentDecayRate(Bs2PhiKKComponent& comp)
{
  return ComponentDecayRate(comp,"");
}
double Bs2PhiKKSignal::ComponentDecayRate(Bs2PhiKKComponent& comp, string option)
{
  return TimeIntegratedDecayRate(comp.Amplitude(mKK, phi, ctheta_1, ctheta_2, option),comp.Amplitude(mKK, -phi, -ctheta_1, -ctheta_2, option));
}
double Bs2PhiKKSignal::TimeIntegratedDecayRate(TComplex A, TComplex Abar)
{
  double GH = (2-dGsGs.value); // Actually Γ/2ΓH but who cares about an overall factor Γ?
  double GL = (2+dGsGs.value); // Actually Γ/2ΓL
//  return (A.Rho2())*(1/GL + 1/GH) + (TComplex::Power(TComplex::Conjugate(A),2)).Re()*(1/GL - 1/GH);
  double termone = (A.Rho2() + Abar.Rho2())*(1./GL + 1./GH);
  double termtwo = 2*(Abar*TComplex::Conjugate(A)).Re()*(1./GL - 1./GH);
  if(!std::isnan(termone) && !std::isnan(termtwo))
    return  (termone + termtwo) * p1stp3() * acceptance;
  else
    return 0;
}
/*****************************************************************************/
void Bs2PhiKKSignal::CalculateAcceptance()
{
  if(acceptance_moments) acceptance = acc_m->Evaluate(mKK*1000, phi, ctheta_1, ctheta_2);
  else if(acceptance_histogram) acceptance = acc_h->Eval({mKK*1000, abs(phi), abs(ctheta_1), abs(ctheta_2)});
  acceptance = acceptance > 0 ? acceptance : 1e-12;
}
/*****************************************************************************/
double Bs2PhiKKSignal::p1stp3()
{
//  _mKK/=1000.;
  const double mK   = Bs2PhiKKComponent::mK;//1000.;
  const double mBs  = Bs2PhiKKComponent::mBs;//1000.;
  const double mPhi = phimass.value;
  double pR = DPHelpers::daughterMomentum(mKK, mK,  mK);
  double pB = DPHelpers::daughterMomentum(mBs, mKK, mPhi);
  double pRpB = pR * pB;
  return TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  return -1;
}

