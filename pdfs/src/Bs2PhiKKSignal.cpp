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
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
// RapidFit
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
  // Variables for resonance parameters
  , fzeroMass(Bs2PhiKKComponent::mfzero)
  , fzerogpipi(Bs2PhiKKComponent::gpipi)
  , fzeroRg(Bs2PhiKKComponent::Rg)
  , phiMass(Bs2PhiKKComponent::mphi)
  , phiWidth(Bs2PhiKKComponent::wphi)
  , ftwoMass(Bs2PhiKKComponent::mftwo)
  , ftwoWidth(Bs2PhiKKComponent::wftwo)
  // Resonance parameter names
  , fzeroMassName(config->getName("fzeroMass"))
  , fzerogpipiName(config->getName("fzerogpipi"))
  , fzeroRgName(config->getName("fzeroRg"))
  , phiMassName(config->getName("phiMass"))
  , phiWidthName(config->getName("phiWidth"))
  , ftwoMassName(config->getName("ftwoMass"))
  , ftwoWidthName(config->getName("ftwoWidth"))
  // Width splitting
  , dGsGs(0.0)
  // Variables for physics amplitudes
  // Scaling factors for each component
  , NonResFrac(0.)
  , SwaveFrac(0.)
  , PwaveFrac(0.)
  , DwaveFrac(0.)
  // Magnitude-squared of amplitudes
  , ASzerosq(0.)
  , APperpsq(0.)
  , APzerosq(0.)
  , APparasq(0.)
  , ADperpsq(0.)
  , ADzerosq(0.)
  , ADparasq(0.)
  // Phases of amplitudes
  , deltaSzero(0.)
  , deltaPperp(0.)
  , deltaPzero(0.)
  , deltaPpara(0.)
  , deltaDperp(0.)
  , deltaDzero(0.)
  , deltaDpara(0.)
  // Width splitting
  , dGsGsName(config->getName("dGsGs"))
  // Fit parameter names for physics amplitudes
  // Scaling factors for each component
  , NonResFracName(config->getName("NonResFrac"))
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
  , floatResPars(config->isTrue("FloatResPars"))
  , acceptance_moments((string)config->getConfigurationValue("CoefficientsFile") != "")
  , acceptance_histogram((string)config->getConfigurationValue("HistogramFile") != "")
{
  string RBs_str = config->getConfigurationValue("RBs");
  RBs = std::atof(RBs_str.c_str());
  string RKK_str = config->getConfigurationValue("RKK");
  RKK = std::atof(RKK_str.c_str());
  MakePrototypes();
  if(acceptance_histogram)
  {
    // The goal here is to get a TKDTreeID object to define the binning scheme, and to fill the bin content array from a tree that's sitting in the same file.
    // Much of this is re-implement from my FourDHist_Adaptive class in the BsPhiKK project
    TFile* accfile = TFile::Open(config->getConfigurationValue("HistogramFile").c_str());
    TIter next = accfile->GetListOfKeys();
    TKey* key;
    TTree* acctree;
    // Look for any TTree and TKDTree object in the file
    cout << "Looking for TTree and TKDTree objects in the file " << accfile->GetName() << endl;
    while((key = (TKey*)next()))
    {
      if(strcmp(key->GetClassName(),"TKDTree<int,double>")==0)
        accbinner = (TKDTreeID*)key->ReadObj();
      if(strcmp(key->GetClassName(),"TTree")==0)
        acctree = (TTree*)key->ReadObj();
    }
    // Throw a wobbly if you can't find either
    if(accbinner == NULL)
      throw runtime_error("Couldn't find a TKDTree object in the file");
    if(acctree == NULL)
      throw runtime_error("Couldn't find a TTree object in the file");
    // Load the binned values of acceptance
    cout << "Loading the binned acceptance values from the tree" << endl;
    int nbins = accbinner->GetNNodes()+1;
    double content;
    acctree->SetBranchAddress("content",&content);
    for(int ibin = 0; ibin < nbins; ibin++)
    {
      acctree->GetEntry(ibin);
      accbincontent.push_back(content);
    }
  }
  bool drawAll = config->isTrue("DrawAll");
  if(drawAll || config->isTrue("DrawNonRes"))
  {
    componentlist.push_back("NonRes");
  }
  if(drawAll || config->isTrue("DrawSwave"))
  {
    componentlist.push_back("Swave");
  }
  if(config->isTrue("DrawPwave"))
  {
    componentlist.push_back("Pwave");
  }
  if(drawAll || config->isTrue("DrawPwaveComponents"))
  {
    componentlist.push_back("Pwave-even");
    componentlist.push_back("Pwave-odd");
  }
  if(config->isTrue("DrawDwave"))
  {
    componentlist.push_back("Dwave");
  }
  if(drawAll || config->isTrue("DrawDwaveComponents"))
  {
    componentlist.push_back("Dwave-even");
    componentlist.push_back("Dwave-odd");
  }
  if(drawAll || config->isTrue("DrawInterference"))
  {
    componentlist.push_back("interference");
  }
  if(acceptance_moments) acc_m = new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile"));
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
  // Variables for resonance parameters
  , fzeroMass(copy.fzeroMass)
  , fzerogpipi(copy.fzerogpipi)
  , fzeroRg(copy.fzeroRg)
  , phiMass(copy.phiMass)
  , phiWidth(copy.phiWidth)
  , ftwoMass(copy.ftwoMass)
  , ftwoWidth(copy.ftwoWidth)
  // Resonance parameter names
  , fzeroMassName(copy.fzeroMassName)
  , fzerogpipiName(copy.fzerogpipiName)
  , fzeroRgName(copy.fzeroRgName)
  , phiMassName(copy.phiMassName)
  , phiWidthName(copy.phiWidthName)
  , ftwoMassName(copy.ftwoMassName)
  , ftwoWidthName(copy.ftwoWidthName)
  // Width splitting
  , dGsGs(copy.dGsGs)
  // Variables for physics amplitudes
  // Scaling factors for each component
  , NonResFrac(copy.NonResFrac)
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
  , deltaSzero(copy.deltaSzero)
  , deltaPperp(copy.deltaPperp)
  , deltaPzero(copy.deltaPzero)
  , deltaPpara(copy.deltaPpara)
  , deltaDperp(copy.deltaDperp)
  , deltaDzero(copy.deltaDzero)
  , deltaDpara(copy.deltaDpara)
  // Width splitting
  , dGsGsName(copy.dGsGsName)
  // Fit parameter names for physics amplitudes
  // Scaling factors for each component
  , NonResFracName(copy.NonResFracName)
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
  , componentlist(copy.componentlist)
  // Options
  , floatResPars(copy.floatResPars)
  , acceptance_moments(copy.acceptance_moments)
  , RBs(copy.RBs)
  , RKK(copy.RKK)
  , acceptance_histogram(copy.acceptance_histogram)
{
  if(acceptance_moments) acc_m = new LegendreMomentShape(*copy.acc_m);
  if(acceptance_histogram)
  {
    accbinner = copy.accbinner;
    accbincontent = copy.accbincontent;
  }
  Initialise();
}
/*****************************************************************************/
// Destructor
Bs2PhiKKSignal::~Bs2PhiKKSignal()
{
  for(auto comp : components)
    delete comp;
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKSignal::Initialise()
{
  // Initialise the signal components
  NonRes = new Bs2PhiKKComponent(0, 0, 0, "NR", RBs, RKK);
  Swave  = new Bs2PhiKKComponent(0, 0, 0, "FT", RBs, RKK);
  Pwave  = new Bs2PhiKKComponent(1, 0, 0, "BW", RBs, RKK);
  Dwave  = new Bs2PhiKKComponent(2, 0, 0, "BW", RBs, RKK);
  components.push_back(NonRes);
  components.push_back(Swave);
  components.push_back(Pwave);
  components.push_back(Dwave);
  SetResonanceParameters();
  SetComponentAmplitudes();
  // Enable numerical normalisation and disable caching
  this->SetNumericalNormalisation( true );
  this->TurnCachingOff();
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
  // Resonance parameters
  if(floatResPars)
  {
    parameterNames.push_back(fzeroMassName);
    parameterNames.push_back(fzerogpipiName);
    parameterNames.push_back(fzeroRgName);
    parameterNames.push_back(phiMassName);
    parameterNames.push_back(phiWidthName);
    parameterNames.push_back(ftwoMassName);
    parameterNames.push_back(ftwoWidthName);
  }
  parameterNames.push_back(dGsGsName);
  // Scaling factors for each component
  parameterNames.push_back(NonResFracName);
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
  // Resonance parameters
  if(floatResPars)
  {
    fzeroMass  = allParameters.GetPhysicsParameter(fzeroMassName )->GetValue();
    fzerogpipi = allParameters.GetPhysicsParameter(fzerogpipiName)->GetValue();
    fzeroRg    = allParameters.GetPhysicsParameter(fzeroRgName   )->GetValue();
    phiMass    = allParameters.GetPhysicsParameter(phiMassName   )->GetValue();
    phiWidth   = allParameters.GetPhysicsParameter(phiWidthName  )->GetValue();
    ftwoMass   = allParameters.GetPhysicsParameter(ftwoMassName  )->GetValue();
    ftwoWidth  = allParameters.GetPhysicsParameter(ftwoWidthName )->GetValue();
    SetResonanceParameters();
  }
  dGsGs   = allParameters.GetPhysicsParameter(dGsGsName  )->GetValue();
  // Scaling factors for each component
  NonResFrac = allParameters.GetPhysicsParameter(NonResFracName)->GetValue();
  SwaveFrac  = allParameters.GetPhysicsParameter(SwaveFracName )->GetValue();
  PwaveFrac  = allParameters.GetPhysicsParameter(PwaveFracName )->GetValue();
  DwaveFrac  = allParameters.GetPhysicsParameter(DwaveFracName )->GetValue();
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
  string compName = component->getComponentName();
  if(compName=="0") // should quicken things up slightly?
    return Evaluate(measurement);
  ReadDataPoint(measurement);
  if(mKK < 2 * Bs2PhiKKComponent::mK) return 0;
  // Evaluate the PDF at this point
  double absMatElSq = 0; // square of matrix element: |M|^2
  if(compName=="NonRes")
    absMatElSq = ComponentDecayRate(NonRes);
  else if(compName=="Swave")
    absMatElSq = ComponentDecayRate(Swave);
  else if(compName=="Pwave")
    absMatElSq = ComponentDecayRate(Pwave);
  else if(compName=="Pwave-odd")
    absMatElSq = ComponentDecayRate(Pwave, "odd");
  else if(compName=="Pwave-even")
    absMatElSq = ComponentDecayRate(Pwave, "even");
  else if(compName=="Dwave")
    absMatElSq = ComponentDecayRate(Dwave);
  else if(compName=="Dwave-odd")
    absMatElSq = ComponentDecayRate(Dwave, "odd");
  else if(compName=="Dwave-even")
    absMatElSq = ComponentDecayRate(Dwave, "even");
  else if(compName=="interference")
  {
    absMatElSq = TotalDecayRate();
    for(auto comp : components)
      absMatElSq -= ComponentDecayRate(comp);
  }
  else
    return Evaluate(measurement);
  return absMatElSq * p1stp3(mKK) * Acceptance();
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
  ReadDataPoint(measurement);
  if(mKK < 2 * Bs2PhiKKComponent::mK) return 0;
  return TotalDecayRate() * p1stp3(mKK) * Acceptance();
}
/*****************************************************************************/
void Bs2PhiKKSignal::ReadDataPoint(DataPoint* measurement)
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
    for(auto comp : components)
      comp->Print();
  }
}
/*****************************************************************************/
double Bs2PhiKKSignal::TotalDecayRate()
{
  std::map<bool,TComplex> TotalAmp;
  for(auto anti : {false,true})
  {
    TotalAmp[anti];
    for(auto comp : components)
      TotalAmp[anti] += anti ? comp->Amplitude(mKK, -phi, -ctheta_1, -ctheta_2) : comp->Amplitude(mKK, phi, ctheta_1, ctheta_2);
  }
  return TimeIntegratedDecayRate(TotalAmp[false],TotalAmp[true]);
}
double Bs2PhiKKSignal::ComponentDecayRate(Bs2PhiKKComponent* comp)
{
  return ComponentDecayRate(comp,"");
}
double Bs2PhiKKSignal::ComponentDecayRate(Bs2PhiKKComponent* comp, string option)
{
  return TimeIntegratedDecayRate(comp->Amplitude(mKK, phi, ctheta_1, ctheta_2, option),comp->Amplitude(mKK, -phi, -ctheta_1, -ctheta_2, option));
}
double Bs2PhiKKSignal::TimeIntegratedDecayRate(TComplex A, TComplex Abar)
{
  double GH = (2-dGsGs); // Actually Γ/2ΓH but who cares about an overall factor Γ?
  double GL = (2+dGsGs); // Actually Γ/2ΓL
//  return (A.Rho2())*(1/GL + 1/GH) + (TComplex::Power(TComplex::Conjugate(A),2)).Re()*(1/GL - 1/GH);
  return (A.Rho2() + Abar.Rho2())*(1/GL + 1/GH) + 2*(Abar*TComplex::Conjugate(A)).Re()*(1/GL - 1/GH);
}
/*****************************************************************************/
void Bs2PhiKKSignal::SetComponentAmplitudes()
{
  TComplex      ASzero(0, 0),
  APperp(0, 0), APzero(0, 0), APpara(0, 0),
  ADperp(0, 0), ADzero(0, 0), ADpara(0, 0);
  // S-wave complex amplitude
  ASzero = TComplex(sqrt(ASzerosq), deltaSzero, true) * sqrt(SwaveFrac);
  // P-wave complex amplitudes
  APperp = TComplex(sqrt(APperpsq), deltaPperp, true) * sqrt(PwaveFrac);
  APzero = TComplex(sqrt(APzerosq), deltaPzero, true) * sqrt(PwaveFrac);
  APpara = TComplex(sqrt(APparasq), deltaPpara+TMath::Pi(), true) * sqrt(PwaveFrac);
  TComplex APplus  = (APpara + APperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
  TComplex APminus = (APpara - APperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
  // D-wave complex amplitudes
  ADperp = TComplex(sqrt(ADperpsq), deltaDperp, true) * sqrt(DwaveFrac);
  ADzero = TComplex(sqrt(ADzerosq), deltaDzero, true) * sqrt(DwaveFrac);
  ADpara = TComplex(sqrt(ADparasq), deltaDpara+TMath::Pi(), true) * sqrt(DwaveFrac);
  TComplex ADplus  = (ADpara + ADperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
  TComplex ADminus = (ADpara - ADperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
  // Set the S-wave helicity amplitude
  NonRes->SetHelicityAmplitudes(0, sqrt(NonResFrac), 0);
  // Set the S-wave helicity amplitude
  Swave->SetHelicityAmplitudes(0, ASzero);
  // Set the P-wave helicity amplitudes
  Pwave->SetHelicityAmplitudes(0, APminus);
  Pwave->SetHelicityAmplitudes(1, APzero );
  Pwave->SetHelicityAmplitudes(2, APplus );
  // Set the D-wave helicity amplitudes
  Dwave->SetHelicityAmplitudes(0, ADminus);
  Dwave->SetHelicityAmplitudes(1, ADzero );
  Dwave->SetHelicityAmplitudes(2, ADplus );
}
void Bs2PhiKKSignal::SetResonanceParameters()
{
  Swave->SetMassCouplings(fzeroMass, fzerogpipi, fzeroRg * fzerogpipi);
  Pwave->SetMassWidth(phiMass, phiWidth);
  Dwave->SetMassWidth(ftwoMass, ftwoWidth);
}
/*****************************************************************************/
double Bs2PhiKKSignal::Acceptance()
{
  double acceptance;
  if(acceptance_moments) acceptance = acc_m->Evaluate(mKK, phi, ctheta_1, ctheta_2);
  else if(acceptance_histogram)
  {
    double point[] = {mKK*1e3,phi,ctheta_1,ctheta_2};
    int bin = accbinner->FindNode(point) - accbinner->GetNNodes();
    // Return the appropriate value for this bin
    acceptance = accbincontent[bin];
  }
  else acceptance = 1.;
  return acceptance>0 ? acceptance : 1e-12;
}
/*****************************************************************************/
double Bs2PhiKKSignal::p1stp3(double _mKK)
{
//  _mKK/=1000.;
  const double mK   = Bs2PhiKKComponent::mK;//1000.;
  const double mBs  = Bs2PhiKKComponent::mBs;//1000.;
  const double mPhi = Bs2PhiKKComponent::mphi;//1000.;
  double pR = DPHelpers::daughterMomentum(_mKK, mK, mK);
  double pB = DPHelpers::daughterMomentum(mBs, _mKK, mPhi);
  double pRpB = pR * pB;
  return TMath::IsNaN(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
  (void)boundary;
  return -1;
}

