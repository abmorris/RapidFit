// Self
#include "Bs2PhiKKSignal.h"
// Std Libraries
#include <iostream>
#include <stdexcept>
#include <complex>
// ROOT Libraries
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
	// Dependent variable names
	 mKKName     (config->getName("mKK"     ))
	,phiName     (config->getName("phi"     ))
	,ctheta_1Name(config->getName("ctheta_1"))
	,ctheta_2Name(config->getName("ctheta_2"))
	,acceptance_moments((std::string)config->getConfigurationValue("CoefficientsFile") != "")
	,acceptance_histogram((std::string)config->getConfigurationValue("HistogramFile") != "")
{
	std::cout << "\nBuilding Bs → ϕ K+ K− signal PDF\n\n";
	std::string phiname = config->getConfigurationValue("phiname");
	phimass = Bs2PhiKKComponent::PhysPar(config,phiname+"_mass");
	dGsGs = Bs2PhiKKComponent::PhysPar(config,"dGsGs");
	std::vector<std::string> reslist = StringProcessing::SplitString(config->getConfigurationValue("resonances"), ' ');
	std::cout << "┏━━━━━━━━━━━━━━━┯━━━━━━━┯━━━━━━━━━━━━━━━┓\n";
	std::cout << "┃ Component\t│ Spin\t│ Lineshape\t┃\n";
	std::cout << "┠───────────────┼───────┼───────────────┨\n";
	for(const auto& name: reslist)
	{
		if(name=="") continue;
		components[name] = ParseComponent(config,phiname,name);
		componentnames.push_back(name);
	}
	if(components.size() > 1) componentnames.push_back("interference");
	std::cout << "┗━━━━━━━━━━━━━━━┷━━━━━━━┷━━━━━━━━━━━━━━━┛" << std::endl;
	if(acceptance_moments) acc_m = std::unique_ptr<LegendreMomentShape>(new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile")));
	else if(acceptance_histogram)
	{
		TFile* histfile = TFile::Open(config->getConfigurationValue("HistogramFile").c_str());
		acc_h = std::shared_ptr<NDHist_Adaptive>(new NDHist_Adaptive(histfile));
		acc_h->LoadFromTree((TTree*)histfile->Get("AccTree"));
		acc_h->SetDimScales({1.,0.1,1.,1.}); // TODO: read this from the file
	}
	Initialise();
	MakePrototypes();
}
/*****************************************************************************/
// Copy constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(const Bs2PhiKKSignal& copy) : BasePDF( (BasePDF) copy)
	// Dependent variable names
	,mKKName     (copy.mKKName     )
	,phiName     (copy.phiName     )
	,ctheta_1Name(copy.ctheta_1Name)
	,ctheta_2Name(copy.ctheta_2Name)
	// Width splitting
	,dGsGs(copy.dGsGs)
	// Phi mass
	,phimass(copy.phimass)
	// PDF components
	,components(copy.components)
	,componentnames(copy.componentnames)
	// Plotting components
	// Options
	,acceptance_moments(copy.acceptance_moments)
	,acceptance_histogram(copy.acceptance_histogram)
{
	if(acceptance_moments) acc_m = std::unique_ptr<LegendreMomentShape>(new LegendreMomentShape(*copy.acc_m));
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
Bs2PhiKKComponent Bs2PhiKKSignal::ParseComponent(PDFConfigurator* config, std::string phiname, std::string option) const
{
	// Syntax: <resonance name>(<spin>,<lineshape>)
	// - the list is space-delimited: no extra spaces, please!
	// - resonance name must be alphanumeric
	// - spin must be a single digit 0, 1 or 2, or you have to implement the code for higher spins
	// - lineshape name must be 2 capital letters
	// - see Bs2PhiKKComponent.cpp for implemented lineshapes (BW, FT, NR... but more can be added)
	// Example: "phi1020(1,BW)"
	std::string KKname = option.substr(0,option.find('('));
	int JKK = atoi(option.substr(option.find('(')+1,1).c_str());
	std::string lineshape = option.substr(option.find(',')+1,2);
	std::cout << "┃ " << KKname << "\t│    " << JKK << "\t│ " << lineshape << " shape \t┃\n";
	return Bs2PhiKKComponent(config, phiname, KKname, JKK, lineshape);
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
	std::vector<std::string> parameterNames;
	// Resonance parameters
	parameterNames.push_back(dGsGs.name);
	parameterNames.push_back(phimass.name);
	for(const auto& comp: components)
		for(std::string par: comp.second.GetPhysicsParameters())
			parameterNames.push_back(par);
	std::sort(parameterNames.begin(),parameterNames.end());
	parameterNames.erase(std::unique(parameterNames.begin(),parameterNames.end()),parameterNames.end());
	for(const auto& name: parameterNames)
		std::cout << name << std::endl;
	allParameters = ParameterSet(parameterNames);
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKSignal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	dGsGs.Update(&allParameters);
	phimass.Update(&allParameters);
	for(auto& comp: components)
		comp.second.SetPhysicsParameters(&allParameters);
	return isOK;
}
/*****************************************************************************/
// List of components
std::vector<std::string> Bs2PhiKKSignal::PDFComponents()
{
	// Avoid redundant plotting for single-component PDFs
	if(components.size() == 1) return {};
	return componentnames;
}
/*****************************************************************************/
// Evaluate a single component
double Bs2PhiKKSignal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
	std::string compName = component->getComponentName();
	if(compName=="0") // should quicken things up slightly?
		return Evaluate(measurement);
	const std::array<double,4> datapoint = ReadDataPoint(measurement);
	// Evaluate the PDF at this point
	double MatrixElementSquared;
	if(compName=="interference")
	{
		MatrixElementSquared = TotalMsq(datapoint);
		for(const auto& comp : components)
			MatrixElementSquared -= TimeIntegratedMsq(comp.second.Amplitude(datapoint));
	}
	else
		MatrixElementSquared = TimeIntegratedMsq(components[compName].Amplitude(datapoint, compName)); // The name can also contain an option like "odd" or "even"
	double jacobian = p1stp3(datapoint[0]);
	double acceptance = Acceptance(datapoint);
	return MatrixElementSquared * jacobian * acceptance;
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
	const std::array<double,4> datapoint = ReadDataPoint(measurement);
	double MatrixElementSquared = TotalMsq(datapoint);
	double jacobian = p1stp3(datapoint[0]);
	double acceptance = Acceptance(datapoint);
	return MatrixElementSquared * jacobian * acceptance;
}
/*****************************************************************************/
std::array<double,4> Bs2PhiKKSignal::ReadDataPoint(DataPoint* measurement) const
{
	// Get values from the datapoint
	double mKK      = measurement->GetObservable(mKKName     )->GetValue();
	double phi      = measurement->GetObservable(phiName     )->GetValue();
	double ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
	double ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
	phi+=M_PI;
	return {mKK, phi, ctheta_1, ctheta_2};
}
/*****************************************************************************/
double Bs2PhiKKSignal::TotalMsq(const std::array<double,4>& datapoint) const
{
	std::array<std::complex<double>,2> TotalAmp = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	for(const auto& comp : components)
	{
		std::array<std::complex<double>,2> CompAmp = comp.second.Amplitude(datapoint);
		if(std::isnan(CompAmp[0].real()) || std::isnan(CompAmp[0].imag())){ std::cerr << comp.first << " amplitude evaluates to " << CompAmp[0] << std::endl; std::exit(1);}
		TotalAmp[false] += CompAmp[false];
		TotalAmp[true] += CompAmp[true];
	}
	return TimeIntegratedMsq(TotalAmp);
}
double Bs2PhiKKSignal::TimeIntegratedMsq(const std::array<std::complex<double>,2>& Amp) const
{
	double GH = (2 - dGsGs.value); // Actually Γ/2ΓH but who cares about an overall factor Γ?
	double GL = (2 + dGsGs.value); // Actually Γ/2ΓL
	double termone = (std::norm(Amp[false]) + std::norm(Amp[true])) * (1./GL + 1./GH);
	double termtwo = 2*std::real(Amp[true] * std::conj(Amp[false])) * (1./GL - 1./GH);
	return termone + termtwo;
}
/*****************************************************************************/
double Bs2PhiKKSignal::Acceptance(const std::array<double,4>& datapoint) const
{
	double acceptance;
	if(acceptance_moments)
		acceptance = acc_m->Evaluate(datapoint);
	else if(acceptance_histogram)
	{
		std::vector<double> absdatapoint;
		for(const auto& val: datapoint)
			absdatapoint.push_back(std::abs(val));
		acceptance = acc_h->Eval(absdatapoint);
	}
	if(std::isnan(acceptance)) std::cerr << "Acceptance evaluates to nan" << std::endl;
	return acceptance;
}
/*****************************************************************************/
double Bs2PhiKKSignal::p1stp3(const double& mKK) const
{
	const double mK   = Bs2PhiKKComponent::mK;
	if(mKK < 2*mK) return 0;
	const double mBs  = Bs2PhiKKComponent::mBs;
	const double mPhi = phimass.value;
	double pR = DPHelpers::daughterMomentum(mKK, mK,  mK);
	double pB = DPHelpers::daughterMomentum(mBs, mKK, mPhi);
	double pRpB = pR * pB;
	if(std::isnan(pRpB)) std::cerr << "p1stp3 evaluates to nan" << std::endl;
	return pRpB;
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
	(void)boundary;
	return -1;
}

