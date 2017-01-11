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
	std::vector<std::string> reslist = StringProcessing::SplitString( config->getConfigurationValue("resonances"), ' ' );
	std::cout << "┏━━━━━━━━━━━━━━━┯━━━━━━━┯━━━━━━━━━━━━━━━┓\n";
	std::cout << "┃ Component\t│ Spin\t│ Lineshape\t┃\n";
	std::cout << "┠───────────────┼───────┼───────────────┨\n";
	for(auto name: reslist)
	{
		if(name=="") continue;
		Bs2PhiKKComponent comp = std::move(ParseComponent(config,phiname,name));
		components.push_back(comp);
		componentnames.push_back(name);
	}
	std::cout << "┗━━━━━━━━━━━━━━━┷━━━━━━━┷━━━━━━━━━━━━━━━┛" << std::endl;
	if(componentnames.size() > 1)
		componentnames.push_back("interference");
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
	// Plotting components
	,componentnames(copy.componentnames)
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
	return std::move(Bs2PhiKKComponent(config, phiname, KKname, JKK, lineshape));
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
		for(std::string par: comp.GetPhysicsParameters())
			if(par!=phimass.name.Name()) parameterNames.push_back(par);
//	allParameters = *( new ParameterSet(parameterNames) );
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
		comp.SetPhysicsParameters(&allParameters);
	return isOK;
}
/*****************************************************************************/
// List of components
std::vector<std::string> Bs2PhiKKSignal::PDFComponents()
{
	// Avoid redundant plotting for single-component PDFs
	if(componentnames.size()>1)
		return componentnames;
	else
		return {};
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
	double MatrixElementSquared = -1; // Important that this initalises to negative value
	if(compName=="interference")
	{
		MatrixElementSquared = TotalMsq(datapoint);
		for(const auto& comp : components)
			MatrixElementSquared -= TimeIntegratedMsq(comp.Amplitude(datapoint));
	}
	else
	{
		for(const auto& comp: components)
			if(compName.find(comp.GetName()) != std::string::npos)
				MatrixElementSquared = TimeIntegratedMsq(comp.Amplitude(datapoint, compName));
		// If none of the component names matched, assume total PDF
		if(MatrixElementSquared<0) return Evaluate(measurement);
	}
	double jacobian = p1stp3(datapoint[0]);
	double acceptance = Acceptance(datapoint);
	double returnVal = MatrixElementSquared * jacobian * acceptance;
	return returnVal;
}
/*****************************************************************************/
// Calculate the function value
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
	const std::array<double,4> datapoint = ReadDataPoint(measurement);
	double MatrixElementSquared = TotalMsq(datapoint);
	double jacobian = p1stp3(datapoint[0]);
	double acceptance = Acceptance(datapoint);
	double returnVal = MatrixElementSquared * jacobian * acceptance;
	return returnVal;
}
/*****************************************************************************/
std::array<double,4> Bs2PhiKKSignal::ReadDataPoint(DataPoint* measurement) const
{
	// Get values from the datapoint
	double mKK       = measurement->GetObservable(mKKName     )->GetValue();
	double phi       = measurement->GetObservable(phiName     )->GetValue();
	double ctheta_1  = measurement->GetObservable(ctheta_1Name)->GetValue();
	double ctheta_2  = measurement->GetObservable(ctheta_2Name)->GetValue();
	// Check if the datapoint makes sense
	if(std::abs(phi) > M_PI || std::abs(ctheta_1) > 1 || std::abs(ctheta_2) > 1)
		std::cerr << "Datapoint outside phasespace!" << std::endl;
	phi+=M_PI;
	return {mKK, phi, ctheta_1, ctheta_2};
}
/*****************************************************************************/
double Bs2PhiKKSignal::TotalMsq(const std::array<double,4>& datapoint) const
{
	std::array<std::complex<double>,2> TotalAmp = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	for(const auto& comp : components)
	{
		std::array<std::complex<double>,2> CompAmp = comp.Amplitude(datapoint);
		TotalAmp[false] += CompAmp[false];
		TotalAmp[true] += CompAmp[true];
	}
	return TimeIntegratedMsq(TotalAmp);
}
double Bs2PhiKKSignal::TimeIntegratedMsq(const std::array<std::complex<double>,2>& Amp) const
{
	double GH = (2-dGsGs.value); // Actually Γ/2ΓH but who cares about an overall factor Γ?
	double GL = (2+dGsGs.value); // Actually Γ/2ΓL
	double termone = (std::norm(Amp[false]) + std::norm(Amp[true]))*(1./GL + 1./GH);
	double termtwo = 2*std::real(Amp[true]*std::conj(Amp[false]))*(1./GL - 1./GH);
	if(!std::isnan(termone) && !std::isnan(termtwo))
		return termone + termtwo;
	else
		return 0;
}
/*****************************************************************************/
double Bs2PhiKKSignal::Acceptance(const std::array<double,4>& datapoint) const
{
	double acceptance;
	if(acceptance_moments)
		acceptance = acc_m->Evaluate(datapoint);
	else if(acceptance_histogram)
		acceptance = acc_h->Eval(std::vector<double>(datapoint.begin(),datapoint.end()));
	return std::isnan(acceptance) ? 0 : acceptance;
}
/*****************************************************************************/
double Bs2PhiKKSignal::p1stp3(const double& mKK) const
{
	const double mK   = Bs2PhiKKComponent::mK;
	const double mBs  = Bs2PhiKKComponent::mBs;
	const double mPhi = phimass.value;
	double pR = DPHelpers::daughterMomentum(mKK, mK,  mK);
	double pB = DPHelpers::daughterMomentum(mBs, mKK, mPhi);
	double pRpB = pR * pB;
	return std::isnan(pRpB) ? 0 : pRpB; // Protect against divide-by-zero
}
/*****************************************************************************/
double Bs2PhiKKSignal::Normalisation(PhaseSpaceBoundary* boundary)
{
	(void)boundary;
	return -1;
}

