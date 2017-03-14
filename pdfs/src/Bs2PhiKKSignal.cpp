// Std Libraries
#include <stdexcept>
// ROOT Libraries
#include "TFile.h"
#include "TTree.h"
// RapidFit
#include "Bs2PhiKKSignal.h"
#include "DPHelpers.hh"
#include "PDFConfigurator.h"
#include "StringProcessing.h"
// GSL
#include <gsl/gsl_randist.h>

PDF_CREATOR( Bs2PhiKKSignal )
/*****************************************************************************/
// Constructor
Bs2PhiKKSignal::Bs2PhiKKSignal(PDFConfigurator* config) : Bs2PhiKK(config)
	, acceptance_moments((std::string)config->getConfigurationValue("CoefficientsFile") != "")
	, acceptance_histogram((std::string)config->getConfigurationValue("HistogramFile") != "")
{
	std::cout << "\nBuilding Bs → ϕ K+ K− signal PDF\n\n";
	std::string phiname = config->getConfigurationValue("phiname");
	thraccscale = Bs2PhiKK::PhysPar(config,"thraccscale");
	phimass = Bs2PhiKK::PhysPar(config,phiname+"_mass");
	dGsGs = Bs2PhiKK::PhysPar(config,"dGsGs");
	if(config->isTrue("convolve"))
	{
		for(std::string name: {"nsteps","nsigma"})
			mKKresconfig[name] = std::stod(config->getConfigurationValue("mKKres_"+name));
		mKKres_sigmazero = Bs2PhiKK::PhysPar(config,"mKKres_sigmazero");
	}
	std::vector<std::string> reslist = StringProcessing::SplitString(config->getConfigurationValue("resonances"), ' ');
	std::cout << "┏━━━━━━━━━━━━━━━┯━━━━━━━┯━━━━━━━━━━━━━━━┓\n";
	std::cout << "┃ Component     │ Spin  │ Lineshape     ┃\n";
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
Bs2PhiKKSignal::Bs2PhiKKSignal(const Bs2PhiKKSignal& copy)
	: BasePDF((BasePDF)copy)
	, Bs2PhiKK((Bs2PhiKK)copy)
	// Width splitting
	, dGsGs(copy.dGsGs)
	// Phi mass
	, phimass(copy.phimass)
	// threshold acceptance scale
	, thraccscale(copy.thraccscale)
	// mass resolution parameters
	, mKKres_sigmazero(copy.mKKres_sigmazero)
	, mKKresconfig(copy.mKKresconfig)
	// PDF components
	, components(copy.components)
	// Plotting components
	, componentnames(copy.componentnames)
	// Options
	, acceptance_moments(copy.acceptance_moments)
	, acceptance_histogram(copy.acceptance_histogram)
{
	if(acceptance_moments) acc_m = std::unique_ptr<LegendreMomentShape>(new LegendreMomentShape(*copy.acc_m));
	else if(acceptance_histogram) acc_h = copy.acc_h;
	Initialise();
}
/*****************************************************************************/
// Code common to the constructors
void Bs2PhiKKSignal::Initialise()
{
	// Enable numerical normalisation and disable caching
	this->SetNumericalNormalisation( true );
//	this->TurnCachingOff();
}
/*****************************************************************************/
// Build a component object from a passed option
Bs2PhiKKSignalComponent Bs2PhiKKSignal::ParseComponent(PDFConfigurator* config, std::string phiname, std::string option) const
{
	// Syntax: <resonance name>(<spin>,<lineshape>)
	// - the list is space-delimited: no extra spaces, please!
	// - resonance name must be alphanumeric
	// - spin must be a single digit 0, 1 or 2, or you have to implement the code for higher spins
	// - lineshape name must be 2 capital letters
	// - see Bs2PhiKKSignalComponent.cpp for implemented lineshapes (BW, FT, NR... but more can be added)
	// Example: "phi1020(1,BW)"
	std::string KKname = option.substr(0,option.find('('));
	int JKK = atoi(option.substr(option.find('(')+1,1).c_str());
	std::string lineshape = option.substr(option.find(',')+1,2);
	std::cout << "┃ " << KKname << "\t│    " << JKK << "\t│ " << lineshape << " shape \t┃\n";
	return Bs2PhiKKSignalComponent(config, phiname, KKname, JKK, lineshape);
}
/*****************************************************************************/
// Make the data point and parameter set
void Bs2PhiKKSignal::MakePrototypes()
{
	std::cout << "Bs2PhiKKSignal: making prototypes\n";
	// Make the DataPoint prototype
	MakePrototypeDataPoint(allObservables);
	// Make the parameter set
	std::vector<std::string> parameterNames;
	// Resonance parameters
	for(const auto& par: {dGsGs, phimass, thraccscale, mKKres_sigmazero})
		parameterNames.push_back(par.name);
	for(const auto& comp: components)
		for(std::string par: comp.second.GetPhysicsParameters())
			parameterNames.push_back(par);
	std::sort(parameterNames.begin(),parameterNames.end());
	parameterNames.erase(std::unique(parameterNames.begin(),parameterNames.end()),parameterNames.end());
	allParameters = ParameterSet(parameterNames);
}
// List of components
std::vector<std::string> Bs2PhiKKSignal::PDFComponents()
{
	// Avoid redundant plotting for single-component PDFs
	if(components.size() == 1) return {};
	return componentnames;
}
/*****************************************************************************/
// Set the physics parameters
bool Bs2PhiKKSignal::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
	UnsetCache();
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	for(auto* par: {&dGsGs, &phimass, &thraccscale, &mKKres_sigmazero})
		par->Update(&allParameters);
	for(auto& comp: components)
			comp.second.SetPhysicsParameters(&allParameters);
	return isOK;
}
/*Calculate the likelihood****************************************************/
// Evaluate a single component
double Bs2PhiKKSignal::EvaluateComponent(DataPoint* measurement, ComponentRef* component)
{
	std::string compName = component->getComponentName();
	// Verification
	if(compName.find_first_not_of("0123456789") == string::npos) return Evaluate(measurement); // If the component name is purely numerical, then we're being asked for the total PDF
	const Bs2PhiKK::datapoint_t datapoint = ReadDataPoint(measurement);
	if(!Bs2PhiKK::IsPhysicalDataPoint(datapoint)) throw std::out_of_range("Unphysical datapoint"); // Check if the point is above threshold and |cos(θ)| <= 1
	// Evaluation
	MsqFunc_t EvaluateMsq = compName=="interference"? &Bs2PhiKKSignal::InterferenceMsq : &Bs2PhiKKSignal::ComponentMsq; // Choose the function to calcualte the |M|²
	double MatrixElementSquared = mKKresconfig.empty()? (this->*EvaluateMsq)(datapoint,compName) : Convolve(EvaluateMsq,datapoint,compName); // Do the convolved or unconvolved |M|² calculation
	return Evaluate_Base(MatrixElementSquared, datapoint);
}
// Evaluate the entire PDF
double Bs2PhiKKSignal::Evaluate(DataPoint* measurement)
{
	// Verification
	const Bs2PhiKK::datapoint_t datapoint = ReadDataPoint(measurement);
	if(!Bs2PhiKK::IsPhysicalDataPoint(datapoint)) throw std::out_of_range("Unphysical datapoint");  // Check if the point is above threshold and |cos(θ)| <= 1
	// Evaluation
	double MatrixElementSquared = mKKresconfig.empty()? TotalMsq(datapoint) : Convolve(&Bs2PhiKKSignal::TotalMsq,datapoint,""); // Do the convolved or unconvolved |M|² calculation
	return Evaluate_Base(MatrixElementSquared, datapoint);
}
// The stuff common to both Evaluate() and EvaluateComponent()
double Bs2PhiKKSignal::Evaluate_Base(const double MatrixElementSquared, const Bs2PhiKK::datapoint_t& datapoint) const
{
	return MatrixElementSquared * p1stp3(datapoint[0]) * Acceptance(datapoint);
}
/*Calculate matrix elements***************************************************/
// Total |M|²: coherent sum of all amplitudes
double Bs2PhiKKSignal::TotalMsq(const Bs2PhiKK::datapoint_t& datapoint, const std::string& dummy) const
{
	(void)dummy;
	Bs2PhiKK::amplitude_t TotalAmp = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	for(const auto& comp : components)
	{
		Bs2PhiKK::amplitude_t CompAmp = comp.second.Amplitude(datapoint);
		if(std::isnan(CompAmp[0].real()) || std::isnan(CompAmp[0].imag())){ std::cerr << comp.first << " amplitude evaluates to " << CompAmp[0] << std::endl; std::exit(1);}
		TotalAmp[false] += CompAmp[false];
		TotalAmp[true] += CompAmp[true];
	}
	return TimeIntegratedMsq(TotalAmp);
}
// Single-component |M|²
double Bs2PhiKKSignal::ComponentMsq(const Bs2PhiKK::datapoint_t& datapoint, const std::string& compName) const
{
	return TimeIntegratedMsq(components.at(compName).Amplitude(datapoint, compName));
}
// Interference |M|²: difference between incoherent and coherent sums of components
double Bs2PhiKKSignal::InterferenceMsq(const Bs2PhiKK::datapoint_t& datapoint, const std::string& dummy) const
{
	(void)dummy;
	double MatrixElementSquared = TotalMsq(datapoint);
		for(const auto& comp : components)
			MatrixElementSquared -= TimeIntegratedMsq(comp.second.Amplitude(datapoint));
	return MatrixElementSquared;
}
// Take the amplitudes of the B and Bbar decay and return the time-integrated |M|²
double Bs2PhiKKSignal::TimeIntegratedMsq(const Bs2PhiKK::amplitude_t& Amp) const
{
	double GH = (2 - dGsGs.value); // Actually Γ/2ΓH but who cares about an overall factor Γ?
	double GL = (2 + dGsGs.value); // Actually Γ/2ΓL
	double termone = (std::norm(Amp[false]) + std::norm(Amp[true])) * (1./GL + 1./GH);
	double termtwo = 2*std::real(Amp[true] * std::conj(Amp[false])) * (1./GL - 1./GH);
	return termone + termtwo;
}
// Convolution of the matrix element function with a Gaussian... the slow integral way
double Bs2PhiKKSignal::Convolve(MsqFunc_t EvaluateMsq, const Bs2PhiKK::datapoint_t& datapoint, const std::string& compName) const
{
	const double nsigma = mKKresconfig.at("nsigma");
	const int nsteps = mKKresconfig.at("nsteps");
	const double resolution = std::sqrt(mKKres_sigmazero.value*(datapoint[0]-2*Bs2PhiKK::mK)); // Mass-dependent Gaussian width
	// If the integration region goes below threshold, don't do the convolution
	if(datapoint[0] - nsigma*resolution > 2*Bs2PhiKK::mK)
	{
		double Msq_conv = 0.;
		const double stepsize = 2.*nsigma*resolution/nsteps;
		// Integrate over range −nσ to +nσ
		for(double x = -nsigma*resolution; x < nsigma*resolution; x += stepsize)
			Msq_conv += gsl_ran_gaussian_pdf(x,resolution) * (this->*EvaluateMsq)({datapoint[0]-x,datapoint[1],datapoint[2],datapoint[3]},compName) * stepsize;
		return Msq_conv;
	}
	else
		return (this->*EvaluateMsq)(datapoint,compName);
}
/*Stuff that factors out of the time integral*********************************/
double Bs2PhiKKSignal::Acceptance(const Bs2PhiKK::datapoint_t& datapoint) const
{
	double acceptance;
	if(acceptance_moments)
	{
		// Get the shape from stored Legendre moments
		acceptance = acc_m->Evaluate(datapoint);
		// Multiply by a switch-on function
		acceptance *= std::erf(thraccscale.value*(datapoint[0]-2*Bs2PhiKK::mK));
//		acceptance *= std::tanh(thraccscale.value*(datapoint[0]-2*Bs2PhiKK::mK));
//		acceptance *= std::atan(thraccscale.value*(datapoint[0]-2*Bs2PhiKK::mK))*2.0/M_PI;
	}
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
double Bs2PhiKKSignal::p1stp3(const double& mKK) const
{
	double pR = DPHelpers::daughterMomentum(mKK, Bs2PhiKK::mK, Bs2PhiKK::mK);
	double pB = DPHelpers::daughterMomentum(Bs2PhiKK::mBs, mKK, phimass.value);
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

