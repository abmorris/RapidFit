// Self
#include "Bs2PhiKKSignalComponent.h"
// RapidFit Dalitz Plot Libraries
#include "DPBWResonanceShape.hh"
#include "DPFlatteShape.hh"
#include "DPNonresonant.hh"
#include "DPComplexSpline.hh"
#include "DPHelpers.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"
// Constructor
Bs2PhiKKSignalComponent::Bs2PhiKKSignalComponent(PDFConfigurator* config, std::string _phiname, std::string KKname, int _JKK, std::string _lineshape)
	: phimass(Bs2PhiKK::PhysPar(config,_phiname+"_mass"))
	, JKK(_JKK)
	, lineshape(_lineshape)
	, Bsbarrier(DPBarrierFactor(std::abs(JKK-1), 1.0, 0)) // Min L_B approximation
	, KKbarrier(DPBarrierFactor(JKK, 3.0, 0))
{
	if(lineshape != "SP")
		fraction = Bs2PhiKK::PhysPar(config,KKname+"_fraction");
	// Barrier factors
	if(lineshape != "NR")
	{
		BsBFradius = Bs2PhiKK::PhysPar(config,"BsBFradius");
		KKBFradius = Bs2PhiKK::PhysPar(config,"KKBFradius");
	}
	// KK resonance parameters
	for(string KKpar: Bs2PhiKK::LineShapeParameterNames(KKname,lineshape))
		KKpars.push_back(Bs2PhiKK::PhysPar(config,KKpar));
	if(KKpars.empty() && lineshape!="NR")
	{
		lineshape = "NR";
		std::cerr << "Bs2PhiKKSignalComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant." << std::endl;
	}
	// Helicity amplitude observables
	int lambda_max = std::min(1, JKK); // Maximum helicity
	std::vector<int> helicities;
	bool usephi = false;
	for(const auto& name: config->GetPhaseSpaceBoundary()->GetAllNames())
		if(name == "phi")
			usephi = true;
	if(lineshape != "NR" && lineshape != "SP")
		for(int lambda = usephi ? -lambda_max : 0; lambda <= lambda_max; lambda++)
			helicities.push_back(lambda);
	long unsigned int n = helicities.size();
	std::vector<std::string> phasesuffices;
	std::vector<std::string> magsqsuffices;
	switch(n)
	{
		case 0:
			break;
		case 1:
			phasesuffices = {"deltazero"};
			break;
		case 2:
			magsqsuffices = {"Azerosq"};
			phasesuffices = {"deltazero","deltaplus"};
			break;
		case 3:
			magsqsuffices = {"Azerosq","Aplussq"};
			phasesuffices = {"deltazero","deltaplus","deltaminus"};
			break;
		default:
			std::cerr << "Bs2PhiKKSignalComponent can't handle this many helicities: " << n << std::endl;
			std::exit(-1);
			break;
	}
	for(std::string suffix: magsqsuffices)
		magsqs.push_back(Bs2PhiKK::PhysPar(config,KKname+"_"+suffix));
	for(std::string suffix: phasesuffices)
		phases.push_back(Bs2PhiKK::PhysPar(config,KKname+"_"+suffix));
	// Make the helicity amplitude vector
	for(const auto& lambda: helicities)
		Ahel[lambda] = std::polar<double>(std::sqrt(1. / (double)n), 0);
	Initialise();
}
// Copy constructor
Bs2PhiKKSignalComponent::Bs2PhiKKSignalComponent(const Bs2PhiKKSignalComponent& other)
	// Floatable parameters
	: fraction(other.fraction)
	, Ahel(other.Ahel)
	, magsqs(other.magsqs)
	, phases(other.phases)
	, phimass(other.phimass)
	, KKpars(other.KKpars)
	, BsBFradius(other.BsBFradius)
	, KKBFradius(other.KKBFradius)
	// Fixed parameters
	, JKK(other.JKK)
	// DP objects
	, Bsbarrier(other.Bsbarrier)
	, KKbarrier(other.KKbarrier)
	// Options
	, lineshape(other.lineshape)
{
	Initialise();
}
/*****************************************************************************/
void Bs2PhiKKSignalComponent::Initialise()
{
	// Breit Wigner
	if(lineshape=="BW")
		KKLineShape = std::make_unique<DPBWResonanceShape>(DPBWResonanceShape(KKpars[0].value, KKpars[1].value, JKK, Bs2PhiKK::mK, Bs2PhiKK::mK, KKpars[2].value));
	// Flatte
	else if(lineshape=="FT")
		KKLineShape = std::make_unique<DPFlatteShape>(DPFlatteShape(KKpars[0].value, KKpars[1].value, Bs2PhiKK::mpi, Bs2PhiKK::mpi, KKpars[1].value*KKpars[2].value, Bs2PhiKK::mK, Bs2PhiKK::mK));
	else if(lineshape=="SP")
	{
		std::vector<double> breakpoints;
		for(unsigned i = 0; i < KKpars.size(); i+=3)
			breakpoints.push_back(KKpars[i].value);
		KKLineShape = std::make_unique<DPComplexSpline>(DPComplexSpline(breakpoints));
	}
	else
		KKLineShape = std::make_unique<DPNonresonant>(DPNonresonant());
	// Build the barrier factor and Wigner function objects
	wignerPhi = std::make_unique<DPWignerFunctionJ1>(DPWignerFunctionJ1());
	switch (JKK) // I hate this but I'd rather it just worked...
	{
		case 0:
			wignerKK  = std::make_unique<DPWignerFunctionJ0>(DPWignerFunctionJ0());
			break;
		case 1:
			wignerKK  = std::make_unique<DPWignerFunctionJ1>(DPWignerFunctionJ1());
			break;
		case 2:
			wignerKK  = std::make_unique<DPWignerFunctionJ2>(DPWignerFunctionJ2());
			break;
		default:
			wignerKK  = std::make_unique<DPWignerFunctionGeneral>(DPWignerFunctionGeneral(JKK)); // This should only happen for the rho_3 (1690)
			break;
	}
	UpdateAmplitudes();
	Bs2PhiKK::UpdateLineshape(lineshape, *KKLineShape, KKpars);
}
/*****************************************************************************/
// Angular part of the amplitude
std::complex<double> Bs2PhiKKSignalComponent::F(const int lambda, const Bs2PhiKK::datapoint_t& datapoint) const
{
	std::complex<double> returnval;
	const double d_phi = wignerPhi->function(datapoint.at(Bs2PhiKK::_ctheta_1_), lambda, 0);
	const double d_KK = wignerKK->function(datapoint.at(Bs2PhiKK::_ctheta_2_), lambda, 0);
	returnval = d_phi * d_KK;
	if(datapoint.find(Bs2PhiKK::_phi_) != datapoint.end())
		returnval *= std::polar<double>(1, lambda*datapoint.at(Bs2PhiKK::_phi_));
	return returnval;
}
std::complex<double> Bs2PhiKKSignalComponent::AngularPart(const Bs2PhiKK::datapoint_t& datapoint) const
{
	std::complex<double> angularPart(0, 0);
	if(Ahel.empty() || datapoint.find(Bs2PhiKK::_ctheta_1_) == datapoint.end() || datapoint.find(Bs2PhiKK::_ctheta_2_) == datapoint.end())
		return std::complex<double>(1, 0); // Either nonresonant or helicity angles not present
	else if(lineshape == "SP") // Phase already accounted-for by complex spline, so just return generic S-wave shape
		angularPart = F(0, datapoint);
	else
		for(const auto& A : Ahel)
		{
			int lambda = A.first;
			if(Ahel.size() == 3)
			angularPart += A.second * F(lambda, datapoint);
		}
	return angularPart;
}
// Orbital and barrier factor
double Bs2PhiKKSignalComponent::OFBF(const double mKK) const
{
	if(mKK < 2*Bs2PhiKK::mK)
		return 0;
	if(lineshape=="NR" || lineshape == "SP")
		return 1;
	// Momenta
	double pBs = DPHelpers::daughterMomentum(Bs2PhiKK::mBs, phimass.value, mKK);
	double pKK = DPHelpers::daughterMomentum(mKK, Bs2PhiKK::mK, Bs2PhiKK::mK);
	// Orbital factor
	double orbitalFactor = std::pow(pBs/Bs2PhiKK::mBs, std::abs(JKK-1))* // Min L_B approximation
	                       std::pow(pKK/KKpars[0].value, JKK);
	// Barrier factors
	double barrierFactor = Bsbarrier.barrier(pBs)*
	                       KKbarrier.barrier(pKK);
	return orbitalFactor*barrierFactor;
}
// Mass-dependent part
std::complex<double> Bs2PhiKKSignalComponent::MassPart(const double mKK) const
{
	std::complex<double> massPart = KKLineShape->massShape(mKK);
	if(lineshape != "SP")
	{
		massPart *= fraction.value * OFBF(mKK);
	}
	return massPart;
}
// The full amplitude.
Bs2PhiKK::amplitude_t Bs2PhiKKSignalComponent::Amplitude(const Bs2PhiKK::datapoint_t& datapoint) const
{
	Bs2PhiKK::amplitude_t angularPart = {AngularPart(datapoint), AngularPart(Bs2PhiKK::Parity(datapoint))};
	std::complex<double> massPart(1, 0);
	if(datapoint.find(Bs2PhiKK::_mKK_) != datapoint.end())
		massPart = MassPart(datapoint.at(Bs2PhiKK::_mKK_));
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
// The full amplitude with an option.
Bs2PhiKK::amplitude_t Bs2PhiKKSignalComponent::Amplitude(const Bs2PhiKK::datapoint_t& datapoint, const std::string option) const
{
	if(Ahel.empty() || option == "" || option.find("odd") == std::string::npos || option.find("even") == std::string::npos || datapoint.find(Bs2PhiKK::_phi_) == datapoint.end())
		return Amplitude(datapoint);
	// Angular part
	Bs2PhiKK::amplitude_t angularPart = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	std::complex<double> Aperp = std::polar<double>(std::sqrt(magsqs[0].value),phases[0].value);
	std::complex<double> Apara = std::polar<double>(std::sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value);
	// Temporary helicity amplitudes
	std::vector<std::complex<double>> HelAmp;
	// CP-odd component
	if(option.find("odd") != std::string::npos)
		HelAmp = {-Aperp/std::sqrt(2.), std::complex<double>(0, 0), Aperp/std::sqrt(2.)};
	// CP-even component
	else
		HelAmp = {Apara/std::sqrt(2.), Ahel.at(0), Apara/std::sqrt(2.)};
	for(const auto& A : Ahel)
	{
		int lambda = A.first;
		angularPart[false] += HelAmp[lambda+1] * F(lambda, datapoint);
		angularPart[true] += HelAmp[lambda+1] * F(lambda, Bs2PhiKK::Parity(datapoint));
	}
	std::complex<double> massPart(1, 0);
	if(datapoint.find(Bs2PhiKK::_mKK_) != datapoint.end())
		massPart = MassPart(datapoint.at(Bs2PhiKK::_mKK_));
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
/*****************************************************************************/
// Update everything from the parameter set
void Bs2PhiKKSignalComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
	// Update the parameters objects first
	if(lineshape != "SP")
		fraction.Update(fitpars);
	phimass.Update(fitpars);
	if(lineshape != "NR")
	{
		KKBFradius.Update(fitpars);
		BsBFradius.Update(fitpars);
	}
	for(auto* set: {&magsqs,&phases,&KKpars})
		for(auto& par: *set)
			par.Update(fitpars);
	UpdateBarriers();
	UpdateAmplitudes();
	Bs2PhiKK::UpdateLineshape(lineshape, *KKLineShape, KKpars);
}
void Bs2PhiKKSignalComponent::UpdateAmplitudes()
{
	// Update the helicity amplitudes
	if(Ahel.empty())
		return;
	else if(Ahel.size() == 1)
		Ahel[0] = std::polar<double>(1.0,phases[0].value);
	else if(Ahel.size() == 2)
	{
		double Aplus_mag = std::sqrt(1. - magsqs[0].value);
		if(std::isnan(Aplus_mag)) Aplus_mag = 0; // Handle this gracefully. Impose external constraints to stop this happening.
		Ahel[ 0] = std::polar<double>(std::sqrt(magsqs[0].value),phases[0].value);
		Ahel[+1] = std::polar<double>(Aplus_mag,phases[1].value);
	}
	else if(Ahel.size() == 3)
	{
		double Aminus_mag = std::sqrt(1. - magsqs[0].value - magsqs[1].value);
		if(std::isnan(Aminus_mag)) Aminus_mag = 0; // Handle this gracefully. Impose external constraints to stop this happening.
		Ahel[ 0] = std::polar<double>(std::sqrt(magsqs[0].value),phases[0].value);
		Ahel[+1] = std::polar<double>(std::sqrt(magsqs[1].value),phases[1].value);
		Ahel[-1] = std::polar<double>(Aminus_mag,phases[2].value);
	}
	else
	{
		std::cerr << "Bs2PhiKKSignalComponent can't handle this many helicities: " << Ahel.size() << std::endl;
		std::exit(-1);
	}
}
void Bs2PhiKKSignalComponent::UpdateBarriers()
{
	if(lineshape == "NR")
		return;
	double mphi = phimass.value;
	double Mres = KKpars[0].value;
	double m_min  = Bs2PhiKK::mK + Bs2PhiKK::mK;
	double m_max  = Bs2PhiKK::mBs - mphi;
	double m0_eff = m_min + (m_max - m_min) * (1 + std::tanh((Mres - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
	double pBs0 = DPHelpers::daughterMomentum(Bs2PhiKK::mBs, mphi, m0_eff);
	double pKK0 = DPHelpers::daughterMomentum(Mres, Bs2PhiKK::mK, Bs2PhiKK::mK);
	Bsbarrier.setparameters(BsBFradius.value,pBs0);
	KKbarrier.setparameters(KKBFradius.value,pKK0);
}
vector<ObservableRef> Bs2PhiKKSignalComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {{phimass},magsqs,phases,KKpars})
		for(const auto& par: set)
			parameters.push_back(par.name);
	// Add barrier factors if this is a resonant component
	if(lineshape != "NR")
		for(const auto& par: {BsBFradius,KKBFradius})
			parameters.push_back(par.name);
	// Add fraction if this isn't a spline shape
	if(lineshape != "SP")
		parameters.push_back(fraction.name);
	return parameters;
}

