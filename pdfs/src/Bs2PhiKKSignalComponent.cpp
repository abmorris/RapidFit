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
	, Bsbarrier(DPBarrierFactor(abs(JKK-1), 1.0, 0)) // Min L_B approximation
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
	if(lineshape != "NR" && lineshape != "SP")
		for(int lambda = -lambda_max; lambda <= lambda_max; lambda++)
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
		case 3:
			magsqsuffices = {"Aperpsq","Azerosq"};
			phasesuffices = {"deltaperp","deltazero","deltapara"};
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
		Ahel[lambda] = std::polar<double>(sqrt(1. / (double)n), 0);
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
		KKLineShape = std::unique_ptr<DPBWResonanceShape>(new DPBWResonanceShape(KKpars[0].value, KKpars[1].value, JKK, Bs2PhiKK::mK, Bs2PhiKK::mK, KKpars[2].value));
	// Flatte
	else if(lineshape=="FT")
		KKLineShape = std::unique_ptr<DPFlatteShape>(new DPFlatteShape(KKpars[0].value, KKpars[1].value, Bs2PhiKK::mpi, Bs2PhiKK::mpi, KKpars[1].value*KKpars[2].value, Bs2PhiKK::mK, Bs2PhiKK::mK));
	else if(lineshape=="SP")
	{
		std::vector<double> breakpoints;
		for(unsigned i = 0; i < KKpars.size(); i+=3)
			breakpoints.push_back(KKpars[i].value);
		KKLineShape = std::unique_ptr<DPComplexSpline>(new DPComplexSpline(breakpoints));
	}
	else
		KKLineShape = std::unique_ptr<DPNonresonant>(new DPNonresonant());
	// Build the barrier factor and Wigner function objects
	wignerPhi = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
	switch (JKK) // I hate this but I'd rather it just worked...
	{
		case 0:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ0>(new DPWignerFunctionJ0());
			break;
		case 1:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ1>(new DPWignerFunctionJ1());
			break;
		case 2:
			wignerKK  = std::unique_ptr<DPWignerFunctionJ2>(new DPWignerFunctionJ2());
			break;
		default:
			wignerKK  = std::unique_ptr<DPWignerFunctionGeneral>(new DPWignerFunctionGeneral(JKK)); // This should only happen for the rho_3 (1690)
			break;
	}
	UpdateAmplitudes();
	Bs2PhiKK::UpdateLineshape(lineshape, *KKLineShape, KKpars);
}
/*****************************************************************************/
// Angular part of the amplitude
std::complex<double> Bs2PhiKKSignalComponent::F(const int lambda, const double Phi, const double ctheta_1, const double ctheta_2) const
{
	const double d_phi = wignerPhi->function(ctheta_1, lambda, 0);
	const double d_KK = wignerKK->function(ctheta_2, lambda, 0);
	return d_phi * d_KK * std::polar<double>(1, lambda*Phi);
}
std::complex<double> Bs2PhiKKSignalComponent::AngularPart(const double phi, const double ctheta_1, const double ctheta_2) const
{
	std::complex<double> angularPart(0, 0);
	if(lineshape == "SP") // Phase already accounted-for by complex spline, so just return generic S-wave shape
		angularPart = F(0, phi, ctheta_1, ctheta_2);
	else if(Ahel.empty())
		return std::complex<double>(1, 0); // Must be non-resonant
	else
		for(const auto& A : Ahel)
		{
			int lambda = A.first;
			angularPart += A.second * F(lambda, phi, ctheta_1, ctheta_2);
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
	double orbitalFactor = std::pow(pBs/Bs2PhiKK::mBs, abs(JKK-1))* // Min L_B approximation
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
	double mKK = datapoint[0];
	double phi = datapoint[1];
	double ctheta_1 = datapoint[2];
	double ctheta_2 = datapoint[3];
	Bs2PhiKK::amplitude_t angularPart = {AngularPart(phi, ctheta_1, ctheta_2), AngularPart(-phi, -ctheta_1, -ctheta_2)};
	std::complex<double> massPart = MassPart(mKK);
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
// The full amplitude with an option.
Bs2PhiKK::amplitude_t Bs2PhiKKSignalComponent::Amplitude(const Bs2PhiKK::datapoint_t& datapoint, const std::string option) const
{
	if(Ahel.empty() || option == "" || option.find("odd") == std::string::npos || option.find("even") == std::string::npos ) return Amplitude(datapoint);
	double mKK = datapoint[0];
	double phi = datapoint[1];
	double ctheta_1 = datapoint[2];
	double ctheta_2 = datapoint[3];
	// Angular part
	Bs2PhiKK::amplitude_t angularPart = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
	std::complex<double> Apara = std::polar(sqrt(1. - magsqs[0].value - magsqs[1].value),phases[2].value);
	// Temporary helicity amplitudes
	std::vector<std::complex<double>> HelAmp;
	// CP-odd component
	if(option.find("odd") != std::string::npos)
		HelAmp = {-Aperp/sqrt(2.), std::complex<double>(0, 0), Aperp/sqrt(2.)};
	// CP-even component
	else
		HelAmp = {Apara/sqrt(2.), Ahel.at(0), Apara/sqrt(2.)};
	for(const auto& A : Ahel)
	{
		int lambda = A.first;
		angularPart[false] += HelAmp[lambda+1] * F(lambda, phi, ctheta_1, ctheta_2);
		angularPart[true] += HelAmp[lambda+1] * F(lambda, -phi, -ctheta_1, -ctheta_2);
	}
	std::complex<double> massPart = MassPart(mKK);
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
	else if(Ahel.size() == 3)
	{
		std::complex<double> Aperp = std::polar(sqrt(magsqs[0].value),phases[0].value);
		double Apara_mag = sqrt(1. - magsqs[0].value - magsqs[1].value);
		if(std::isnan(Apara_mag)) Apara_mag = 0; // Handle this gracefully. Impose external constraints to stop this happening.
		std::complex<double> Apara = std::polar(Apara_mag,phases[2].value);
		Ahel[-1] = (Apara - Aperp)/sqrt(2.); // A− = (A‖ − A⊥)/sqrt(2)
		Ahel[ 0] = std::polar(sqrt(magsqs[1].value),phases[1].value);;
		Ahel[+1] = (Apara + Aperp)/sqrt(2.); // A+ = (A‖ + A⊥)/sqrt(2)
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

