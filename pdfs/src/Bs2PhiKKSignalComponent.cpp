// Std
#include <algorithm>
// Self
#include "Bs2PhiKKSignalComponent.h"
// RapidFit Dalitz Plot Libraries
#include "DPBWResonanceShape.hh"
#include "DPFlatteShape.hh"
#include "DPNonresonant.hh"
#include "DPHelpers.hh"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ2.hh"
// Constructor
Bs2PhiKKSignalComponent::Bs2PhiKKSignalComponent(PDFConfigurator* config, std::string phiname, std::string KKname, int _LBs, int _Lphi, int _LKK, std::string _philineshapename, std::string _KKlineshapename, const std::vector<bool>& _UseObservable)
	: LBs(_LBs)
	, Lphi(_Lphi)
	, LKK(_LKK)
	, philineshapename(_philineshapename)
	, KKlineshapename(_KKlineshapename)
	, Bsbarrier(DPBarrierFactor(LBs, 1.0, 0))
	, KKbarrier(DPBarrierFactor(LKK, 3.0, 0))
	, UseObservable(_UseObservable)
	, AngularPart(&Bs2PhiKKSignalComponent::AngularPartDefault)
{
	// sort phiname and KKname alphabetically so that Bs->AB has the same fraction and amplitues as Bs->BA
	std::vector<std::string> daughternames = {phiname, KKname};
	std::sort(daughternames.begin(), daughternames.end());
	std::string decayname = daughternames[0] + "_" + daughternames[1];
	fraction = Bs2PhiKK::PhysPar(config,decayname+"_fraction");
	// Barrier factors
	BsBFradius = Bs2PhiKK::PhysPar(config,"BsBFradius");
	KKBFradius = Bs2PhiKK::PhysPar(config,"KKBFradius");
	// phi resonance parameters
	for(string phipar: Bs2PhiKK::LineShapeParameterNames(phiname,philineshapename))
	{
		phipars.push_back(Bs2PhiKK::PhysPar(config,phipar));
	}
	// KK resonance parameters
	for(string KKpar: Bs2PhiKK::LineShapeParameterNames(KKname,KKlineshapename))
	{
		KKpars.push_back(Bs2PhiKK::PhysPar(config,KKpar));
	}
	if(phipars.empty() && philineshapename!="NR")
	{
		philineshapename = "NR";
		std::cerr << "Bs2PhiKKSignalComponent WARNING: unknown lineshape '" << philineshapename << "'. Treating this component non-resonant." << std::endl;
	}
	if(KKpars.empty() && KKlineshapename!="NR")
	{
		KKlineshapename = "NR";
		std::cerr << "Bs2PhiKKSignalComponent WARNING: unknown lineshape '" << KKlineshapename << "'. Treating this component non-resonant." << std::endl;
	}
	// Helicity amplitude observables
	int lambda_max = std::min(Lphi, LKK); // Maximum helicity
	std::vector<int> helicities;
	bool usephi = false, usect1 = false, usect2 = false;
	for(const auto& name: config->GetPhaseSpaceBoundary()->GetAllNames())
	{
		if(name == "phi")
		{
			usephi = true;
		}
		else if(name == "ctheta_1")
		{
			usect1 = true;
		}
		else if(name == "ctheta_2")
		{
			usect2 = true;
		}
	}
	for(int lambda = usephi ? -lambda_max : 0; lambda <= lambda_max; lambda++)
	{
		helicities.push_back(lambda);
	}
	if(usect1 && usect2)
	{
		long unsigned int n = helicities.size();
		std::vector<std::string> phasesuffices;
		std::vector<std::string> magsuffices;
		switch(n)
		{
			case 0:
				break;
			case 1:
				phasesuffices = {"deltazero"};
				break;
			case 2:
				magsuffices = {"Azero"};
				phasesuffices = {"deltazero","deltaplus"};
				break;
			case 3:
				magsuffices = {"Azero","Aplus"};
				phasesuffices = {"deltazero","deltaplus","deltaminus"};
				break;
			default:
				std::cerr << "Bs2PhiKKSignalComponent can't handle this many helicities: " << n << std::endl;
				std::exit(-1);
				break;
		}
		for(std::string suffix: magsuffices)
		{
			mags.push_back(Bs2PhiKK::PhysPar(config,decayname+"_"+suffix));
		}
		for(std::string suffix: phasesuffices)
		{
			phases.push_back(Bs2PhiKK::PhysPar(config,decayname+"_"+suffix));
		}
		// Make the helicity amplitude vector
		for(const auto& lambda: helicities)
		{
			Ahel[lambda] = std::polar<double>(std::sqrt(1. / (double)n), 0);
		}
	}
	Initialise();
}
// Copy constructor
Bs2PhiKKSignalComponent::Bs2PhiKKSignalComponent(const Bs2PhiKKSignalComponent& other)
	// Floatable parameters
	: fraction(other.fraction)
	, Ahel(other.Ahel)
	, mags(other.mags)
	, phases(other.phases)
	, phipars(other.phipars)
	, KKpars(other.KKpars)
	, BsBFradius(other.BsBFradius)
	, KKBFradius(other.KKBFradius)
	// Fixed parameters
	, LBs(other.LBs)
	, Lphi(other.Lphi)
	, LKK(other.LKK)
	// DP objects
	, Bsbarrier(other.Bsbarrier)
	, KKbarrier(other.KKbarrier)
	// Options
	, philineshapename(other.philineshapename)
	, KKlineshapename(other.KKlineshapename)
	, UseObservable(other.UseObservable)
	, AngularPart(&Bs2PhiKKSignalComponent::AngularPartDefault)
{
	Initialise();
}
/*****************************************************************************/
void Bs2PhiKKSignalComponent::Initialise()
{
	if(KKlineshapename=="BW")
	{
		KKLineShape = std::unique_ptr<DPBWResonanceShape>(new DPBWResonanceShape(KKpars[0].value, KKpars[1].value, LKK, Bs2PhiKK::mK, Bs2PhiKK::mK, KKpars[2].value));
	}
	// Flatte
	else if(KKlineshapename=="FT")
	{
		KKLineShape = std::unique_ptr<DPFlatteShape>(new DPFlatteShape(KKpars[0].value, KKpars[1].value, Bs2PhiKK::mpi, Bs2PhiKK::mpi, KKpars[1].value*KKpars[2].value, Bs2PhiKK::mK, Bs2PhiKK::mK));
	}
	else
	{
		KKLineShape = std::unique_ptr<DPNonresonant>(new DPNonresonant());
	}
	auto phistuff = std::make_tuple(&Lphi, &wignerPhi);
	auto KKstuff = std::make_tuple(&LKK, &wignerKK);
	for(auto& stuff: {phistuff, KKstuff})
	{
		const auto L = *std::get<0>(stuff);
		auto& wigner = *std::get<1>(stuff);
		switch (L)
		{
			case 0:
				wigner = &DPWignerFunctionJ0::function;
				break;
			case 1:
				wigner = &DPWignerFunctionJ1::function;
				break;
			case 2:
				wigner = &DPWignerFunctionJ2::function;
				break;
			default:
				std::cerr << "Can't construct Wigner function for spin " << L << std::endl;
				break;
		}
	}
	if(!UseObservable[Bs2PhiKK::_ctheta_1_] || !UseObservable[Bs2PhiKK::_ctheta_2_] || KKlineshapename == "NR" || philineshapename == "NR")
	{
		AngularPart = &Bs2PhiKKSignalComponent::AngularPartNonRes;
	}
	else if(!UseObservable[Bs2PhiKK::_phi_])
	{
		AngularPart = &Bs2PhiKKSignalComponent::AngularPartNoPhi;
	}
	UpdateAmplitudes();
	Bs2PhiKK::UpdateLineshape(KKlineshapename, *KKLineShape, KKpars);
}
/*****************************************************************************/
// Angular part of the amplitude
double Bs2PhiKKSignalComponent::F(const int& lambda, const double ctheta_1, const double ctheta_2) const
{
	const double d_phi = wignerPhi(ctheta_1, lambda);
	const double d_KK  = wignerKK (ctheta_2, lambda);
	return d_phi * d_KK;
}
std::complex<double> Bs2PhiKKSignalComponent::AngularPartNonRes(const double ctheta_1, const double ctheta_2, const double phi) const
{
	(void)ctheta_1;
	(void)ctheta_2;
	(void)phi;
	return Ahel.at(0); // Either nonresonant or helicity angles not present
}
std::complex<double> Bs2PhiKKSignalComponent::AngularPartNoPhi(const double ctheta_1, const double ctheta_2, const double phi) const
{
	(void)phi;
	std::complex<double> angularPart(0, 0);
	for(const auto& A : Ahel)
	{
		const int lambda = A.first;
		std::complex<double> thisterm = A.second * F(lambda, ctheta_1, ctheta_2);
		if(lambda != 0)
		{
			thisterm *= M_PI/3.;
		}
		angularPart += thisterm;
	}
	return angularPart;
}
std::complex<double> Bs2PhiKKSignalComponent::AngularPartDefault(const double ctheta_1, const double ctheta_2, const double phi) const
{
	std::complex<double> angularPart(0, 0);
	for(const auto& A : Ahel)
	{
		const int lambda = A.first;
		angularPart += A.second * F(lambda, ctheta_1, ctheta_2) * std::polar<double>(1.0, lambda*phi);
	}
	return angularPart;
}
// Orbital and barrier factor
double Bs2PhiKKSignalComponent::OFBF(const double mKK, const double mKK0, const double mphi0, const std::string& lineshape, const int _LKK, const DPBarrierFactor& barrier) const
{
	if(mKK < 2*Bs2PhiKK::mK)
	{
		return 0;
	}
	if(lineshape == "NR")
	{
		return 1;
	}
	// Momenta
	const double pBs = DPHelpers::daughterMomentum(Bs2PhiKK::mBs, mphi0, mKK);
	const double pKK = DPHelpers::daughterMomentum(mKK, Bs2PhiKK::mK, Bs2PhiKK::mK);
	// Masses
	const double mBs0 = Bs2PhiKK::mBs;
	// Orbital factor
	const double orbitalFactor = std::pow(pBs/mBs0, LBs) * std::pow(pKK/mKK0, _LKK);
	// Barrier factors
	const double barrierFactor = Bsbarrier.barrier(pBs) * barrier.barrier(pKK);
	return orbitalFactor*barrierFactor*std::sqrt(pBs*pKK);
}
// Mass-dependent part
std::complex<double> Bs2PhiKKSignalComponent::MassPart(const double mKK) const
{
	if(!UseObservable[Bs2PhiKK::_mKK_])
	{
		return std::complex<double>(1,0);
	}
	const double mphi0 = phipars[0].value; // pole mass of the "phi"
	const double mKK0 = KKpars[0].value; // pole mass of the "KK"
//	double norm =  std::abs(KKmasspartintegral);
//	std::cout << KKlineshapename << " norm: " << norm << "\n";
	std::complex<double> KKmasspart = KKLineShape->massShape(mKK) * OFBF(mKK, mKK0, mphi0, KKlineshapename, LKK, KKbarrier);
	return fraction.value * KKmasspart;
}
// The full amplitude.
Bs2PhiKK::amplitude_t Bs2PhiKKSignalComponent::Amplitude(const double mKK, const double ctheta_1, const double ctheta_2, const double phi) const
{
	std::complex<double> angularPartBs    = (this->*AngularPart)( ctheta_1,  ctheta_2,  phi);
	std::complex<double> angularPartBsbar = (this->*AngularPart)(-ctheta_1, -ctheta_2, -phi);
	std::complex<double> massPart = MassPart(mKK);
	return {massPart*angularPartBs, massPart*angularPartBsbar};
}
// The full amplitude with an option.
Bs2PhiKK::amplitude_t Bs2PhiKKSignalComponent::Amplitude(const double mKK, const double ctheta_1, const double ctheta_2, const double phi, const std::string option) const
{
	if(Ahel.empty() || option == "" || option.find("odd") == std::string::npos || option.find("even") == std::string::npos || !UseObservable[Bs2PhiKK::_phi_])
	{
		return Amplitude(mKK, ctheta_1, ctheta_2, phi);
	}
	// Angular part
	Bs2PhiKK::amplitude_t angularPart = {std::complex<double>(0, 0), std::complex<double>(0, 0)};
	std::complex<double> Aperp = std::polar<double>(mags[0].value,phases[0].value);
	std::complex<double> Apara = std::polar<double>(std::sqrt(1. - std::pow(mags[0].value,2) - std::pow(mags[1].value,2)),phases[2].value);
	// Temporary helicity amplitudes
	std::vector<std::complex<double>> HelAmp;
	// CP-odd component
	if(option.find("odd") != std::string::npos)
	{
		HelAmp = {-Aperp/std::sqrt(2.), std::complex<double>(0, 0), Aperp/std::sqrt(2.)};
	}
	// CP-even component
	else
	{
		HelAmp = {Apara/std::sqrt(2.), Ahel.at(0), Apara/std::sqrt(2.)};
	}
	for(const auto& A : Ahel)
	{
		int lambda = A.first;
		angularPart[false] += HelAmp[lambda+1] * F(lambda, ctheta_1, ctheta_2) * std::polar<double>(1, lambda*phi);
		angularPart[true] += HelAmp[lambda+1] * F(lambda, -ctheta_1, -ctheta_2) * std::polar<double>(1, -lambda*phi);
	}
	std::complex<double> massPart(1, 0);
	massPart = MassPart(mKK);
	return {massPart*angularPart[false], massPart*angularPart[true]};
}
/*****************************************************************************/
// Update everything from the parameter set
void Bs2PhiKKSignalComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
	// Update the parameters objects first
	for(auto* par: {&fraction, &KKBFradius, &BsBFradius})
	{
		par->Update(fitpars);
	}
	for(auto* set: {&mags,&phases,&phipars,&KKpars})
	{
		for(auto& par: *set)
		{
			par.Update(fitpars);
		}
	}
	UpdateBarriers();
	UpdateAmplitudes();
	Bs2PhiKK::UpdateLineshape(KKlineshapename, *KKLineShape, KKpars);
}
void Bs2PhiKKSignalComponent::UpdateAmplitudes()
{
	// Update the helicity amplitudes
	if(Ahel.empty())
	{
		return;
	}
	else if(Ahel.size() == 1)
	{
		Ahel[0] = std::polar<double>(1.0,phases[0].value);
	}
	else if(Ahel.size() == 2)
	{
		double Aplus_mag = std::sqrt(1. - std::pow(mags[0].value,2));
		if(std::isnan(Aplus_mag))
		{
			Aplus_mag = 0; // Handle this gracefully. Impose external constraints to stop this happening.
		}
		Ahel[ 0] = std::polar<double>(mags[0].value,phases[0].value);
		Ahel[+1] = std::polar<double>(Aplus_mag    ,phases[1].value);
	}
	else if(Ahel.size() == 3)
	{
		double Aminus_mag = std::sqrt(1. - std::pow(mags[0].value,2) - std::pow(mags[1].value,2));
		if(std::isnan(Aminus_mag))
		{
			Aminus_mag = 0; // Handle this gracefully. Impose external constraints to stop this happening.
		}
		Ahel[ 0] = std::polar<double>(mags[0].value,phases[0].value);
		Ahel[+1] = std::polar<double>(mags[1].value,phases[1].value);
		Ahel[-1] = std::polar<double>(Aminus_mag   ,phases[2].value);
	}
	else
	{
		std::cerr << "Bs2PhiKKSignalComponent can't handle this many helicities: " << Ahel.size() << std::endl;
		std::exit(-1);
	}
}
void Bs2PhiKKSignalComponent::UpdateBarriers()
{
	double mphi0 = phipars[0].value;
	double mKK0 = KKpars[0].value;
	double m_min  = Bs2PhiKK::mK + Bs2PhiKK::mK;
	double m_max  = Bs2PhiKK::mBs - mphi0;
	double m0_eff = m_min + (m_max - m_min) * (1 + std::tanh((mKK0 - (m_min + m_max) / 2) / (m_max - m_min))) / 2;
	double pBs0 = DPHelpers::daughterMomentum(Bs2PhiKK::mBs, mphi0, m0_eff);
	double pKK0 = DPHelpers::daughterMomentum(mKK0, Bs2PhiKK::mK, Bs2PhiKK::mK);
	Bsbarrier.setparameters(BsBFradius.value,pBs0);
	KKbarrier.setparameters(KKBFradius.value,pKK0);
}
vector<ObservableRef> Bs2PhiKKSignalComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {mags,phases,phipars,KKpars})
	{
		for(const auto& par: set)
		{
			parameters.push_back(par.name);
		}
	}
	for(const auto& par: {BsBFradius,KKBFradius,fraction})
	{
		parameters.push_back(par.name);
	}
	return parameters;
}

