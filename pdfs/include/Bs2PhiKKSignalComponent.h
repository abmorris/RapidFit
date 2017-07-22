#ifndef Bs2PhiKKSignalComponent_H
#define Bs2PhiKKSignalComponent_H
// Std
#include <vector>
#include <memory>
#include <map>
// RapidFit
#include "Bs2PhiKK.h"
#include "DPBarrierFactor.hh"

class Bs2PhiKKSignalComponent
{
	public:
		Bs2PhiKKSignalComponent(PDFConfigurator* config, std::string phiname, std::string KKname, int _LBs, int _Lphi, int _LKK, std::string _lineshape, const std::vector<bool>& _UseObservable); 
		Bs2PhiKKSignalComponent(const Bs2PhiKKSignalComponent&);
		~Bs2PhiKKSignalComponent() {}
		void SetPhysicsParameters(ParameterSet* pars);
		std::vector<ObservableRef> GetPhysicsParameters() const;
		// These Amplitude functions return a 2-element array of the complex amplitudes of the B and Bbar decays
		Bs2PhiKK::amplitude_t Amplitude(const double&, const double&, const double&, const double&) const; // KK_M, Phi_angle, cos_theta1, cos_theta2
		Bs2PhiKK::amplitude_t Amplitude(const double&, const double&, const double&, const double&, const std::string) const; // Same but with an option "even" or "odd"
	private:
		const std::vector<bool> UseObservable; // Cache whether or not we use each dimension. Saves a lot of cycles compared to repeatedly calling ObservableNames.count()
		// Floatable parameters
		Bs2PhiKK::PhysPar fraction; // Unnormalised variable to control the relative contribution of each resonance. Do not use as the fit fraction!!
		Bs2PhiKK::PhysPar BsBFradius; // Bs barrier factor radius
		Bs2PhiKK::PhysPar KKBFradius; // KK barrier factor radius
		// Polarisation amplitude components (perp, zero, para)
		std::vector<Bs2PhiKK::PhysPar> mags; // Magnitudes
		std::vector<Bs2PhiKK::PhysPar> phases; // Phases
		std::map<int,std::complex<double>> Ahel; // Helicity amplitudes as complex numbers
		// Resonance parameters
		Bs2PhiKK::PhysPar phimass; // Mass of the "phi"
		std::vector<Bs2PhiKK::PhysPar> KKpars; // Mass and width of Breit Wigner, or mass, g_pipi and R=(g_KK/g_pipi) of Flatte. Empty for non-resonant
		int Lphi; // Spin of the "phi" (might actually be f0)
		int LKK; // Spin of the KK resonance (0, 1 or 2)
		int LBs; // Orbital angular momentum of the Bs (default min LB = |LKK-1| for LKK = 0, 1, 2)
		std::string lineshape; // Choose the resonance shape: "BW", "FT" or "NR"
		bool firstnonres; // If the "first" KK pair is nonresonant
		// Helper functions
		void Initialise();
		void UpdateAmplitudes();
		void UpdateBarriers();
		double F(const int&, const double&, const double&) const; // Angular distribution: helicity, datapoint
		std::complex<double> (Bs2PhiKKSignalComponent::*AngularPart)(const double&, const double&, const double&) const;
		std::complex<double> AngularPartNonRes(const double&, const double&, const double&) const;
		std::complex<double> AngularPartNoPhi(const double&, const double&, const double&) const;
		std::complex<double> AngularPartDefault(const double&, const double&, const double&) const;
		std::complex<double> (Bs2PhiKKSignalComponent::*MassPart)(const double&) const;
		std::complex<double> MassPartNonRes(const double&) const;
		std::complex<double> MassPartResonant(const double&) const;
		double OFBF(const double&) const; // Product of orbital and barrier factors
		// Wigner d-functions for the angular-dependent part
		double (*wignerKK)(double, int);
		double (*wignerPhi)(double, int);
		// Blatt-Weisskopf barrier penetration factors
		DPBarrierFactor Bsbarrier;
		DPBarrierFactor KKbarrier;
		// Resonance lineshape function for the mass-dependent part
		std::unique_ptr<DPMassShape> KKLineShape {};
};
#endif

