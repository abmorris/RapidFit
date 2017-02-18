#ifndef Bs2PhiKKSignalComponent_H
#define Bs2PhiKKSignalComponent_H
// Std
#include <vector>
#include <memory>
#include <map>
// RapidFit
#include "Bs2PhiKK.h"
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPWignerFunctionGeneral.hh"

class Bs2PhiKKSignalComponent
{
	public:
		Bs2PhiKKSignalComponent() {}
		Bs2PhiKKSignalComponent(PDFConfigurator*, std::string, std::string, int, std::string); // config, phi name, resonance name, spin
		Bs2PhiKKSignalComponent(const Bs2PhiKKSignalComponent&);
		~Bs2PhiKKSignalComponent() {}
		Bs2PhiKKSignalComponent& operator=(const Bs2PhiKKSignalComponent&);
		void SetPhysicsParameters(ParameterSet* pars);
		std::vector<ObservableRef> GetPhysicsParameters() const;
		// These Amplitude functions return a 2-element array of the complex amplitudes of the B and Bbar decays
		Bs2PhiKK::amplitude_t Amplitude(const Bs2PhiKK::datapoint_t&) const; // {KK_M, Phi_angle, cos_theta1, cos_theta2}
		Bs2PhiKK::amplitude_t Amplitude(const Bs2PhiKK::datapoint_t&, const std::string) const; // Same but with an option "even" or "odd"
	private:
		// Floatable parameters
		Bs2PhiKK::PhysPar fraction; // Unnormalised variable to control the relative contribution of each resonance. Do not use as the fit fraction!!
		Bs2PhiKK::PhysPar phimass; // For use when calculating the orbital factor. Assume the PDF has already declared the need for this
		Bs2PhiKK::PhysPar BsBFradius; // Barrier factor radii
		Bs2PhiKK::PhysPar KKBFradius;
		std::map<int,std::complex<double>> Ahel;
		// Polarisation amplitude components (perp, zero, para)
		std::vector<Bs2PhiKK::PhysPar> magsqs; // Square of magnitudes: para will be calculated from the other two
		std::vector<Bs2PhiKK::PhysPar> phases; // Phases
		// Resonance parameters
		std::vector<Bs2PhiKK::PhysPar> KKpars; // Mass and width of Breit Wigner, or mass, g_pipi and R=(g_KK/g_pipi) of Flatte. Empty for non-resonant
		// Fixed parameters
		int Jphi; // Spin of the phi (P-wave, 1)
		int JKK; // Spin of the KK resonance (0, 1 or 2)
		std::string lineshape; // Choose the resonance shape: "BW", "FT" or "NR"
		void Initialise();
		void UpdateAmplitudes();
		void UpdateLineshape();
		void UpdateBarriers();
		std::complex<double> F(const int, const double, const double, const double) const; // Angular distribution: helicity, phi, costheta1, costheta2
		std::complex<double> AngularPart(const double, const double, const double) const; // index, phi, costheta1, costheta2
		double OFBF(const double) const; // Product of orbital and barrier factors
		// Wigner d-functions for the angular-dependent part
		std::unique_ptr<DPWignerFunction> wignerKK {};
		std::unique_ptr<DPWignerFunction> wignerPhi {};
		// Blatt-Weisskopf barrier penetration factors
		DPBarrierFactor Bsbarrier;
		DPBarrierFactor KKbarrier;
		// Resonance lineshape function for the mass-dependent part
		std::unique_ptr<DPMassShape> KKLineShape {};
};
#endif

