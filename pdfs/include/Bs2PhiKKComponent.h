#ifndef __BS2PHIKKCOMPONENT_H__
#define __BS2PHIKKCOMPONENT_H__
// Std
#include <complex>
#include <string>
#include <vector>
#include <memory>
// RapidFit
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPWignerFunctionGeneral.hh"

class Bs2PhiKKComponent
{
	public:
		Bs2PhiKKComponent(PDFConfigurator*, std::string, std::string, int, std::string); // config, phi name, resonance name, spin
		Bs2PhiKKComponent(const Bs2PhiKKComponent&);
		~Bs2PhiKKComponent();
		std::string GetName() const {return KKname;}
		bool SetPhysicsParameters(ParameterSet* pars); // Update all the parameters. Return whether any have been modified
		std::vector<ObservableRef> GetPhysicsParameters() const;
		std::complex<double> Amplitude(const double, const double, const double, const double) const; // KK_M, Phi_angle, cos_theta1, cos_theta2
		std::complex<double> Amplitude(const double, const double, const double, const double, const std::string) const; // Same but with an option "even" or "odd"
		static double mBs;
		static double mK;
		static double mpi;
		// Simplify the case where a value and a name correspond 1:1
		struct PhysPar
		{
			// Construct this however you want
			PhysPar() {}
			PhysPar(ObservableRef _name) : name(_name), value(0) {}
			PhysPar(ObservableRef _name, double _value) : name(_name), value(_value) {}
			PhysPar(PDFConfigurator* config, std::string _name) : name(config->getName(_name)), value(0) {}
			PhysPar(PDFConfigurator* config, std::string _name, double _value) : name(config->getName(_name)), value(_value) {}
			PhysPar(const PhysPar& other) : value(other.value), name(other.name) {}
			bool Update(const ParameterSet* pars);
			double value;
			ObservableRef name;
		};
	private:
		// Floatable parameters
		PhysPar fraction; // Unnormalised variable to control the relative contribution of each resonance. Do not use at the fit fraction!!
		std::vector<std::complex<double>> Ahel;  // Helicity amplitude(s)
		std::vector<int> helicities; // Store the possible values of helicity to enable looping over A(helicities[i])
		// Polarisation amplitude components (perp, zero, para)
		std::vector<PhysPar> magsqs; // Square of magnitudes: para will be calculated from the other two
		std::vector<PhysPar> phases; // Phases
		// Resonance parameters
		std::vector<PhysPar> KKpars; // Mass and width of Breit Wigner, or mass, g_pipi and R=(g_KK/g_pipi) of Flatte. Empty for non-resonant
		double mphi; // For orbital/barrier factor calculations. Assume the PDF already has this (for p1*p3 calculation)
		// Fixed parameters
		int Jphi; // Spin of the phi (P-wave, 1)
		int JKK; // Spin of the KK resonance (0, 1 or 2)
		double RBs; // Bs barrier factor radius
		double RKK; // KK barrier factor radius
		std::string lineshape; // Choose the resonance shape: "BW", "FT" or "NR"
		ObservableRef phiMassname; // The name decides which set of PhysicsParameters it will look for in the RapidFit XML
		std::string KKname;
		bool fixedamplitudes; // Save some time updating parameters
		bool fixedlineshape;
		bool fixedphimass;
		bool init; // Have we had our first SetPhysicsParameters?
		void Initialise();
		void UpdateAmplitudes();
		void UpdateLineshape();
		std::complex<double> A(const int) const; // Polarisation amplitude coefficients
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

