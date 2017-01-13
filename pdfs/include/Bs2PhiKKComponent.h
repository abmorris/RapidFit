#ifndef __BS2PHIKKCOMPONENT_H__
#define __BS2PHIKKCOMPONENT_H__
// Std
#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <memory>
#include <array>
#include <map>
// RapidFit
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPWignerFunctionGeneral.hh"

class Bs2PhiKKComponent
{
	public:
		Bs2PhiKKComponent() {}
		Bs2PhiKKComponent(PDFConfigurator*, std::string, std::string, int, std::string); // config, phi name, resonance name, spin
		Bs2PhiKKComponent(const Bs2PhiKKComponent&);
		Bs2PhiKKComponent(Bs2PhiKKComponent&&);
		~Bs2PhiKKComponent();
		Bs2PhiKKComponent& operator=(const Bs2PhiKKComponent&);
		Bs2PhiKKComponent& operator=(Bs2PhiKKComponent&&);
		void SetPhysicsParameters(ParameterSet* pars);
		std::vector<ObservableRef> GetPhysicsParameters() const;
		// These Amplitude functions return a 2-element array of the complex amplitudes of the B and Bbar decays
		std::array<std::complex<double>,2> Amplitude(const std::array<double,4>&) const; // {KK_M, Phi_angle, cos_theta1, cos_theta2}
		std::array<std::complex<double>,2> Amplitude(const std::array<double,4>&, const std::string) const; // Same but with an option "even" or "odd"
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
			PhysPar(PhysPar&& other): value(std::move(other.value)), name(std::move(other.name)) {}
			PhysPar& operator=(const Bs2PhiKKComponent::PhysPar& other) {name = other.name;return *this;}
			PhysPar& operator=(Bs2PhiKKComponent::PhysPar&& other) {name = std::move(other.name);return *this;}
			void Update(const ParameterSet* pars) {value = pars->GetPhysicsParameter(name)->GetValue(); if(std::isnan(value)) std::cerr << name.Name() << " has been given a nan value!" << std::endl; }
			double value;
			ObservableRef name;
		};
	private:
		// Floatable parameters
		PhysPar fraction; // Unnormalised variable to control the relative contribution of each resonance. Do not use as the fit fraction!!
		PhysPar phimass; // For use when calculating the orbital factor. Assume the PDF has already declared the need for this
		std::map<int,std::complex<double>> Ahel;
		// Polarisation amplitude components (perp, zero, para)
		std::vector<PhysPar> magsqs; // Square of magnitudes: para will be calculated from the other two
		std::vector<PhysPar> phases; // Phases
		// Resonance parameters
		std::vector<PhysPar> KKpars; // Mass and width of Breit Wigner, or mass, g_pipi and R=(g_KK/g_pipi) of Flatte. Empty for non-resonant
		// Fixed parameters
		int Jphi; // Spin of the phi (P-wave, 1)
		int JKK; // Spin of the KK resonance (0, 1 or 2)
		double RBs; // Bs barrier factor radius
		double RKK; // KK barrier factor radius
		std::string lineshape; // Choose the resonance shape: "BW", "FT" or "NR"
		void Initialise();
		void UpdateAmplitudes();
		void UpdateLineshape();
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

