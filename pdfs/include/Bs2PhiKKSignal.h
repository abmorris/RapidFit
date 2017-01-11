#ifndef Bs2PhiKKSignal_H
#define Bs2PhiKKSignal_H

#include <memory>
#include <map>
#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKComponent.h"
#include "LegendreMomentShape.h"
#include "NDHist_Adaptive.h"

class Bs2PhiKKSignal : public BasePDF
{
	public:
		// *structors
		Bs2PhiKKSignal(PDFConfigurator*);
		Bs2PhiKKSignal(const Bs2PhiKKSignal&);
		~Bs2PhiKKSignal();
		// Required methods
		double Evaluate(DataPoint*);
		double Normalisation(PhaseSpaceBoundary*);
		bool SetPhysicsParameters(ParameterSet*);
		// Extra stuff
		double EvaluateComponent(DataPoint*, ComponentRef* );
		std::vector<std::string> PDFComponents();
	private:
		std::vector<Bs2PhiKKComponent> components;
		int index; // Current datapoint index for caching
		// K+K− mass and helicity angles
		ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
		// Bs width splitting
		Bs2PhiKKComponent::PhysPar dGsGs;
		// phi(1020) mass
		Bs2PhiKKComponent::PhysPar phimass;
		// Options
		std::vector<std::string> componentnames; // List of components to be drawn
		bool acceptance_moments; // Use Legendre moments for acceptance
		bool acceptance_histogram; // Use adaptively-binned histogram for acceptance
		// Acceptance objects
		std::unique_ptr<LegendreMomentShape> acc_m;
		std::shared_ptr<NDHist_Adaptive> acc_h;
		// Calculation
		double TotalMsq(const std::array<double,4>&) const;
		double ComponentMsq(const Bs2PhiKKComponent&, const std::array<double,4>&) const; // For plotting individual components
		double ComponentMsq(const Bs2PhiKKComponent&, const std::array<double,4>&, const std::string) const; // Pass option "odd" or "even"
		double TimeIntegratedMsq(const std::array<std::complex<double>,2>&) const;
		std::array<double,4> ReadDataPoint(DataPoint*) const;
		double p1stp3(const double&) const;
		double Acceptance(const std::array<double,4>&) const;
		// Stuff to do on creation
		void Initialise();
		void MakePrototypes();
		Bs2PhiKKComponent ParseComponent(PDFConfigurator*, std::string, std::string) const;
};
#endif

