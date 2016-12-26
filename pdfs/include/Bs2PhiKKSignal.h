#ifndef Bs2PhiKKSignal_H
#define Bs2PhiKKSignal_H

#include <memory>
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
		double EvaluateComponent( DataPoint*, ComponentRef* );
		std::vector<std::string> PDFComponents();
	protected:
		std::vector<Bs2PhiKKComponent> components;
		// K+K− mass and helicity angles
		double        mKK    , ctheta_1    , ctheta_2    , phi    ;
		ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
		// Bs width splitting
		PhysPar dGsGs;
		// phi(1020) mass
		PhysPar phimass;
		// Options
		std::vector<std::string> componentnames;
		bool acceptance_moments;
		bool acceptance_histogram;
		// Acceptance
		std::unique_ptr<LegendreMomentShape> acc_m;
		std::shared_ptr<NDHist_Adaptive> acc_h;
		double acceptance; // Evaluate during ReadDataPoint()
	private:
		// Calculation
		double TotalDecayRate();
		double ComponentDecayRate(Bs2PhiKKComponent&); // For plotting individual components
		double ComponentDecayRate(Bs2PhiKKComponent&, std::string); // Pass option "odd" or "even"
		double TimeIntegratedDecayRate(std::complex<double>,std::complex<double>);
		void ReadDataPoint(DataPoint*);
		double p1stp3();
		void CalculateAcceptance();
		// Stuff to do on creation
		void Initialise();
		void MakePrototypes();
		Bs2PhiKKComponent ParseComponent(PDFConfigurator*, std::string, std::string);
};
#endif

