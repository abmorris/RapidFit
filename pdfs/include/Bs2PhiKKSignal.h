#ifndef Bs2PhiKKSignal_H
#define Bs2PhiKKSignal_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKSignalComponent.h"
#include "LegendreMomentShape.h"

class Bs2PhiKKSignal : public BasePDF, public Bs2PhiKK
{
	public:
		// *structors
		Bs2PhiKKSignal(PDFConfigurator*);
		Bs2PhiKKSignal(const Bs2PhiKKSignal&);
		~Bs2PhiKKSignal() {}
		// Required methods
		double Evaluate(DataPoint*);
		double Normalisation(PhaseSpaceBoundary*);
		bool SetPhysicsParameters(ParameterSet*);
		// Extra stuff
		double EvaluateComponent(DataPoint*, ComponentRef* );
		std::vector<std::string> PDFComponents();
	private:
		std::map<std::string,Bs2PhiKKSignalComponent> components; // Iterable list of amplitude components
		std::vector<std::string> componentnames; // List of names for plotting purposes only
		// Parameters used outside the amplitude calculation
		Bs2PhiKK::PhysPar GH, GL; // Bs mass eigenstate widths
		std::array<Bs2PhiKK::PhysPar,2> thraccscale; // threshold mass acceptance
		// Mass resolution variables
		Bs2PhiKK::PhysPar mKKres_sigmazero; // variable that parameterises the mass resolution
		std::map<std::string,double> mKKresconfig; // configuration parameters
		// Options
		bool acceptance_moments; // Use Legendre moments for acceptance
		// Acceptance
		std::array<std::unique_ptr<LegendreMomentShape>,2> acc_m;
		double phiwindowweight;
		// Calculation of the matrix element
		double TimeIntegratedMsq(const Bs2PhiKK::amplitude_t&) const; // Receive a complex amplitude and turn it into a |M|²
		typedef double (Bs2PhiKKSignal::*MsqFunc_t)(const double&, const double&, const double&, const double&, const std::string&) const;
		double TotalMsq(const double&, const double&, const double&, const double&, const std::string& dummy = "") const; // Calculate the total |M|². An MsqFunc_t object can point to this
		double ComponentMsq(const double&, const double&, const double&, const double&, const std::string&) const; // Calculate the |M|² of a single component. An MsqFunc_t object can point to this
		double InterferenceMsq(const double&, const double&, const double&, const double&, const std::string& dummy = "") const; // Calculate the difference between the total |M|² and the sum of individual |M|²s. An MsqFunc_t object can point to this
		// Turn the matrix element into the PDF
		double Evaluate_Base(const double&, const unsigned&, const double&, const double&, const double&, const double&) const;
		double p1stp3(const double&) const;
		double Convolve(MsqFunc_t, const double&, const double&, const double&, const double&, const std::string&) const; // Take one of the three above functions and convolve it with a double Gaussian for m(KK) resolution
		double Acceptance(const unsigned&, const double&, const double&, const double&, const double&) const;
		// Stuff to do on creation
		void MakePrototypes();
};
#endif

