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
		typedef double (Bs2PhiKKSignal::*MsqFunc_t)(const Bs2PhiKK::datapoint_t&, const std::string&) const;
		std::map<std::string,Bs2PhiKKSignalComponent> components; // Iterable list of amplitude components
		std::vector<std::string> componentnames; // List of names for plotting purposes only
		// Parameters used outside the amplitude calculation
		Bs2PhiKK::PhysPar dGsGs; // Bs width splitting
		Bs2PhiKK::PhysPar phimass; // phi(1020) mass
		Bs2PhiKK::PhysPar thraccscale; // threshold mass acceptance
		// Mass resolution variables
		Bs2PhiKK::PhysPar mKKres_sigmazero; // variable that parameterises the mass resolution
		std::map<std::string,double> mKKresconfig; // configuration parameters
		// Options
		bool acceptance_moments; // Use Legendre moments for acceptance
		// Acceptance objects
		std::array<std::unique_ptr<LegendreMomentShape>,2> acc_m;
		// Calculation of the matrix element
		double TimeIntegratedMsq(const Bs2PhiKK::amplitude_t&) const; // Receive a complex amplitude and turn it into a |M|²
		double TotalMsq(const Bs2PhiKK::datapoint_t&, const std::string& dummy = "") const; // Calculate the total |M|². An MsqFunc_t object can point to this
		double ComponentMsq(const Bs2PhiKK::datapoint_t&, const std::string&) const; // Calculate the |M|² of a single component. An MsqFunc_t object can point to this
		double InterferenceMsq(const Bs2PhiKK::datapoint_t&, const std::string& dummy = "") const; // Calculate the difference between the total |M|² and the sum of individual |M|²s. An MsqFunc_t object can point to this
		double Convolve(MsqFunc_t, const Bs2PhiKK::datapoint_t&, const std::string&) const; // Take one of the three above functions and convolve it with a double Gaussian for m(KK) resolution
		// Turn the matrix element into the PDF
		double Evaluate_Base(const double, const Bs2PhiKK::datapoint_t&) const;
		double p1stp3(const double&) const;
		double Acceptance(const Bs2PhiKK::datapoint_t&) const;
		// Stuff to do on creation
		void MakePrototypes();
		Bs2PhiKKSignalComponent ParseComponent(PDFConfigurator*, std::string, std::string) const;
};
#endif

