#ifndef Bs2PhiKKBackground_H
#define Bs2PhiKKBackground_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKHelpers.h"
#include "Bs2PhiKKBackgroundComponent.h"
#include "LegendreMomentShape.h"

class Bs2PhiKKBackground : public BasePDF
{
	public:
		// *structors
		Bs2PhiKKBackground(PDFConfigurator*);
		Bs2PhiKKBackground(const Bs2PhiKKBackground&);
		~Bs2PhiKKBackground() {}
		// Required methods
		double Evaluate(DataPoint*);
		double Normalisation(PhaseSpaceBoundary*);
		bool SetPhysicsParameters(ParameterSet*);
		// Extra stuff
		double EvaluateComponent(DataPoint*, ComponentRef* );
		std::vector<std::string> PDFComponents();
	private:
		std::map<std::string,Bs2PhiKKBackgroundComponent> components; // Iterable list of amplitude components
		ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName; // Datapoint stuff: K+K− mass and helicity angles
		// Retrieve an array of doubles from a RapidFit Datapoint object
		Bs2PhiKK::datapoint_t ReadDataPoint(DataPoint*) const;
		// Stuff to do on creation
		void Initialise();
		void MakePrototypes();
		Bs2PhiKKBackgroundComponent ParseComponent(PDFConfigurator*, std::string) const;
};
#endif

