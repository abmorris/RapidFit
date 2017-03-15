#ifndef Bs2PhiKKBackground_H
#define Bs2PhiKKBackground_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKBackgroundComponent.h"
#include "LegendreMomentShape.h"

class Bs2PhiKKBackground : public BasePDF, public Bs2PhiKK
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
		// Stuff to do on creation
		void MakePrototypes();
		Bs2PhiKKBackgroundComponent ParseComponent(PDFConfigurator*, std::string) const;
};
#endif

