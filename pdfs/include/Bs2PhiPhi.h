#ifndef Bs2PhiPhi_H
#define Bs2PhiPhi_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2PhiPhi : public BasePDF
{
	public:
		// *structors
		Bs2PhiPhi(PDFConfigurator*);
		Bs2PhiPhi(const Bs2PhiPhi&);
		~Bs2PhiPhi();
		// Required methods
		double Evaluate(DataPoint*);
		double Normalisation(PhaseSpaceBoundary*);
		bool SetPhysicsParameters(ParameterSet*);
		// Extra stuff
		double EvaluateComponent( DataPoint*, ComponentRef* );
		vector<string> PDFComponents();
	protected:
		// K+K− mass and helicity angles
		double        ctheta_1,     ctheta_2,     phi;
		ObservableRef ctheta_1Name, ctheta_2Name, phiName;
		// Bs width splitting
		double        dGsGs    ;
		ObservableRef dGsGsName;
		// Amplitude parameters
		double        Aperpsq      , Azerosq      , Aparasq      ;
		ObservableRef AperpsqName  , AzerosqName                 ;
		double        deltaperp    , deltazero    , deltapara    ;
		ObservableRef deltaperpName, deltazeroName, deltaparaName;
		// Options
		vector<string> componentlist;
		bool numerical;
	private:
		// Stuff to do on creation
		void Initialise();
		void MakePrototypes();
		// Options
		string compName;
};
#endif

