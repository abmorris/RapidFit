/** @class Bs2PhiKKSignal Bs2PhiKKSignal.h
 *
 *  RapidFit PDF for Bs2PhiKKSignal
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
#ifndef Bs2PhiKKSignal_H
#define Bs2PhiKKSignal_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKComponent.h"
#include "LegendreMomentShape.h"
#include "TKDTree.h"

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
      vector<string> PDFComponents();
    protected:
      // K+K− mass and helicity angles
      double        mKK    , ctheta_1    , ctheta_2    , phi    ;
      ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
      // Amplitude Parameters
      double        ASzerosq      , APperpsq      , APzerosq      , APparasq      , ADperpsq      , ADzerosq      , ADparasq      ;
      ObservableRef                 APperpsqName  , APzerosqName                  , ADperpsqName  , ADzerosqName                  ;
      double        deltaSzero    , deltaPperp    , deltaPzero    , deltaPpara    , deltaDperp    , deltaDzero    , deltaDpara    ;
      ObservableRef deltaSzeroName, deltaPperpName, deltaPzeroName, deltaPparaName, deltaDperpName, deltaDzeroName, deltaDparaName;
      // Control the size of the components
      double        NonResFrac    , SwaveFrac    , PwaveFrac    , DwaveFrac    ;
      ObservableRef NonResFracName, SwaveFracName, PwaveFracName, DwaveFracName;
      // Resonance parameters
      double        fzeroMass    , fzerogpipi    , fzeroRg    , phiMass    , phiWidth    , ftwoMass    , ftwoWidth    ;
      ObservableRef fzeroMassName, fzerogpipiName, fzeroRgName, phiMassName, phiWidthName, ftwoMassName, ftwoWidthName;
      // Barrier factor radius. Pull these from the options
      double RBs, RKK;
      // Options
      vector<string> componentlist;
      bool floatResPars;
      bool acceptance_moments;
      bool acceptance_histogram;
      // Acceptance
      LegendreMomentShape* acc_m;
      // Acceptance object
      TKDTreeID* accbinner;
      vector <double> accbincontent;
    private:
      // Calculation
      double TotalDecayRate();
      double ComponentDecayRate(Bs2PhiKKComponent*);
      double ComponentDecayRate(Bs2PhiKKComponent*, string);
      void ReadDataPoint(DataPoint*);
      Bs2PhiKKComponent* Swave;
      Bs2PhiKKComponent* Pwave;
      Bs2PhiKKComponent* Dwave;
      Bs2PhiKKComponent* NonRes;
      vector<Bs2PhiKKComponent*> components;
      void SetComponentAmplitudes();
      void SetResonanceParameters();
      double p1stp3(double);
      double Acceptance();
      // Stuff to do on creation
      void Initialise();
      void MakePrototypes();
};
#endif

