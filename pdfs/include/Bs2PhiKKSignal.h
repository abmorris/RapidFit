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
      // Other model parameters
      double        SwaveFrac    , PwaveFrac    , DwaveFrac    ;
      ObservableRef SwaveFracName, PwaveFracName, DwaveFracName;
      // Options
      vector<string> componentlist;
    private:
      // Calculation
      TComplex TotalAmplitude(bool);
      Bs2PhiKKComponent* Swave;
      Bs2PhiKKComponent* Pwave;
      Bs2PhiKKComponent* Dwave;
      void SetComponentAmplitudes();
      double PhaseSpace(double);
      // Stuff to do on creation
      void Initialise();
      void MakePrototypes();
      // Options
      string compName;
};
#endif

