/** @class Bs2PhiKKPwave Bs2PhiKKPwave.h
 *
 *  RapidFit PDF for Bs2PhiKKPwave
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
#ifndef Bs2PhiKKPwave_H
#define Bs2PhiKKPwave_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "Bs2PhiKKComponent.h"

class Bs2PhiKKPwave : public BasePDF
{
    public:
      // *structors
      Bs2PhiKKPwave(PDFConfigurator*);
      Bs2PhiKKPwave(const Bs2PhiKKPwave&);
      ~Bs2PhiKKPwave();
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
      // Parameters
      double        Aperpsq      , Azerosq      , Aparasq      ;
      ObservableRef AperpsqName  , AzerosqName                 ;
      double        deltaperp    , deltazero    , deltapara    ;
      ObservableRef deltaperpName, deltazeroName, deltaparaName;
      // Options
      vector<string> componentlist;
    private:
      // Calculation
      Bs2PhiKKComponent* Pwave;
      void SetComponentAmplitudes();
      double PhaseSpace(double);
      // Stuff to do on creation
      void Initialise();
      void MakePrototypes();
      // Options
      string compName;
};
#endif

