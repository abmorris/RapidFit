/** @class Bs2PhiPhi_Adam Bs2PhiPhi_Adam.h
 *
 *  RapidFit PDF for Bs2PhiPhi_Adam
 *
 *  @author Adam Morris
 *  @date Aug 2016
 */
#ifndef Bs2PhiPhi_Adam_H
#define Bs2PhiPhi_Adam_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2PhiPhi_Adam : public BasePDF
{
    public:
      // *structors
      Bs2PhiPhi_Adam(PDFConfigurator*);
      Bs2PhiPhi_Adam(const Bs2PhiPhi_Adam&);
      ~Bs2PhiPhi_Adam();
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
      // Parameters
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

