/** @class Bs2PhiKKTotal Bs2PhiKKTotal.h
 *
 *  RapidFit PDF for Bs2PhiKKTotal
 *
 *  @author Adam Morris
 *  @date Feb 2016
 */
#ifndef Bs2PhiKKTotal_H
#define Bs2PhiKKTotal_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
// Self
#include "Bs2PhiKKComponent.h"
#include "LegendreMomentShape.h"

class Bs2PhiKKTotal : public BasePDF
{
    public:
      // *structors
      Bs2PhiKKTotal(PDFConfigurator*);
      Bs2PhiKKTotal(const Bs2PhiKKTotal&);
      ~Bs2PhiKKTotal();
      // Required methods
      double Evaluate(DataPoint*);
      double Normalisation(PhaseSpaceBoundary*);
      bool SetPhysicsParameters(ParameterSet*);
      // Extra stuff
      double EvaluateComponent( DataPoint*, ComponentRef* );
      vector<string> PDFComponents();
    protected:
      // K+K− mass and helicity angles
      double        mKK,     ctheta_1,     ctheta_2,     phi;
      ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
      // Non-resonant component
      double        ANonRes;
      ObservableRef ANonResName;
      // Magnitude-squared of helicity amplitudes
      double        ASsq,     APsq[3],     ADsq[3];
      ObservableRef ASsqName, APsqName[3], ADsqName[3];
      // Phase of helicity amplitudes
      double        deltaS,     deltaP[3],     deltaD[3];
      ObservableRef deltaSName, deltaPName[3], deltaDName[3];
      // m(KK) boundaries
      double mKKmin, mKKmax;
      // Acceptance object
      LegendreMomentShape* acc;
    private:
      // The m(KK) components
      Bs2PhiKKComponent* Swave;
      Bs2PhiKKComponent* Pwave;
      Bs2PhiKKComponent* Dwave;
      Bs2PhiKKComponent* NonRes;
      vector<string> componentlist;
      // Stuff to do on creation
      void Initialise();
      bool init;
      void MakePrototypes();
      // Calculation
      TComplex TotalAmplitude(double, double, double, double);
      double Acceptance(double, double, double, double);
      void ReadDataPoint(DataPoint*);
      void SetComponentAmplitudes();
      int compIndex;
      bool debug = false;
};
#endif

