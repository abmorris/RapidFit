#ifndef Bs2PhiKKPeakingBackground_H
#define Bs2PhiKKPeakingBackground_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
// C++
#include <string>
// Self
#include "LegendreMomentShape.h"

class Bs2PhiKKPeakingBackground : public BasePDF
{
    public:
      // *structors
      Bs2PhiKKPeakingBackground(PDFConfigurator*);
      Bs2PhiKKPeakingBackground(const Bs2PhiKKPeakingBackground&);
      ~Bs2PhiKKPeakingBackground();
      // Required methods
      double Evaluate(DataPoint*);
      double EvaluateComponent(DataPoint*,ComponentRef*);
      bool SetPhysicsParameters(ParameterSet*);
      vector<string> GetDoNotIntegrateList();
      vector<string> PDFComponents();
    protected:
      // Shape parameters
      double        mean,     sigma,     alpha,     n;
      ObservableRef meanName, sigmaName, alphaName, nName;
      // K+K− mass and helicity angles
      double        mKK,     ctheta_1,     ctheta_2,     phi;
      ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
      // Acceptance object
      LegendreMomentShape* shape;
    private:
      string mode;
      // Stuff to do on creation
      void Initialise();
      void MakePrototypes();
};
#endif

