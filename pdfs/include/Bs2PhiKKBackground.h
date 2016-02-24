/** @class Bs2PhiKKBackground Bs2PhiKKBackground.h
 *
 *  RapidFit PDF for Bs2PhiKKBackground
 *
 *  @author Adam Morris
 *  @date Feb 2016
 */
#ifndef Bs2PhiKKBackground_H
#define Bs2PhiKKBackground_H

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

class Bs2PhiKKBackground : public BasePDF
{
    public:
      // *structors
      Bs2PhiKKBackground(PDFConfigurator*);
      Bs2PhiKKBackground(const Bs2PhiKKBackground&);
      ~Bs2PhiKKBackground();
      // Required methods
      double Evaluate(DataPoint*);
      double EvaluateComponent(DataPoint*,ComponentRef*);
      bool SetPhysicsParameters(ParameterSet*);
      vector<string> GetDoNotIntegrateList();
    protected:
      // Shape parameters
      double        A,     B,     C,     M;
      ObservableRef AName, BName, CName, MName;
      // K+K− mass and helicity angles
      double        mKK,     ctheta_1,     ctheta_2,     phi;
      ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName;
      // Acceptance object
      LegendreMomentShape* shape;
    private:
      string filename;
      // Stuff to do on creation
      void Initialise();
      void MakePrototypes();
};
#endif

