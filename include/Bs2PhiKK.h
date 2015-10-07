/** @class Bs2PhiKK Bs2PhiKK.h
 *
 *  RapidFit PDF for Bs2PhiKK
 *
 *  @author Adam Morris
 *  @date Sept 2015
 */
#ifndef Bs2PhiKK_H
#define Bs2PhiKK_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "TMath.h"

class Bs2PhiKK : public BasePDF
{
    public:
      Bs2PhiKK(PDFConfigurator*);
      Bs2PhiKK(const Bs2PhiKK&);
      ~Bs2PhiKK();

      //Calculate the PDF value
      virtual double Evaluate(DataPoint*);
      virtual bool SetPhysicsParameters(ParameterSet*);
      virtual vector<string> GetDoNotIntegrateList();
    protected:
      virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);
      double Azero2, Apara2, GammaL, GammaH, deltapara;
      bool init;
    private:
      Bs2PhiKK& operator=( const Bs2PhiKK& );
      void MakePrototypes();
      ObservableRef Apara2Name, Azero2Name, GammaLName, GammaHName, deltaparaName;
      // These contain the ObservableRefs that correspond
      // to the observable names that are used in the
      // PDF.
      ObservableRef ctheta_1Name;  // angle between Phi1 p vector and its K+ daughter's p vector
      ObservableRef ctheta_2Name;  // angle between Phi2 p vector and its K+ daughter's p vector
      ObservableRef phiName;       // angle between the decay planes of the 2 Phi resonances
      // Observables as member variables
      double ctheta_1, ctheta_2, phi;
      // Functions required for helicity
      inline double cHt1sq() const { return (ctheta_1*ctheta_1) ; }
      inline double cHt2sq() const { return (ctheta_2*ctheta_2) ; }
      inline double sHt1sq() const { return (sin(acos(ctheta_1))*sin(acos(ctheta_1))) ; }
      inline double sHt2sq() const { return (sin(acos(ctheta_2))*sin(acos(ctheta_2))) ; }
      inline double c2Ht1()  const { return (cos(2.*acos(ctheta_1))) ; }
      inline double c2Ht2()  const { return (cos(2.*acos(ctheta_2))) ; }
      inline double s2Ht1()  const { return (sin(2.*acos(ctheta_1))) ; }
      inline double s2Ht2()  const { return (sin(2.*acos(ctheta_2))) ; }
      inline double cHphi()  const { return (cos(phi)); }
      inline double sHphi()  const { return (sin(phi)); }
      inline double c2Hphi() const { return (cos(2.*phi)); }
      inline double s2Hphi() const { return (sin(2.*phi)); }
      inline double GF()     const { return 9./32./TMath::Pi() ; }
      // Angle factor for three angle PDFs in helicity basis
      inline double HangleFactorA0A0()   const { return 4. * cHt1sq() * cHt2sq() * GF(); }
      inline double HangleFactorAPAP()   const { return sHt1sq() * sHt2sq() * (1.+c2Hphi()) * GF(); }
      inline double HangleFactorATAT()   const { return sHt1sq() * sHt2sq() * (1.-c2Hphi()) * GF(); }
      inline double HangleFactorImAPAT() const { return -2. * sHt1sq() * sHt2sq() * s2Hphi() * GF(); }
      inline double HangleFactorReA0AP() const { return sqrt(2.) * s2Ht1() * s2Ht2() * cHphi() * GF(); }
      inline double HangleFactorImA0AT() const { return -1. * sqrt(2.) * s2Ht1() * s2Ht2() * sHphi() * GF(); }
      // Amplitudes
      double Aperp2() { return 1.0 - Azero2 - Apara2 ; }
      double K1() { return Azero2/GammaL   ; }
      double K2() { return Apara2/GammaL   ; }
      double K3() { return Aperp2()/GammaH ; }
      double K4() { return 0; }
      double K5() { return sqrt(Azero2)*sqrt(Apara2)*cos(deltapara)/GammaL ; }
      double K6() { return 0; }
};
#endif
