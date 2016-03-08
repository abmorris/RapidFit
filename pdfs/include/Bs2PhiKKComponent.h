/** @class Bs2PhiKKComponent Bs2PhiKKComponent.cpp
 *
 *  RapidFit PDF for Bs2PhiKKComponent
 *
 *  @author Adam Morris
 *  @date Nov-Dec 2015
 */
#ifndef __BS2PHIKKCOMPONENT_H__
#define __BS2PHIKKCOMPONENT_H__
// ROOT
#include "TComplex.h"
// Std
#include <string>
#include <vector>
// RapidFit
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPWignerFunctionGeneral.hh"
using std::string;
using std::vector;
class Bs2PhiKKComponent
{
  public:
    Bs2PhiKKComponent(int, double, double, string, double, double); // J2, M2, W2, shape, RBs, RKK
    Bs2PhiKKComponent(const Bs2PhiKKComponent&);
    ~Bs2PhiKKComponent();
    void SetHelicityAmplitudes(int, double, double); 
    TComplex Amplitude(double, double, double, double); // KK_M, Phi_angle, cos_theta1, cos_theta2
    void Print();
    static double mBs ;
    static double mphi;
    static double mK  ;
    static double mpi ;
  protected:
    double*           _Amag;  // Arrays to contain the amplitudes.
    double*           _Aphase;
    int               _J1; // Spin of the phi (P-wave, 1)
    int               _J2; // Spin of the KK resonance (0, 1 or 2)
    double            _M1; // Mass of the phi
    double            _M2; // Mass of the KK resonance
    double            _W1; // Width of the phi
    double            _W2; // Width of the KK resonance
    double            _RBs; // Bs barrier factor radius
    double            _RKK; // KK barrier factor radius
    string            _shape; // Choose the resonance shape
    DPWignerFunction* wigner;
    DPWignerFunction* wignerPhi;
  private:
    void              Initialise();
    TComplex          A(int);                    // Polarisation amplitude coefficients
    TComplex          F(double, double, double); // Angular part
    TComplex          M(double);                 // Mass-dependent part (KK resonance shape)
    double            OFBF(double);
    DPBarrierFactor*  Bsbarrier; // Blatt-Weisskopf barrier penetration factor for the Bs
    DPBarrierFactor*  KKbarrier; // Barrier factor for the KK resonance
    DPMassShape*      _M; // Pointer to resonance shape function
    int               _lambda_max; // Keep track of the max. helicity value
    bool debug = false;
};
#endif

