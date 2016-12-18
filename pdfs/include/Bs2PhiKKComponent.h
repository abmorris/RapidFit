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
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include "DPWignerFunctionGeneral.hh"
using std::string;
using std::vector;

// Simplify the case where a value and a name correspond 1:1
struct PhysPar
{
  // Construct this however you want
  PhysPar() {}
  PhysPar(ObservableRef _name) : name(_name), value(0) {}
  PhysPar(ObservableRef _name, double _value) : name(_name), value(_value) {}
  PhysPar(PDFConfigurator* config, string _name) : name(config->getName(_name)), value(0) { cout << "Will look for parameter " << _name << endl; }
  PhysPar(PDFConfigurator* config, string _name, double _value) : name(config->getName(_name)), value(_value) {}
  PhysPar(const PhysPar& other) : value(other.value), name(other.name) {}
  void Update(ParameterSet& pars) { value = pars.GetPhysicsParameter(name)->GetValue(); }
  void Update(ParameterSet* pars) { value = pars->GetPhysicsParameter(name)->GetValue(); }
  double value;
  ObservableRef name;
};
class Bs2PhiKKComponent
{
  public:
    Bs2PhiKKComponent(PDFConfigurator*, string, string, int, string); // config, phi name, resonance name, spin
    Bs2PhiKKComponent(const Bs2PhiKKComponent&);
    ~Bs2PhiKKComponent();
    string GetName() {return KKname;}
    void SetPhysicsParameters(ParameterSet* pars) { SetPhysicsParameters(*pars); } // Update all the parameters
    void SetPhysicsParameters(ParameterSet&); // Update all the parameters
    vector<ObservableRef> GetPhysicsParameters();
    TComplex Amplitude(double, double, double, double); // KK_M, Phi_angle, cos_theta1, cos_theta2
    TComplex Amplitude(double, double, double, double, string);
    static double mBs;
    static double mK;
    static double mpi;
  protected:
    // Floatable parameters
    PhysPar fraction; // Unnormalised variable to control the relative contribution of each resonance. Do not use at the fit fraction!!
    vector<TComplex> Ahel;  // Helicity amplitude(s)
    vector<int>      helicities; // Store the possible values of helicity to enable looping over A(helicities[i])
    // Polarisation amplitude components (perp, zero, para)
    vector<PhysPar> magsqs; // Square of magnitudes: para will be calculated from the other two
    vector<PhysPar> phases; // Phases
    // Resonance parameters
    vector<PhysPar> phipars; // Mass and width of Breit Wigner
    vector<PhysPar> KKpars; // Mass and width of Breit Wigner, or mass, g_pipi and R=(g_KK/g_pipi) of Flatte. Empty for non-resonant
    // Fixed parameters
    int    Jphi; // Spin of the phi (P-wave, 1)
    int    JKK; // Spin of the KK resonance (0, 1 or 2)
    double RBs; // Bs barrier factor radius
    double RKK; // KK barrier factor radius
    string lineshape; // Choose the resonance shape: "BW", "FT" or "NR"
    string phiname; // The name decides which set of PhysicsParameters it will look for in the RapidFit XML
    string KKname;
    // Resonance lineshape function for the mass-dependent part
    DPMassShape* KKLineShape;
  private:
    void     Initialise();
    TComplex A(int); // Polarisation amplitude coefficients
    TComplex F(int, double, double, double); // Angular distribution: helicity, phi, costheta1, costheta2
    double   OFBF(double); // Product of orbital and barrier factors
    // Wigner d-functions for the angular-dependent part
    DPWignerFunction* wignerKK;
    DPWignerFunction* wignerPhi;
    // Blatt-Weisskopf barrier penetration factors
    DPBarrierFactor*  Bsbarrier;
    DPBarrierFactor*  KKbarrier;
};
#endif

