/** @class Bs2PhiKKBackground Bs2PhiKKBackground.cpp
 *
 *  RapidFit PDF for Bs2PhiKKBackground
 *
 *  @author Adam Morris
 *  @date Feb 2016
 */
#include "Bs2PhiKKBackground.h"
PDF_CREATOR( Bs2PhiKKBackground )

Bs2PhiKKBackground::Bs2PhiKKBackground(PDFConfigurator* config) :
  // Dependent variables
    mKK(0.0)
  , phi(0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Dependent variable names
  , mKKName       ( config->getName("mKK"     ) )
  , ctheta_1Name  ( config->getName("ctheta_1") )
  , ctheta_2Name  ( config->getName("ctheta_2") )
  , phiName       ( config->getName("phi"     ) )
  // Physics parameters
  , A(0.0)
  , B(0.0)
  , C(0.0)
  // Physics parameter names
  , AName ( config->getName("A") )
  , BName ( config->getName("B") )
  , CName ( config->getName("C") )
{
  shape = new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile"));
}
Bs2PhiKKBackground::Bs2PhiKKBackground(const Bs2PhiKKBackground& copy) :
    BasePDF( (BasePDF) copy)
  // Dependent variables
  , mKK(copy.mKK)
  , phi(copy.phi)
  , ctheta_1(copy.ctheta_1)
  , ctheta_2(copy.ctheta_2)
  // Dependent variable names
  , mKKName       (copy.mKKName)
  , ctheta_1Name  (copy.ctheta_1Name)
  , ctheta_2Name  (copy.ctheta_2Name)
  , phiName       (copy.phiName)
  // Physics parameters
  , A(copy.A)
  , B(copy.B)
  , C(copy.C)
  // Physics parameter names
  , AName(copy.AName)
  , BName(copy.BName)
  , CName(copy.CName)
{
  shape = new LegendreMomentShape(*copy.shape);
}
Bs2PhiKKBackground::~Bs2PhiKKBackground()
{
//  delete shape;
}
void Bs2PhiKKBackground::Initialise()
{

}
void Bs2PhiKKBackground::MakePrototypes()
{
  cout << "Making prototypes" << endl;
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
}
double Bs2PhiKKBackground::Evaluate(DataPoint* measurement)
{
  mKK      = measurement->GetObservable(mKKName     )->GetValue();
  ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
  phi      = measurement->GetObservable(phiName     )->GetValue();
  double result(0);
  result = shape->Evaluate(mKK, phi, ctheta_1, ctheta_2); // TODO: mass-dependent part!
  return result;
}
double Bs2PhiKKBackground::Normalisation(PhaseSpaceBoundary* boundary)
{
  double result(0);
  double mKKhi      = boundary->GetConstraint(mKKName)->GetMaximum(); 
  double mKKlo      = boundary->GetConstraint(mKKName)->GetMinimum();
  double phihi      = boundary->GetConstraint(phiName)->GetMaximum();
  double philo      = boundary->GetConstraint(phiName)->GetMinimum();
  double ctheta_1hi = boundary->GetConstraint(ctheta_1Name)->GetMaximum();
  double ctheta_1lo = boundary->GetConstraint(ctheta_1Name)->GetMinimum();
  double ctheta_2hi = boundary->GetConstraint(ctheta_2Name)->GetMaximum();
  double ctheta_2lo = boundary->GetConstraint(ctheta_2Name)->GetMinimum();
  result = shape->Integral(mKKhi, mKKlo, phihi, philo, ctheta_1hi, ctheta_1lo, ctheta_2hi, ctheta_2lo); // TODO: mass-dependent part!
  return result;
}
bool Bs2PhiKKBackground::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  A = allParameters.GetPhysicsParameter(AName)->GetValue();
  B = allParameters.GetPhysicsParameter(BName)->GetValue();
  C = allParameters.GetPhysicsParameter(CName)->GetValue();
  return isOK;
}
vector<string> Bs2PhiKKBackground::GetDoNotIntegrateList()
{
  vector<string> list;
  return list;
}
