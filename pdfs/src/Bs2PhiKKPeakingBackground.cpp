/** @class Bs2PhiKKPeakingBackground Bs2PhiKKPeakingBackground.cpp
 *
 *  RapidFit PDF for Bs2PhiKKPeakingBackground
 *
 *  @author Adam Morris
 *  @date Mar 2016
 */
#include "Bs2PhiKKPeakingBackground.h"
PDF_CREATOR( Bs2PhiKKPeakingBackground )

Bs2PhiKKPeakingBackground::Bs2PhiKKPeakingBackground(PDFConfigurator* config) :
  // Dependent variables
    mKK     (0.0)
  , phi     (0.0)
  , ctheta_1(0.0)
  , ctheta_2(0.0)
  // Physics parameters
  , mean    (0.0)
  , sigma   (0.0)
  , alpha   (0.0)
  , n       (0.0)
  // Dependent variable names
  , mKKName     (config->getName("mKK"     ))
  , ctheta_1Name(config->getName("ctheta_1"))
  , ctheta_2Name(config->getName("ctheta_2"))
  , phiName     (config->getName("phi"     ))
{
  MakePrototypes();
  // Get the prefix for the physics parameters in the configuration file
  // This allows for the sum of multiple instances of the same PDF with different parameters
  mode = config->getConfigurationValue("DecayMode");
  // Physics parameter names
  meanName  = config->getName(mode+"mean" )
  sigmaName = config->getName(mode+"sigma")
  alphaName = config->getName(mode+"alpha")
  nName     = config->getName(mode+"n"    )
  // Create the angular function
  shape = new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile"));
  this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
}
Bs2PhiKKPeakingBackground::Bs2PhiKKPeakingBackground(const Bs2PhiKKPeakingBackground& copy) :
    BasePDF( (BasePDF) copy)
  // Dependent variables
  , mKK           (copy.mKK         )
  , phi           (copy.phi         )
  , ctheta_1      (copy.ctheta_1    )
  , ctheta_2      (copy.ctheta_2    )
  // Dependent variable names
  , mKKName       (copy.mKKName     )
  , ctheta_1Name  (copy.ctheta_1Name)
  , ctheta_2Name  (copy.ctheta_2Name)
  , phiName       (copy.phiName     )
  // Physics parameters
  , mean          (copy.mean        )
  , sigma         (copy.sigma       )
  , alpha         (copy.alpha       )
  , n             (copy.n           )
  // Physics parameter names
  , meanName      (copy.meanName    )
  , sigmaName     (copy.sigmaName   )
  , alphaName     (copy.alphaName   )
  , nName         (copy.nName       )
{
  shape = new LegendreMomentShape(*copy.shape);
  this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
}
Bs2PhiKKPeakingBackground::~Bs2PhiKKPeakingBackground()
{
//  delete shape;
}
void Bs2PhiKKPeakingBackground::Initialise()
{

}
void Bs2PhiKKPeakingBackground::MakePrototypes()
{
  cout << "Making prototypes" << endl;
  // Make the DataPoint prototype
  // The ordering here matters. It has to be the same as the XML file, apparently.
  allObservables.push_back(mKKName     );
  allObservables.push_back(phiName     );
  allObservables.push_back(ctheta_1Name);
  allObservables.push_back(ctheta_2Name);
  // Make the parameter set
  vector<string> parameterNames;
  parameterNames.push_back(meanName    );
  parameterNames.push_back(sigmaName   );
  parameterNames.push_back(alphaName   );
  parameterNames.push_back(nName       );
  allParameters = *( new ParameterSet(parameterNames) );
}
double Bs2PhiKKPeakingBackground::EvaluateComponent(DataPoint* measurement,ComponentRef* component)
{
  (void)component;
  return Evaluate(measurement);
}
vector<string> Bs2PhiKKPeakingBackground::PDFComponents()
{
  vector<string> list;
  return list;
}
double Bs2PhiKKPeakingBackground::Evaluate(DataPoint* measurement)
{
  mKK      = measurement->GetObservable(mKKName     )->GetValue();
  phi      = measurement->GetObservable(phiName     )->GetValue();
  ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
  ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
  double result = 1; // TODO: insert CrystalBall(mKK | mean, sigma, alpha, n)
  result *= shape->Evaluate(mKK, phi, ctheta_1, ctheta_2);
  return result;
}
bool Bs2PhiKKPeakingBackground::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
  bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  mean  = allParameters.GetPhysicsParameter(meanName )->GetValue();
  sigma = allParameters.GetPhysicsParameter(sigmaName)->GetValue();
  alpha = allParameters.GetPhysicsParameter(alphaName)->GetValue();
  n     = allParameters.GetPhysicsParameter(nName    )->GetValue();
  return isOK;
}
vector<string> Bs2PhiKKPeakingBackground::GetDoNotIntegrateList()
{
  vector<string> list;
  return list;
}
