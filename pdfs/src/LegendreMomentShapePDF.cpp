#include "LegendreMomentShapePDF.h"
PDF_CREATOR( LegendreMomentShapePDF )

LegendreMomentShapePDF::LegendreMomentShapePDF(PDFConfigurator* config) :
	// Dependent variables
	  mKK(0.0)
	, phi(0.0)
	, ctheta_1(0.0)
	, ctheta_2(0.0)
	// Dependent variable names
	, mKKName       ( config->getName("mKK"     ) )
	, phiName       ( config->getName("phi"     ) )
	, ctheta_1Name  ( config->getName("ctheta_1") )
	, ctheta_2Name  ( config->getName("ctheta_2") )
{
	shape = new LegendreMomentShape(config->getConfigurationValue("CoefficientsFile"));
}
LegendreMomentShapePDF::LegendreMomentShapePDF(const LegendreMomentShapePDF& copy) : BasePDF( (BasePDF) copy)
	// Dependent variables
	, mKK(copy.mKK)
	, phi(copy.phi)
	, ctheta_1(copy.ctheta_1)
	, ctheta_2(copy.ctheta_2)
	// Dependent variable names
	, mKKName       (copy.mKKName)
	, phiName       (copy.phiName)
	, ctheta_1Name  (copy.ctheta_1Name)
	, ctheta_2Name  (copy.ctheta_2Name)
{
	shape = new LegendreMomentShape(*copy.shape);
}
LegendreMomentShapePDF::~LegendreMomentShapePDF()
{
//  delete shape;
}
void LegendreMomentShapePDF::Initialise()
{
	// Enable numerical normalisation and disable caching
	// There should be an analytical integral, but it apparently somehow involves gamma functions of non-positive integers???
	this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
}
void LegendreMomentShapePDF::MakePrototypes()
{
	cout << "Making prototypes" << endl;
	// Make the DataPoint prototype
	// The ordering here matters. It has to be the same as the XML file, apparently.
	allObservables.push_back(mKKName     );
	allObservables.push_back(phiName     );
	allObservables.push_back(ctheta_1Name);
	allObservables.push_back(ctheta_2Name);
}
double LegendreMomentShapePDF::Evaluate(DataPoint* measurement)
{
	mKK      = measurement->GetObservable(mKKName     )->GetValue();
	phi      = measurement->GetObservable(phiName     )->GetValue();
	ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
	ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
	double result(0);
	result = shape->Evaluate(mKK, phi, ctheta_1, ctheta_2);
	return result;
}
double LegendreMomentShapePDF::Normalisation(PhaseSpaceBoundary* boundary)
{
	(void)boundary;
	return -1;
}
bool LegendreMomentShapePDF::SetPhysicsParameters(ParameterSet* parameters)
{
	(void)parameters;
	return true;
}
vector<string> LegendreMomentShapePDF::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}
