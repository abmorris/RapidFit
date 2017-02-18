// Std
#include <algorithm>
// RapidFit
#include "Bs2PhiKKBackground.h"
#include "StringProcessing.h"

PDF_CREATOR( Bs2PhiKKBackground )
/*****************************************************************************/
// Constructor
Bs2PhiKKBackground::Bs2PhiKKBackground(PDFConfigurator* config) : Bs2PhiKK(config)
{
	std::cout << "\nBuilding Bs → ϕ K+ K− background PDF\n\n";
	std::vector<std::string> componentlist = StringProcessing::SplitString(config->getConfigurationValue("components"), ' ');
	std::cout << "┏━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━┓\n";
	std::cout << "┃ Component     │ Type          ┃\n";
	std::cout << "┠───────────────┼───────────────┨\n";
	for(const auto& name: componentlist)
	{
		if(name=="") continue;
		components[name] = ParseComponent(config,name);
	}
	std::cout << "┗━━━━━━━━━━━━━━━┷━━━━━━━━━━━━━━━┛" << std::endl;
	MakePrototypes();
	Initialise();
}
// Copy constructor
Bs2PhiKKBackground::Bs2PhiKKBackground(const Bs2PhiKKBackground& copy)
	: BasePDF((BasePDF)copy)
	, Bs2PhiKK((Bs2PhiKK)copy)
	// PDF components
	, components(copy.components)
{
	Initialise();
}
void Bs2PhiKKBackground::Initialise()
{
	// Enable numerical normalisation and disable caching
	this->SetNumericalNormalisation( true );
//	this->TurnCachingOff();
}
/*****************************************************************************/
// Build a component object from a passed option
Bs2PhiKKBackgroundComponent Bs2PhiKKBackground::ParseComponent(PDFConfigurator* config, std::string option) const
{
	// Syntax: <background name>(<type>)
	// - the list is space-delimited: no extra spaces, please!
	// - background name must be alphanumeric
	// - type determines the shape of the mass-dependent part
	// - see Bs2PhiKKBackgroundComponent.cpp for implemented types (combinatorial, peaking... but more can be added)
	// Example: "BdPhiKstar(peaking)"
	size_t openbracket  = option.find('(');
	size_t closebracket = option.find(')');
	std::string name = option.substr(0,openbracket);
	std::string type = option.substr(openbracket+1,closebracket-openbracket-1);
	std::cout << "┃ " << name << "\t│ " << type << "\t┃\n";
	return Bs2PhiKKBackgroundComponent(config, name, type);
}
/*****************************************************************************/
void Bs2PhiKKBackground::MakePrototypes()
{
	// Make the DataPoint prototype
	// The ordering here matters. It has to be the same as the XML file, apparently.
	allObservables.push_back(mKKName);
	allObservables.push_back(phiName);
	allObservables.push_back(ctheta_1Name);
	allObservables.push_back(ctheta_2Name);
	// Make the parameter set
	std::vector<std::string> parameterNames;
	for(const auto& comp: components)
		for(std::string par: comp.second.GetPhysicsParameters())
			parameterNames.push_back(par);
	std::sort(parameterNames.begin(),parameterNames.end());
	parameterNames.erase(std::unique(parameterNames.begin(),parameterNames.end()),parameterNames.end());
	allParameters = ParameterSet(parameterNames);
}
/*****************************************************************************/
bool Bs2PhiKKBackground::SetPhysicsParameters(ParameterSet* NewParameterSet)
{
	UnsetCache();
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	for(auto& comp: components)
		comp.second.SetPhysicsParameters(&allParameters);
	return isOK;
}
/*****************************************************************************/
vector<string> Bs2PhiKKBackground::PDFComponents()
{
	// Avoid redundant plotting for single-component PDFs
	if(components.size() == 1) return {};
	std::vector<std::string> componentnames;
	for(const auto& comp : components)
		componentnames.push_back(comp.first);
	return componentnames;
}
/*****************************************************************************/
double Bs2PhiKKBackground::EvaluateComponent(DataPoint* measurement,ComponentRef* component)
{
	std::string compName = component->getComponentName();
	if(compName=="0") // should quicken things up slightly?
		return Evaluate(measurement);
	const Bs2PhiKK::datapoint_t datapoint = ReadDataPoint(measurement);
	return components.at(compName).Evaluate(datapoint);
}
double Bs2PhiKKBackground::Evaluate(DataPoint* measurement)
{
	const Bs2PhiKK::datapoint_t datapoint = ReadDataPoint(measurement);
	double totalshape(0);
	for(const auto& comp : components)
		totalshape+=comp.second.Evaluate(datapoint);
	return totalshape;
}
/*****************************************************************************/
double Bs2PhiKKBackground::Normalisation(PhaseSpaceBoundary* boundary)
{
	(void)boundary;
	return -1;
}

