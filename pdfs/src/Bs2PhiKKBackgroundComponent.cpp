#include "Bs2PhiKKBackgroundComponent.h"
Bs2PhiKKBackgroundComponent::Bs2PhiKKBackgroundComponent(PDFConfigurator* config, std::string name, std::string _type)
	: fraction(Bs2PhiKK::PhysPar(config,name+"_fraction"))
	, type(_type)
	, angulardistribution(LegendreMomentShape(config->getConfigurationValue(name+"_CoefficientsFile")))
{
	if(type == "peaking")
	{
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_mean"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_sigma"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_alpha"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_n"));
	}
	else
	{
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_A"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_B"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_C"));
		shapepars.push_back(Bs2PhiKK::PhysPar(config,name+"_M"));
	}
}
Bs2PhiKKBackgroundComponent::Bs2PhiKKBackgroundComponent(const Bs2PhiKKBackgroundComponent& other)
	: fraction(other.fraction)
	, type(other.type)
	, shapepars(other.shapepars)
	, angulardistribution(other.angulardistribution)
{
}
//Bs2PhiKKBackgroundComponent& Bs2PhiKKBackgroundComponent::operator=(const Bs2PhiKKBackgroundComponent& other)
//{
//	fraction = other.fraction;
//	type = other.type;
//	shapepars = other.shapepars;
//	angulardistribution = other.angulardistribution;
//	return *this;
//}
double Bs2PhiKKBackgroundComponent::Evaluate(const Bs2PhiKK::datapoint_t& datapoint) const
{
	double massPart(0);
//	double result(0);
//	double arg = mKK - M;
//	if(arg <= 0) return 0;
//	double ratio = mKK/M;
//	double val = (1- exp(-arg/C))* pow(ratio, A) + B*(ratio-1);
//	result = val > 0 ? val : 0;
//	result *= shape->Evaluate(mKK, phi, ctheta_1, ctheta_2);
//	return result/231.634;
	
	
	double angularPart = angulardistribution.Evaluate(datapoint);
	return massPart * angularPart;
}
// Update everything from the parameter set
void Bs2PhiKKBackgroundComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
	fraction.Update(fitpars);
	for(auto& par: shapepars)
		par.Update(fitpars);
}
vector<ObservableRef> Bs2PhiKKBackgroundComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {{fraction},shapepars})
		for(const auto& par: set)
			parameters.push_back(par.name);
	return parameters;
}

