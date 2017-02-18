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
double Bs2PhiKKBackgroundComponent::Evaluate(const Bs2PhiKK::datapoint_t& datapoint) const
{
	double massPart(0);
	double mKK = datapoint[0];
	if(type == "peaking")
	{
		double mean  = shapepars[0].value;
		double sigma = shapepars[1].value;
		double alpha = shapepars[2].value;
		if((mKK-mean)/sigma>-alpha)
		{
			massPart = std::exp(-std::pow(mKK-mean,2)/(2*sigma*sigma));
		}
		else
		{
			double n = shapepars[3].value;
			double A = std::pow(n/std::abs(alpha),n) * std::exp(-alpha*alpha/2);
			double B = n/std::abs(alpha) - std::abs(alpha);
			massPart = A * std::pow(B - (mKK-mean)/sigma,-n);
		}
	}
	else
	{
		double A = shapepars[0].value;
		double B = shapepars[1].value;
		double C = shapepars[2].value;
		double M = shapepars[3].value;
		double arg = mKK - M;
		if(arg <= 0) return 0;
		double ratio = mKK/M;
		double val = (1- exp(-arg/C))* pow(ratio, A) + B*(ratio-1);
		massPart = val > 0 ? val : 0;
	}
	
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

