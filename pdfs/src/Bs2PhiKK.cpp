#include "Bs2PhiKK.h"
#include "DPHelpers.hh"
#include "StringProcessing.h"

Bs2PhiKK::Bs2PhiKK(PDFConfigurator* config)
{
	std::vector<std::pair<dim, std::string>> wantednames = {{_mKK_, "mKK"}, {_phi_, "phi"}, {_ctheta_1_, "ctheta_1"}, {_ctheta_2_, "ctheta_2"}, {_trigger_, "trigger"}};
	// Dependent variable names
	for(const auto& psbname: config->GetPhaseSpaceBoundary()->GetAllNames())
		for(const auto& name: wantednames)
			if(psbname == name.second)
				ObservableNames.emplace(name.first, config->getName(name.second));
}
Bs2PhiKK::Bs2PhiKK(const Bs2PhiKK& copy)
	// Dependent variable names
	: ObservableNames(copy.ObservableNames)
{
}
void Bs2PhiKK::PhysPar::Update(const ParameterSet* pars)
{
	value = pars->GetPhysicsParameter(name)->GetValue();
	if(std::isnan(value))
		std::cerr << name.Name() << " has been given a nan value!" << std::endl;
}
void Bs2PhiKK::MakePrototypeDataPoint(std::vector<std::string>& allObservables)
{
	// The ordering here matters. It has to be the same as the XML file, apparently.
	// The order should be as in the enum declaration in this class.
	std::cout << "Prototype datapoint: ";
	for(const auto& name: ObservableNames)
	{
		std::cout << name.second.Name();
		allObservables.push_back(name.second.Name());
	}
	std::cout << std::endl;
}
Bs2PhiKK::datapoint_t Bs2PhiKK::ReadDataPoint(DataPoint* measurement) const
{
	// Get values from the datapoint
	datapoint_t dp;
	for(const auto& name: ObservableNames)
		dp[name.first] = measurement->GetObservable(name.second)->GetValue();
	if(ObservableNames.find(_phi_) != ObservableNames.end())
		dp[_phi_]+=M_PI; // TODO: get rid of the need for this
	return dp;
}

bool Bs2PhiKK::IsPhysicalDataPoint(const Bs2PhiKK::datapoint_t& datapoint)
{
	bool isphysical = true;
	if(datapoint.find(_mKK_) != datapoint.end())
		isphysical = isphysical && datapoint.at(_mKK_) > 2*Bs2PhiKK::mK;
	for(const auto& angle : {_ctheta_1_, _ctheta_2_})
		if(datapoint.find(angle) != datapoint.end())
			isphysical = isphysical && std::abs(datapoint.at(angle)) <= 1;
	return isphysical;
}

Bs2PhiKK::datapoint_t Bs2PhiKK::Parity(const Bs2PhiKK::datapoint_t& datapoint)
{
	datapoint_t returnval = datapoint;
	for(const auto& angle : {_phi_, _ctheta_1_, _ctheta_2_})
		if(returnval.find(angle) != returnval.end())
			returnval[angle] *= -1;
	return returnval;
}

std::vector<std::string> Bs2PhiKK::LineShapeParameterNames(std::string name, std::string lineshape)
{
	// Breit Wigner
	if(lineshape == "BW")
		return {name+"_mass", name+"_width", "KKBFradius"};
	// Flatte
	else if(lineshape == "FT")
		return {name+"_mass", name+"_gpipi", name+"_Rg"};
	else if(lineshape == "SP")
	{
		std::vector<std::string> parnames;
		for(const auto& KKname: StringProcessing::SplitString(name, ':'))
		{
			if(KKname == "")
				continue;
			parnames.push_back(KKname+"_mass"); // x-value of spline knot
			parnames.push_back(KKname+"_fraction"); // magnitude of y-value of spline knot
			parnames.push_back(KKname+"_deltazero"); // phase of y-value of spline knot
		}
		return parnames;
	}
	// Assume it's non-resonant
	else
		return {};
}
void Bs2PhiKK::UpdateLineshape(const std::string& lineshape, DPMassShape& KKLineShape, const std::vector<PhysPar>& KKpars)
{
	// Update the resonance line shape
	vector<double> respars;
	if(lineshape == "BW")
	{
		respars.push_back(KKpars[0].value); // mass
		respars.push_back(KKpars[1].value); // width
		respars.push_back(KKpars[2].value); // barrier factor radius
	}
	else if(lineshape == "FT")
	{
		respars.push_back(KKpars[0].value); // mass
		respars.push_back(KKpars[1].value); // gpipi
		respars.push_back(KKpars[1].value*KKpars[2].value); // gKK = gpipi*Rg
	}
	else if(lineshape == "SP")
	{
		for(const auto& par: KKpars)
			respars.push_back(par.value);
	}
	KKLineShape.setParameters(respars);
}

