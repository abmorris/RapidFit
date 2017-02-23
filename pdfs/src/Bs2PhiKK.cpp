#include "Bs2PhiKK.h"
#include "DPHelpers.hh"

double Bs2PhiKK::mBs  = 5.36677;
double Bs2PhiKK::mK   = 0.493677;
double Bs2PhiKK::mpi  = 0.139570;

Bs2PhiKK::Bs2PhiKK(PDFConfigurator* config)
	// Dependent variable names
	: mKKName(config->getName("mKK"))
	, phiName(config->getName("phi"))
	, ctheta_1Name(config->getName("ctheta_1"))
	, ctheta_2Name(config->getName("ctheta_2"))
{
}
Bs2PhiKK::Bs2PhiKK(const Bs2PhiKK& copy)
	// Dependent variable names
	: mKKName(copy.mKKName)
	, phiName(copy.phiName)
	, ctheta_1Name(copy.ctheta_1Name)
	, ctheta_2Name(copy.ctheta_2Name)
{
}
void Bs2PhiKK::PhysPar::Update(const ParameterSet* pars)
{
	value = pars->GetPhysicsParameter(name)->GetValue();
	if(std::isnan(value))
		std::cerr << name.Name() << " has been given a nan value!" << std::endl;
}
Bs2PhiKK::datapoint_t Bs2PhiKK::ReadDataPoint(DataPoint* measurement) const
{
	// Get values from the datapoint
	double mKK      = measurement->GetObservable(mKKName     )->GetValue();
	double phi      = measurement->GetObservable(phiName     )->GetValue();
	double ctheta_1 = measurement->GetObservable(ctheta_1Name)->GetValue();
	double ctheta_2 = measurement->GetObservable(ctheta_2Name)->GetValue();
	phi+=M_PI;
	return {mKK, phi, ctheta_1, ctheta_2};
}

bool Bs2PhiKK::IsPhysicalDataPoint(const Bs2PhiKK::datapoint_t& datapoint)
{
	return datapoint[0] > 2*Bs2PhiKK::mK && std::abs(datapoint[2]) <= 1 && std::abs(datapoint[3]) <= 1;
}

std::vector<std::string> Bs2PhiKK::LineShapeParameterNames(std::string name, std::string lineshape)
{
	// Breit Wigner
	if(lineshape=="BW")
	{
		return {name+"_mass", name+"_width", "KKBFradius"};
	}
	// Flatte
	else if(lineshape=="FT")
	{
		return {name+"_mass", name+"_gpipi", name+"_Rg"};
	}
	// Assume it's non-resonant
	else
	{
		return {};
	}
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
	KKLineShape.setParameters(respars);
}

