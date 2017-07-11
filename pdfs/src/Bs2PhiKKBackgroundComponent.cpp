#include <regex>
#include <utility>
#include "Mathematics.h"
#include "Bs2PhiKKBackgroundComponent.h"
#include "DPBWResonanceShape.hh"
#include "DPFlatteShape.hh"
#include "DPNonresonant.hh"
#include "TH1D.h"

// Constructor
Bs2PhiKKBackgroundComponent::Bs2PhiKKBackgroundComponent(PDFConfigurator* config, std::string name, std::string _type)
	: fraction(Bs2PhiKK::PhysPar(config,name+"_fraction"))
	, type(_type)
	, angulardistribution(LegendreMomentShape(config->getConfigurationValue(name+"_CoefficientsFile")))
{
	if(type == "peaking")
		for(std::string suffix: {"mean", "sigma", "alpha", "n"})
			KKpars.push_back(Bs2PhiKK::PhysPar(config,name+"_"+suffix));
	else if(type == "combinatorial")
		for(std::string suffix: {"A", "B", "C"})
			KKpars.push_back(Bs2PhiKK::PhysPar(config,name+"_"+suffix));
	else if(type == "histogram")
	{
		std::string histfilename = config->getConfigurationValue(name+"_HistFile");
		std::unique_ptr<TFile> histfile(TFile::Open(histfilename.c_str()));
		mKKhist = NDHist_Fixed(*(TH1D*)histfile->Get("mKKhist")); // Copy the histogram to the member object.
	}
	else
		std::cerr << "Bs2PhiKKBackgroundComponent WARNING: unrecognised background type:" << type << ". This will be treated as flat in mass.\n";
}
// Copy constructor
Bs2PhiKKBackgroundComponent::Bs2PhiKKBackgroundComponent(const Bs2PhiKKBackgroundComponent& other)
	: fraction(other.fraction)
	, type(other.type)
	, JKK(other.JKK)
	, KKpars(other.KKpars)
	, mKKhist(other.mKKhist)
	, angulardistribution(other.angulardistribution)
{
}
/*****************************************************************************/
double Bs2PhiKKBackgroundComponent::Evaluate(const Bs2PhiKK::datapoint_t& datapoint) const
{
	double massPart(1.0);
	double mKK = datapoint[Bs2PhiKK::_mKK_];
	if(type == "peaking")
	{
		// Crystal Ball function
		double mean  = KKpars[0].value;
		double sigma = KKpars[1].value;
		double alpha = KKpars[2].value;
		double n     = KKpars[3].value;
		double arg   = (mKK-mean)/sigma;
		if(alpha<0) arg *= -1;
		double absalpha = std::abs(alpha);
		if(arg >= -absalpha)
			massPart = std::exp(-0.5*arg*arg);
		else
		{
			double A = std::pow(n/absalpha,n) * std::exp(-0.5*absalpha*absalpha);
			double B = n/absalpha - absalpha;
			massPart = A / std::pow(B - arg,n);
		}
	}
	else if(type == "combinatorial")
	{
		// Threshold function
		double A = KKpars[0].value;
		double B = KKpars[1].value;
		double C = KKpars[2].value;
		double arg = mKK - 2*Bs2PhiKK::mK;
		if(arg <= 0) return 0;
		double ratio = mKK/2*Bs2PhiKK::mK;
		double val = (1- exp(-arg/C))* pow(ratio, A) + B*(ratio-1);
		massPart = val > 0 ? val : 0;
	}
	else if(type == "histogram")
		massPart = mKKhist.Eval({mKK})/mKKhist.Integral();
	double angularPart = angulardistribution.Evaluate({datapoint[Bs2PhiKK::_mKK_],datapoint[Bs2PhiKK::_phi_],datapoint[Bs2PhiKK::_ctheta_1_],datapoint[Bs2PhiKK::_ctheta_2_]});
	if(std::isnan(massPart)) std::cerr << "Mass part is nan" << std::endl;
	if(massPart < 0) std::cerr << "Mass part is negative" << std::endl;
	if(std::isnan(angularPart)) std::cerr << "Angular part is nan" << std::endl;
	if(angularPart < 0) std::cerr << "Angular part is negative" << std::endl;
	return fraction.value * massPart * angularPart;
}
/*****************************************************************************/
void Bs2PhiKKBackgroundComponent::SetPhysicsParameters(ParameterSet* fitpars)
{
	fraction.Update(fitpars);
	for(auto& par: KKpars)
	{
		par.Update(fitpars);
		if(std::isnan(par.value))
			std::cerr << par.name.Name() << " is nan" << std::endl;
	}
}
vector<ObservableRef> Bs2PhiKKBackgroundComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {{fraction},KKpars})
		for(const auto& par: set)
			parameters.push_back(par.name);
	return parameters;
}

