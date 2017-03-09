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
		for(std::string suffix: {"A", "B", "C", "M"})
			KKpars.push_back(Bs2PhiKK::PhysPar(config,name+"_"+suffix));
	else if(type.find("resonant")!=std::string::npos)
	{
		std::regex pattern("resonant,([0-3]),([A-Z][A-Z])");
		std::smatch result;
		bool match = std::regex_match(type, result, pattern);
		if(!match) std::cerr << "The option for this resonant background is malformed: " << type << std::endl;
		JKK = std::stoi(result[1].str());
		lineshape = result[2].str();
		string prefix = "combinatorial";
		for(string par: Bs2PhiKK::LineShapeParameterNames(name.substr(name.find(prefix)+prefix.length()),lineshape))
			KKpars.push_back(Bs2PhiKK::PhysPar(config,par));
		if(KKpars.empty() && lineshape!="NR")
		{
			lineshape = "NR";
			std::cerr << "Bs2PhiKKBackgroundComponent WARNING: unknown lineshape '" << lineshape << "'. Treating this component non-resonant.\n";
		}
	}
	else if(type == "histogram")
	{
		std::string histfilename = config->getConfigurationValue(name+"_HistFile");
		std::unique_ptr<TFile> histfile(TFile::Open(histfilename.c_str()));
		mKKhist = NDHist_Fixed(*(TH1D*)histfile->Get("mKKhist")); // Copy the histogram to the member object.
	}
	else
		std::cerr << "Bs2PhiKKBackgroundComponent WARNING: unrecognised background type:" << type << ". This will be treated as flat in mass.\n";
	Initialise();
}
// Copy constructor
Bs2PhiKKBackgroundComponent::Bs2PhiKKBackgroundComponent(const Bs2PhiKKBackgroundComponent& other)
	: fraction(other.fraction)
	, type(other.type)
	, lineshape(other.lineshape)
	, JKK(other.JKK)
	, KKpars(other.KKpars)
	, mKKhist(other.mKKhist)
	, angulardistribution(other.angulardistribution)
{
	Initialise();
}
// Copy by assignment
Bs2PhiKKBackgroundComponent& Bs2PhiKKBackgroundComponent::operator=(const Bs2PhiKKBackgroundComponent& other)
{
	fraction = other.fraction;
	type = other.type;
	lineshape = other.lineshape;
	JKK = other.JKK;
	KKpars = other.KKpars;
	mKKhist = other.mKKhist;
	angulardistribution = other.angulardistribution;
	Initialise();
	return *this;
}
void Bs2PhiKKBackgroundComponent::Initialise()
{
	// Breit Wigner
	if(lineshape=="BW")
		KKLineShape = std::unique_ptr<DPBWResonanceShape>(new DPBWResonanceShape(KKpars[0].value, KKpars[1].value, JKK, Bs2PhiKK::mK, Bs2PhiKK::mK, KKpars[2].value));
	// Flatte
	else if(lineshape=="FT")
		KKLineShape = std::unique_ptr<DPFlatteShape>(new DPFlatteShape(KKpars[0].value, KKpars[1].value, Bs2PhiKK::mpi, Bs2PhiKK::mpi, KKpars[1].value*KKpars[2].value, Bs2PhiKK::mK, Bs2PhiKK::mK));
	else if(lineshape=="NR")
		KKLineShape = std::unique_ptr<DPNonresonant>(new DPNonresonant());
	else return; // Don't update the lineshape if it doesn't exist
	Bs2PhiKK::UpdateLineshape(lineshape, *KKLineShape, KKpars);
}
/*****************************************************************************/
double Bs2PhiKKBackgroundComponent::Evaluate(const Bs2PhiKK::datapoint_t& datapoint) const
{
	double massPart(1.0);
	double mKK = datapoint[0];
	if(type == "peaking")
	{
		// Crystal Ball function
		double mean  = KKpars[0].value;
		double sigma = KKpars[1].value;
		double alpha = KKpars[2].value;
		double n     = KKpars[3].value;
		double arg = (mKK-mean)/sigma;
		if(alpha<0) arg *= -1;
		double absalpha = std::abs(alpha);
		double A,B;
		if(arg >= -absalpha)
		{
			massPart = std::exp(-0.5*arg*arg);
		}
		else
		{
			A = std::pow(n/absalpha,n) * std::exp(-0.5*absalpha*absalpha);
			B = n/absalpha - absalpha;
			massPart = A / std::pow(B - arg,n);
		}
	}
	else if(type == "combinatorial")
	{
		// Threshold function
		double A = KKpars[0].value;
		double B = KKpars[1].value;
		double C = KKpars[2].value;
		double M = KKpars[3].value;
		double arg = mKK - M;
		if(arg <= 0) return 0;
		double ratio = mKK/M;
		double val = (1- exp(-arg/C))* pow(ratio, A) + B*(ratio-1);
		massPart = val > 0 ? val : 0;
	}
	else if(type == "histogram")
		massPart = mKKhist.Eval({mKK});
	else if(!lineshape.empty())
	{
		massPart = std::norm(KKLineShape->massShape(mKK))*1e-4;
	}
	double angularPart = angulardistribution.Evaluate(datapoint);
	if(std::isnan(massPart)) std::cerr << "Mass part is nan" << std::endl;
	if(massPart < 0) std::cerr << "Mass part is negative" << std::endl;
	if(std::isnan(angularPart)) std::cerr << "Angular part is nan" << std::endl;
	if(angularPart < 0) std::cerr << "Angular part is negative" << std::endl;
	return massPart * angularPart;
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
	if(!lineshape.empty()) Bs2PhiKK::UpdateLineshape(lineshape, *KKLineShape, KKpars);
}
vector<ObservableRef> Bs2PhiKKBackgroundComponent::GetPhysicsParameters() const
{
	vector<ObservableRef> parameters;
	for(const auto& set: {{fraction},KKpars})
		for(const auto& par: set)
			parameters.push_back(par.name);
	return parameters;
}

