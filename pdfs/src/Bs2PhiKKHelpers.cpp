#include "Bs2PhiKKHelpers.h"

double Bs2PhiKK::mBs  = 5.36677;
double Bs2PhiKK::mK   = 0.493677;
double Bs2PhiKK::mpi  = 0.139570;

void Bs2PhiKK::PhysPar::Update(const ParameterSet* pars)
{
	value = pars->GetPhysicsParameter(name)->GetValue();
	if(std::isnan(value))
		std::cerr << name.Name() << " has been given a nan value!" << std::endl;
}
