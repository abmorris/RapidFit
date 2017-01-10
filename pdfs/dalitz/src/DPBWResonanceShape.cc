#include "DPBWResonanceShape.hh"
#include "DPHelpers.hh"
#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mRR, double gammaRR, int L, double mm1, double mm2, double RR):
	 mR(mRR)
	,gammaR(gammaRR)
	,LR(L)
	,m1(mm1)
	,m2(mm2)
{
	barrier = DPBarrierFactor(LR,RR);
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
}

DPBWResonanceShape::DPBWResonanceShape( const DPBWResonanceShape& other ) : DPMassShape(other)
	,mR(other.mR)
	,gammaR(other.gammaR)
	,LR(other.LR)
	,m1(other.m1)
	,m2(other.m2)
	,barrier(other.barrier)
{
}

std::complex<double> DPBWResonanceShape::massShape(const double m) const
{
	std::complex<double> result(1,0);
	std::complex<double> denominator(mR*mR-m*m,-mR*gamma(m));
	result/=denominator;
	return result;
}

double DPBWResonanceShape::gamma(const double m) const
{
	double pp=DPHelpers::daughterMomentum(m,m1,m2);;  // momentum of daughter at the actual mass
	double bb=barrier.barrier(pR0,pp);  // Barrier factor
	double gg=gammaR*mR/m*bb*bb*std::pow(pp/pR0,2*LR+1);
	return gg;
}

void DPBWResonanceShape::setParameters(const std::vector<double>& pars)
{
	setResonanceParameters(pars[0],pars[1]);
}

void DPBWResonanceShape::setResonanceParameters(const double mass, const double sigma )
{
	mR = mass;
	gammaR = sigma;
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
}

