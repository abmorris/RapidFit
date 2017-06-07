#include "DPBWResonanceShape.hh"
#include "DPHelpers.hh"
#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mRR, double gammaRR, int L, double mm1, double mm2, double RR) : DPMassShape()
	,mR(mRR)
	,gammaR(gammaRR)
	,LR(L)
	,m1(mm1)
	,m2(mm2)
	,pR0(DPHelpers::daughterMomentum(mR,m1,m2))
	,barrier(DPBarrierFactor(LR,RR,pR0))
{
}

DPBWResonanceShape::DPBWResonanceShape( const DPBWResonanceShape& other ) : DPMassShape(other)
	,mR(other.mR)
	,gammaR(other.gammaR)
	,LR(other.LR)
	,m1(other.m1)
	,m2(other.m2)
	,pR0(other.pR0)
	,barrier(other.barrier)
{
}

std::complex<double> DPBWResonanceShape::massShape(const double m) const
{
	double width = gamma(m);
	std::complex<double> result(1, 0);
	const double re = mR*mR - m*m;
	const double im = -mR * width;
	const std::complex<double> denominator(re, im);
	result /= denominator;
	return result;
}

double DPBWResonanceShape::gamma(const double m) const
{
	double pR = DPHelpers::daughterMomentum(m,m1,m2);;  // momentum of daughter at the actual mass
	double BF_sq = barrier.barrier_sq(pR);  // Barrier factor squared
	return gammaR * (mR/m) * BF_sq * std::pow(pR/pR0, 2*LR + 1);
}

void DPBWResonanceShape::setParameters(const std::vector<double>& pars)
{
	setResonanceParameters(pars[0], pars[1]);
	pR0 = DPHelpers::daughterMomentum(mR,m1,m2);
	barrier.setparameters(pars[2], pR0);
}

void DPBWResonanceShape::setResonanceParameters(const double mass, const double sigma )
{
	mR = mass;
	gammaR = sigma;
}

