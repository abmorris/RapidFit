#include "DPBWResonanceShape.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"
#include "DPBarrierL4.hh"
#include "DPBarrierL5.hh"
#include "DPHelpers.hh"
#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mRR, double gammaRR, int L, double mm1, double mm2, double RR):
	 mR(mRR)
	,gammaR(gammaRR)
	,LR(L)
	,m1(mm1)
	,m2(mm2)
	,R(RR)
{
	Init();
}

DPBWResonanceShape::DPBWResonanceShape( const DPBWResonanceShape& other ) : DPMassShape(other)
	,mR(other.mR)
	,gammaR(other.gammaR)
	,LR(other.LR)
	,m1(other.m1)
	,m2(other.m2)
	,R(other.R)
{
	Init();
}

void DPBWResonanceShape::Init()
{
	switch (LR)
	{
		case 0: barrier=new DPBarrierL0(R);
		        break;
		case 1: barrier=new DPBarrierL1(R);
		        break;
		case 2: barrier=new DPBarrierL2(R);
		        break;
		case 3: barrier=new DPBarrierL3(R);
		        break;
		case 4: barrier=new DPBarrierL4(R);
		        break;
		case 5: barrier=new DPBarrierL5(R);
		        break;
		default: std::cerr<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
		         barrier=new DPBarrierL0(R);
		         break;
	}
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
}

DPBWResonanceShape::~DPBWResonanceShape()
{
	delete barrier;
}

std::complex<double> DPBWResonanceShape::massShape(double m)
{
	std::complex<double> result(1,0);
	std::complex<double> denominator(mR*mR-m*m,-mR*gamma(m));
	result/=denominator;
	return result;
}

double DPBWResonanceShape::gamma(double m)
{
	double pp=DPHelpers::daughterMomentum(m,m1,m2);;  // momentum of daughter at the actual mass
	double bb=barrier->barrier(pR0,pp);  // Barrier factor
	double gg=gammaR*mR/m*bb*bb*std::pow(pp/pR0,2*LR+1);
	return gg;
}

void DPBWResonanceShape::setParameters(double* pars)
{
	setResonanceParameters(pars[0],pars[1]);
}

void DPBWResonanceShape::setResonanceParameters( double mass, double sigma )
{
	mR = mass;
	gammaR = sigma;
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
}

