#include "DPBarrierFactor.hh"
#include <iostream>
#include <cmath>
DPBarrierFactor::DPBarrierFactor() : radius(3.1)
{
	function = &DPBarrierFactor::FunctionL0;
}
DPBarrierFactor::DPBarrierFactor(unsigned spin, double _radius) : radius(_radius)
{
	if(std::isnan(radius)) std::cerr << "\t\tDPBarrier radius is nan " << std::endl;
	switch(spin)
	{
		case 0:
			function = &DPBarrierFactor::FunctionL0;
			break;
		case 1:
			function = &DPBarrierFactor::FunctionL1;
			break;
		case 2:
			function = &DPBarrierFactor::FunctionL2;
			break;
		case 3:
			function = &DPBarrierFactor::FunctionL3;
			break;
		default:
			std::cerr << "Only up to spin-3 barrier factors are currently implemented. Defaulting to spin-0." << std::endl;
			function = &DPBarrierFactor::FunctionL0;
	}
}
// Static functions for different spins
double DPBarrierFactor::FunctionL0(const double z)
{
	(void)z;
	return 1;
}
double DPBarrierFactor::FunctionL1(const double z)
{
	return std::sqrt(1+z);
}
double DPBarrierFactor::FunctionL2(const double z)
{
	return std::sqrt(z*z+3*z+9);
}
double DPBarrierFactor::FunctionL3(const double z)
{
	return std::sqrt(z*z*z+6*z*z+45*z+225);
}
// Evaluate the barrier factor
double DPBarrierFactor::barrier(const double p0, const double p) const
{
	double z  = p*p*radius*radius;
	double z0 = p0*p0*radius*radius;
	double factor = function(z0)/function(z);
	if(std::isnan(factor)) std::cerr << "\t\t\tB(" << p << ", " << p0 << " | " << radius << ") = nan" << std::endl;
	return factor;
}

