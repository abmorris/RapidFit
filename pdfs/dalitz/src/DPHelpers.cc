#include "DPHelpers.hh"
#include <cmath>

double DPHelpers::daughterMomentum(double m, double m1, double m2)
{
	double momentum;
	momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
	momentum=sqrt(momentum);
	momentum/=2*m;
	return momentum;
}

