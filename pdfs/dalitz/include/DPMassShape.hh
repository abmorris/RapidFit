#ifndef DP_MASS_SHAPE
#define DP_MASS_SHAPE

#include <complex>
#include <iostream>

class DPMassShape
{
	public:
		DPMassShape() {};
		DPMassShape( const DPMassShape& ) {};
		virtual ~DPMassShape() {};
		virtual std::complex<double> massShape(double m) = 0;
		virtual void setResonanceParameters(double mass, double sigma) = 0;
		virtual void setParameters(double* pars) = 0;
};

#endif
