#ifndef DP_NONRESONANT_SHAPE
#define DP_NONRESONANT_SHAPE

#include "DPMassShape.hh"
#include <complex>

class DPNonresonant: public virtual DPMassShape
{
	public:
		DPNonresonant();
		DPNonresonant(const DPNonresonant& other);
		~DPNonresonant();
		std::complex<double> massShape(double m);
		void setParameters(double* pars) { (void)pars; }
	private:
		void setResonanceParameters( double mass, double sigma ) { (void)mass; (void)sigma; }
};

#endif
