#ifndef DP_LASS_SHAPE
#define DP_LASS_SHAPE

#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include <complex>

class DPLassShape: public virtual DPMassShape
{
	public:
		DPLassShape(double mR, double gammaR, int L, double m1,double m2, double R, double a, double r);
		DPLassShape( const DPLassShape& );
		~DPLassShape();
		std::complex<double> massShape(double m);
		void setResonanceParameters(double a, double r);
		void setParameters(double* pars);
	private:
		void Init();
		double mR;
		double gammaR;
		int LR;
		double m1;
		double m2;
		double a;
		double r;
		double R; // Blatt-Weisskopf radius
		DPBarrierFactor* barrier;
		double pR0;  // Momentum of daughters at mR
		double gamma(double m);
};

#endif

