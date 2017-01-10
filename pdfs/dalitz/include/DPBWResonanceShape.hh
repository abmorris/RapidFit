#ifndef DP_BW_RESONANCE_SHAPE
#define DP_BW_RESONANCE_SHAPE
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include <complex>
class DPBWResonanceShape: public virtual DPMassShape
{
	public:
		DPBWResonanceShape(double mR, double gammaR, int L, double m1, double m2, double R);
		DPBWResonanceShape(const DPBWResonanceShape&);
		~DPBWResonanceShape() {}
		std::complex<double> massShape(const double m) const;
		void setParameters(const std::vector<double>& pars);
	private:
		void Init();
		double mR;
		double gammaR;
		int LR;
		double m1;
		double m2;
		DPBarrierFactor barrier;
		double pR0;  // Momentum of daughters at mR
		double gamma(const double m) const;
		void setResonanceParameters(const double mass, const double sigma );
};
#endif
