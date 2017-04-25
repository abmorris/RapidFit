#ifndef DP_Spline_SHAPE
#define DP_Spline_SHAPE
#include "TSpline.h"
#include "DPMassShape.hh"
#include <complex>
class DPComplexSpline: public virtual DPMassShape
{
	public:
		DPComplexSpline(std::vector<double> breakpoints);
		DPComplexSpline(const DPComplexSpline&);
		~DPComplexSpline() {}
		typedef std::pair<double, std::complex<double>> complex_knot;
		std::complex<double> massShape(const double) const;
		void setParameters(const std::vector<double>&); // alternating triplets of masses, magnitudes and phases
		void setParameters(std::vector<complex_knot>);
	private:
		static bool compare_knot(const DPComplexSpline::complex_knot&, const DPComplexSpline::complex_knot&);
		TSpline3 ReSpline;
		TSpline3 ImSpline;
};
#endif

