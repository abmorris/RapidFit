#ifndef DP_Spline_SHAPE
#define DP_Spline_SHAPE
#include "TSpline.h"
#include "DPMassShape.hh"
#include <complex>
class DPComplexSpline: public virtual DPMassShape
{
	public:
		typedef std::pair<double, std::complex<double>> complex_knot;
		DPComplexSpline(std::vector<double> breakpoints);
		DPComplexSpline(std::vector<complex_knot> knots);
		DPComplexSpline(const DPComplexSpline&);
		~DPComplexSpline() {}
		std::complex<double> massShape(const double) const;
		void setParameters(const std::vector<double>&); // alternating triplets of masses, magnitudes and phases
		void setParameters(std::vector<complex_knot>);
		static bool compare_knot(const DPComplexSpline::complex_knot&, const DPComplexSpline::complex_knot&);
	private:
		TSpline3 ReSpline;
		TSpline3 ImSpline;
};
#endif

