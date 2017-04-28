#include "DPComplexSpline.hh"
#include <algorithm>

// Constructor creates the two spline objects and sets the y-values to zero
DPComplexSpline::DPComplexSpline(std::vector<double> breakpoints) : DPMassShape()
{
	std::vector<double> dummy_yvalues(breakpoints.size(),0);
	std::sort(breakpoints.begin(),breakpoints.end());
	ReSpline = TSpline3("real",breakpoints.data(),dummy_yvalues.data(),breakpoints.size());
	ImSpline = TSpline3("imag",breakpoints.data(),dummy_yvalues.data(),breakpoints.size());
}
// Copy constructor
DPComplexSpline::DPComplexSpline(const DPComplexSpline& other) : DPMassShape(other)
	, ReSpline(other.ReSpline)
	, ImSpline(other.ImSpline)
{
}
// Return the amplitude as calculated by the splines
std::complex<double> DPComplexSpline::massShape(const double m) const
{
	double re = ReSpline.Eval(m);
	double im = ImSpline.Eval(m);
	return std::polar<double>(re,im);
}
// Construct the complex knots {mass, magnitude*exp(i*phase)} from triplets of real numbers
void DPComplexSpline::setParameters(const std::vector<double>& pars)
{
	if(pars.size() % 3 != 0 )
		throw std::out_of_range("DPComplexSpline ERROR: number of real parameters (" + std::to_string(pars.size()) + ") is not a multiple of 3");
	std::vector<complex_knot> complex_params;
	for(int i = 0; i < (int)pars.size(); i+=3)
		complex_params.push_back({pars[i],std::polar<double>(pars[i+1],pars[i+2])});
	setParameters(complex_params);
}
// Set the positions and values of the complex knots
// If you want to merge this function with the other setParameters function, be careful to reproduce the conversion from std::polar to std::complex
void DPComplexSpline::setParameters(std::vector<complex_knot> pars)
{
	const int n = pars.size();
	// Make sure the vector of complex knot values has the same size as there are knots in the real and imaginary splines
	const int npr = ReSpline.GetNp(), npi = ImSpline.GetNp();
	if(n != npr || n != npi)
		throw std::out_of_range("DPComplexSpline ERROR: " + std::to_string(n) + " pars, "
		                                                  + std::to_string(npr) + " real knots, "
		                                                  + std::to_string(npi) + " imaginary knots.");
	// Sort the knots from lowest to highest mass
	std::sort(pars.begin(), pars.end(), compare_knot);
	// Set the knot values in the component splines
	for(int i = 0; i < n; i++)
	{
		ReSpline.SetPoint(i, pars[i].first, pars[i].second.abs());
		ImSpline.SetPoint(i, pars[i].first, pars[i].second.arg());
	}
}
// Comparison function to aid in sorting knots from lowest to highest mass
bool DPComplexSpline::compare_knot(const DPComplexSpline::complex_knot& left, const DPComplexSpline::complex_knot& right)
{
	return left.first < right.first;
}

