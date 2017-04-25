#include "DPComplexSpline.hh"
#include <iostream>

// Constructor creates the two spline objects and sets the y-values to zero
DPComplexSpline::DPComplexSpline(std::vector<double>& breakpoints) : DPMassShape()
{
	std::vector<double> dummy_yvalues(breakpoints.size(),0);
	ReSpline = TSpline3("real",breakpoints.data(),dummy_yvalues.data(),breakpoints.size());
	ImSpline = TSpline3("imag",breakpoints.data(),dummy_yvalues.data(),breakpoints.size());
}
// Copy constructor
DPComplexSpline::DPComplexSpline(const DPComplexSpline& other) : DPMassShape(other)
	, ReSpline(other.ReSpline)
	, ImSpline(other.ImSpline)
{
}
// 
std::complex<double> DPComplexSpline::massShape(const double m) const
{
	double re = ReSpline.Eval(m);
	double im = ImSpline.Eval(m);
	return std::complex<double>(re,im);
}
void DPComplexSpline::setParameters(const std::vector<double>& pars)
{
	if(pars.size() % 3 != 0 )
		throw std::out_of_range("DPComplexSpline ERROR: number of real parameters (" + std::to_string(pars.size()) + ") is not a multiple of 3");
	std::vector<complex_knot> complex_params;
	for(int i = 0; i < (int)pars.size(); i+=3)
		complex_params.push_back({pars[i],std::polar<double>(pars[i+1],pars[i+2])});
	setParameters(complex_params);
}
void DPComplexSpline::setParameters(const std::vector<complex_knot>& pars)
{
	int n = check(pars);
	for(int i = 0; i < n; i++)
	{
		ReSpline.SetPoint(i, pars[i].first, pars[i].second.real());
		ImSpline.SetPoint(i, pars[i].first, pars[i].second.imag());
	}
}
int DPComplexSpline::check(const std::vector<complex_knot>& pars)
{
	const int n = pars.size();
	const int npr = ReSpline.GetNp();
	const int npi = ImSpline.GetNp();
	if(n != npr || n != npi)
		throw std::out_of_range("DPComplexSpline ERROR: " + std::to_string(n) + " pars, "
		                                           + std::to_string(npr) + " real knots, "
		                                           + std::to_string(npi) + " imaginary knots.");
	else
		return n;
}

