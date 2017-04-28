#include "DPComplexSpline.hh"
#include "TGraph.h"
#include "TGraphPolar.h"
#include "TCanvas.h"
#include <algorithm>
#include <iostream>
int main()
{
	std::vector<DPComplexSpline::complex_knot> knots =
	{
		{0.98, std::polar<double>(0.22,0.52)},
		{1.00, std::polar<double>(0.15,0.50)},
		{1.20, std::polar<double>(0.66,0.34)},
		{1.40, std::polar<double>(0.43,0.30)},
		{1.60, std::polar<double>(0.8,0.42)},
		{1.80, std::polar<double>(1.0,0.50)}
	};
	DPComplexSpline SplineShape(knots);
	std::sort(knots.begin(), knots.end(), DPComplexSpline::compare_knot);
	std::vector<double> k_mass, k_radius, k_theta;
	for(const auto& knot: knots)
	{
		k_mass.push_back(knot.first);
		k_radius.push_back(std::abs(knot.second));
		k_theta.push_back(std::arg(knot.second));
	}
	std::vector<double> mass, radius, theta;
	for(double x = knots.begin()->first; x < (knots.end()-1)->first; x+=0.01)
	{
		std::complex<double> y = SplineShape.massShape(x);
		mass.push_back(x);
		radius.push_back(std::abs(y));
		theta.push_back(std::arg(y));
	}
	TCanvas can("canvas","",500,500);
	// Argand diagram
	TGraphPolar k_polargraph(k_mass.size(), k_radius.data(), k_theta.data());
	k_polargraph.SetTitle("");
	k_polargraph.SetMarkerStyle(8);
	k_polargraph.SetMarkerSize(0.5);
	k_polargraph.Draw("AP");
	TGraphPolar polargraph(mass.size(), radius.data(), theta.data());
	polargraph.SetTitle("");
	polargraph.SetLineColor(kRed);
	polargraph.SetLineWidth(2);
	polargraph.Draw("Lsame");
	can.SaveAs("argand.pdf");
	// Magnitude plot
	TGraph k_absgraph(k_mass.size(), k_mass.data(), k_radius.data());
	k_absgraph.SetTitle("");
	k_absgraph.SetMarkerStyle(8);
	k_absgraph.SetMarkerSize(0.5);
	k_absgraph.Draw("AP");
	TGraph absgraph(mass.size(), mass.data(), radius.data());
	absgraph.SetLineColor(kRed);
	absgraph.SetLineWidth(2);
	absgraph.Draw("Lsame");
	can.SaveAs("abs.pdf");
	// Phase plot
	TGraph k_arggraph(k_mass.size(), k_mass.data(), k_theta.data());
	k_arggraph.SetTitle("");
	k_arggraph.SetMarkerStyle(8);
	k_arggraph.SetMarkerSize(0.5);
	k_arggraph.Draw("AP");
	TGraph arggraph(mass.size(), mass.data(), theta.data());
	arggraph.SetLineColor(kRed);
	arggraph.SetLineWidth(2);
	arggraph.Draw("Lsame");
	can.SaveAs("arg.pdf");
	return 0;
}
