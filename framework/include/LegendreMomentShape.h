#ifndef __LEGENDREMOMENTSHAPE_H__
#define __LEGENDREMOMENTSHAPE_H__
#include <vector>
#include <string>
#include <array>
#include <random>
#include "IDataSet.h"
#include "TFile.h"
#include "TTree.h"
class LegendreMomentShape
{
	public:
		LegendreMomentShape(); // Declare without coefficients (e.g. if you want to call Generate() or Open() later)
		LegendreMomentShape(std::string); // Immediately call Open() on the passed string
		LegendreMomentShape(const LegendreMomentShape&);
		~LegendreMomentShape();
		void Open(const std::string, const int entry = 0); // Load the coefficients from a file
		void Save(const std::string); // Save generated coefficients to a file
		void Generate(std::vector<DataPoint>&, const PhaseSpaceBoundary*, const std::string, const std::string, const std::string, const std::string, const bool); // strings are variable names: mass, phi, cosθ1, cosθ2. boolean is "mass_dependent": don't divide by phase space function if true
		void SetMax(const double _l_max, const double _i_max, const double _k_max, const double _j_max) // Only needed when generating coefficients; loaded from file otherwise
		{
			l_max = _l_max;
			i_max = _i_max;
			k_max = _k_max;
			j_max = _j_max;
		}
		static double Moment(const int,const int,const int,const int,const double,const double,const double,const double); // l, i, k, j, mass_mapped, phi, cosθ1, cosθ2
		double Evaluate(const std::array<double,4>&) const;
		double Evaluate(const double,const double,const double,const double) const; // mass, phi, cosθ1, cosθ2
		double mKK_min;
		double mKK_max;
	private:
		struct coefficient
		{
			coefficient(int _l, int _i, int _k, int _j, double _val, double _err)
				: l(_l), i(_i), k(_k), j(_j), val(_val), err(_err)
				, randval(std::normal_distribution<double>(val, err)) {}
			coefficient(const coefficient& other)
				: l(other.l), i(other.i), k(other.k), j(other.j), val(other.val), err(other.err)
				, randval(std::normal_distribution<double>(val, err)) {}
			int l,i,j,k;
			double val;
			double err;
			std::normal_distribution<double> randval;
			void print() const
			{
				printf("c[%d][%d][%d][%d] = %f ± %f with significance %f\n", l, i, k, j, val, err, std::abs(val/err));
			}
		};
		std::vector<coefficient> coeffs;
		bool init;
		void printcoefficients() const;
		int l_max;
		int i_max;
		int k_max;
		int j_max;
		bool copied;
};
#endif
