#ifndef __LEGENDREMOMENTSHAPE_H__
#define __LEGENDREMOMENTSHAPE_H__
#include <vector>
#include <string>
#include "IDataSet.h"
#include "TFile.h"
#include "TTree.h"
using std::string;
using std::vector;
class LegendreMomentShape
{
	public:
		LegendreMomentShape(); // Declare without coefficients (e.g. if you want to call Generate() or Open() later)
		LegendreMomentShape(string); // Immediately call Open() on the passed string
		LegendreMomentShape(const LegendreMomentShape&);
		~LegendreMomentShape();
		void Open(string); // Load the coefficients from a file
		void Save(string); // Save generated coefficients to a file
		void Generate(IDataSet*, PhaseSpaceBoundary*, string, string, string, string); // strings are variable names: mass, phi, cosθ1, cosθ2
		void SetMax(double _l_max, double _i_max, double _k_max, double _j_max) // Only needed when generating coefficients; loaded from file otherwise
		{
			l_max = _l_max;
			i_max = _i_max;
			k_max = _k_max;
			j_max = _j_max;
		}
		static double Moment(int,int,int,int,double,double,double,double); // l, i, k, j, mass_mapped, phi, cosθ1, cosθ2
		double Evaluate(double,double,double,double); // mass, phi, cosθ1, cosθ2
		double mKK_min;
		double mKK_max;
	protected:
		struct coefficient
		{
			int l,i,j,k;
			double val;
			void print()
			{
				printf("c[%d][%d][%d][%d] = %f\n", l, i, k, j, val);
			}
		};
		vector<coefficient> coeffs;
		bool init;
	private:
		double**** newcoefficients();
		void deletecoefficients(double****);
		void storecoefficients(double****);
		void printcoefficients();
		int l_max;
		int i_max;
		int k_max;
		int j_max;
		bool copied;
};
#endif
