#ifndef __LEGENDREMOMENTSHAPE_H__
#define __LEGENDREMOMENTSHAPE_H__
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
using std::string;
using std::vector;
class LegendreMomentShape
{
	public:
		LegendreMomentShape(string);
		LegendreMomentShape(const LegendreMomentShape&);
		~LegendreMomentShape();
		double Evaluate(double,double,double,double);
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
		string filename;
		double**** c;
		void createcoefficients();
		void storecoefficients();
		void printcoefficients();
		int l_max;
		int i_max;
		int k_max;
		int j_max;
		TTree* tree;
		TFile* file;
		bool copied;
};
#endif
