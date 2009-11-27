/**
        @class StatisticsFunctions

        A collection of static methods for statistics calculations

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef STATISTICS_FUNCTIONS_H
#define STATISTICS_FUNCTIONS_H

#include <string>
#include <vector>
#include "PhaseSpaceBoundary.h"
#include "IPDF.h"
#include "IDataSet.h"

using namespace std;

class StatisticsFunctions
{
	public:
		static double Mean( vector<double> );
		static double Variance( vector<double> );
		static int OptimumBinNumber( vector<double> );
		static double Maximum( vector<double> );
		static double Minimum( vector<double> );
		static vector< vector<double> > DiscreteCombinations( vector<string>*, PhaseSpaceBoundary*, vector<string>&, vector<string>&, vector< vector<double> >& );
		static void DoDontIntegrateLists( IPDF*, PhaseSpaceBoundary*, vector<string>*, vector<string>&, vector<string>& );
		static vector<DataPoint> DataAverage( IDataSet*, vector< vector<double> >, vector< vector<double> >, vector<string>, vector<string>, vector<string>&, vector<double>& );
};

#endif
