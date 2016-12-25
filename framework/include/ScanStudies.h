
#pragma once
#ifndef ScanStudies_H
#define ScanStudies_H

//	RapidFit Headers
#include "FitResultVector.h"
#include "MinimiserConfiguration.h"
#include "ConstraintFunction.h"
#include "OutputConfiguration.h"
#include "ParameterSet.h"
#include "PDFWithData.h"
#include "FitFunctionConfiguration.h"
#include "ResultFormatter.h"
#include "DebugClass.h"
//	System Headers
#include <vector>
#include <string>

using std::vector;
using std::string;

class ScanStudies
{

	public:

		//	New Interface to Scanning Code

		/*!
		 * @brief Function to call to perform a ContourScan
		 *
		 * @Param inputMinimiser
		 *
		 * @Param inputParam
		 *
		 * @Param inputPDFWithData
		 *
		 * @Param inputConstraints
		 *
		 * @Param inputConfig
		 *
		 * @Param param1
		 *
		 * @Param param2
		 *
		 * @Param output
		 *
		 * @Param forceContinue
		 */
		static vector<FitResultVector*> ContourScan( MinimiserConfiguration* inputMinimiser, FitFunctionConfiguration* inputFunction, ParameterSet* inputParam,
				const vector< PDFWithData* > inputPDFWithData, const vector< ConstraintFunction* > inputConstraints, OutputConfiguration* inputConfig,
				const string param1, const string param2, const int output=-999, bool forceContinue=false );

		/*!
		 * @brief Function to call to perform a SingleScan
		 *
		 * @Param inputMinimiser
		 *
		 * @Param inputParam
		 *
		 * @Param inputPDFWithData
		 *
		 * @Param inputConstraints
		 *
		 * @Param inputConfig
		 *
		 * @Param param
		 *
		 * @Param output
		 *
		 * @Param forceContinue
		 */
		static FitResultVector* SingleScan(  MinimiserConfiguration* inputMinimiser, FitFunctionConfiguration* inputFunction, ParameterSet* inputParam,
				const vector< PDFWithData* > inputPDFWithData, const vector< ConstraintFunction* > inputConstraints, OutputConfiguration* inputConfig,
				const string param, const int output=-999, bool forceContinue=false );

	private:

		static void DoScan( MinimiserConfiguration *, FitFunctionConfiguration *, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, ScanParam*, FitResultVector*, const int, bool forceContinue=false );

		static void DoScan2D( MinimiserConfiguration*, FitFunctionConfiguration*, ParameterSet*, const vector< PDFWithData* >,
				const vector< ConstraintFunction* >, const pair<ScanParam*, ScanParam* >, vector<FitResultVector*>*, const int, bool forceContinue=false );

};

#endif

