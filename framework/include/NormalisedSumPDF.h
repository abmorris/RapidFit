/**
  @class NormalisedSumPDF

  An implementation of IPDF for adding the values of two other IPDFs, normalising them relative to each other.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-11-12
  */

#pragma once
#ifndef NORMALISED_SUM_PDF_H
#define NORMALISED_SUM_PDF_H

//	RapidFit Headers
#include "IPDF.h"
#include "BasePDF.h"
#include "RapidFitIntegrator.h"
#include "ComponentRef.h"
//	System Headers
#include <array>
#include <vector>
#include <string>
#include <map>

using std::vector;
using std::string;
using std::map;
using std::array;

class NormalisedSumPDF : public BasePDF
{
	public:
		NormalisedSumPDF( PDFConfigurator* config );
		//NormalisedSumPDF( IPDF*, IPDF*, PhaseSpaceBoundary*, string="fractionName" );
		NormalisedSumPDF( const NormalisedSumPDF& );
		~NormalisedSumPDF();

		void SetComponentStatus( const bool input );

		bool GetComponentStatus() const;

		vector<string> PDFComponents();

		void TurnCachingOff();

		virtual void SetUpIntegrator( const RapidFitIntegratorConfig* thisConfig );

		//Return the integral of the function over the given boundary
		double Normalisation( DataPoint*, PhaseSpaceBoundary* );

		//Return the function value at the given point
		double Evaluate( DataPoint* );
		double EvaluateForNumericIntegral( DataPoint* );

		//Set the function parameters
		bool SetPhysicsParameters( ParameterSet* );

		//Return a prototype data point
		vector<string> GetPrototypeDataPoint();

		//Return a prototype set of physics parameters
		vector<string> GetPrototypeParameterSet();

		//Return a list of parameters not to be integrated
		vector<string> GetDoNotIntegrateList();

		bool GetNumericalNormalisation() const;

		//      Return components, component 0 by default
		double EvaluateComponent( DataPoint*, ComponentRef* = NULL );

		IPDF* GetFirstPDF() const;

		IPDF* GetSecondPDF() const;

		string GetFractionName() const;

		bool GetCachingEnabled() const;

		void SetCachingEnabled( bool Input );

		array<double,2> GetCachedIntegrals(DataPoint*, PhaseSpaceBoundary* ) const;

		void SetCachedIntegrals(array<double,2>, DataPoint*, PhaseSpaceBoundary* );

		string XML() const;

		void SetDebugMutex( pthread_mutex_t* Input, bool =true );

		virtual string GetComponentName( ComponentRef* = NULL );

		void ChangePhaseSpace( PhaseSpaceBoundary * InputBoundary );

		PhaseSpaceBoundary* GetPhaseSpace() const;

		vector<IPDF*> GetChildren() const;

	private:
		//	Uncopyable!
		NormalisedSumPDF& operator=(const NormalisedSumPDF&);
		void MakePrototypes( PhaseSpaceBoundary* );

		double GetFirstIntegral( DataPoint* ) const;
		double GetSecondIntegral( DataPoint* ) const;
		map<unsigned int, double> firstIntegral, secondIntegral;

		vector<string> prototypeDataPoint, prototypeParameterSet, doNotIntegrateList;
		IPDF * firstPDF;
		IPDF * secondPDF;
		double firstFraction, firstIntegralCorrection, secondIntegralCorrection;
		ObservableRef fractionName;
		PhaseSpaceBoundary * integrationBoundary;

		bool _plotComponents;
		bool inEvaluate = false;
};

#endif

