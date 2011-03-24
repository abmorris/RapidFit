/**
        @interface IMinimiser

        Interface for all function minimisers

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef I_MINIMISER_H
#define I_MINIMISER_H

#include "FitFunction.h"
#include "FitResult.h"

class IMinimiser
{
	public:
		virtual void Minimise( FitFunction* ) = 0;
		virtual FitResult * GetFitResult() = 0;
		virtual void ContourPlots( vector< pair< string, string > > ) = 0;
		virtual void SetOutputLevel( int ) = 0;
};

#endif
