#ifndef AccCorr_v3_H
#define AccCorr_v3_H

#include <iostream>
#include "PDFConfigurator.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH3F.h"
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Mathematics.h"
#include "XMLConfigReader.h"
#include <vector>
#include "Bs2PhiPhi_AccCorr.h"
#include <string>

class AccCorr_v3
{

	public:
		AccCorr_v3();
		AccCorr_v3( const AccCorr_v3& );
		virtual void AccCorr_v3Init(PDFConfigurator*);
		virtual ~AccCorr_v3();
		virtual double AccEval(vector<double>);
		void MakeWeights(PDFConfigurator*);
		vector<double> ReturnWeights();
		void WeightsFromFile(string);	
		double CleverInterpolate(TH3F*,double,double,double);
		double histoArea(TH3F*);
		double histoArea(TH2F*);
		double getHelicityWeight(PDFConfigurator*,IPDF*,IPDF*,DataPoint*);
		double getBWeight(string,double,double,double,double);
		TH2F makeRWHistKK(string,string,PDFConfigurator*,IDataSet*,IPDF*,IDataSet*,IPDF*);

	private:
		AccCorr_v3& operator = ( const AccCorr_v3& );		
		bool useAccFit;
		int AccLoc;
		string Dim;
		TH3F* histo;
		//TH1F* histo1D;
                //TAxis* xaxis;
                //TAxis* yaxis;
                //TAxis* zaxis;
		int nxbins, nybins, nzbins;
                double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz, total_num_entries;
		double angNorm;
		double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15;

};

#endif
