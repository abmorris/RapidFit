/** @class Bs2PhiPhi_v4 Bs2PhiPhi_v4.h
 *
 *  RapidFit PDF for Bs2PhiPhi_v4
 *
 *  @author Sean Benson
 *  @date 11 Oct 2011
 */

#ifndef Bs2PhiPhi_v4_H
#define Bs2PhiPhi_v4_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif
#include "SlicedAcceptance.h"
#include "TMath.h"
#include "AccCorr_v3.h"
#include "AccCorr_v3.h"
#include "AngularAcceptance.h"
#include "MistagCalib3fb.h"

class Bs2PhiPhi_v4 : public BasePDF
{
    public:
        Bs2PhiPhi_v4(PDFConfigurator*);
        Bs2PhiPhi_v4(const Bs2PhiPhi_v4&);
	~Bs2PhiPhi_v4();

        //Calculate the PDF value
        virtual double Evaluate(DataPoint*);
        double EvaluateBase() const ;
        virtual double EvaluateTimeOnly(DataPoint*);
        double EvaluateBaseTimeOnly() const ;
        virtual bool SetPhysicsParameters(ParameterSet*);
        virtual vector<string> GetDoNotIntegrateList();
	void AssignTAI(PhaseSpaceBoundary *);
	//double EvaluateForNumericGeneration( DataPoint* );
	vector<string> PDFComponents();
	double EvaluateComponent( DataPoint*, ComponentRef*);

    protected:
	int componentIndex;
	bool _usePlotComponents;
	bool performingComponentProjection;
        //Calculate the PDF normalisation
        virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);
        double NormalisationBase(PhaseSpaceBoundary*);
	void preCalculateTimeFactors() const ;
	void preCalculateTimeIntegrals() const ;
	void prepareCDS();
        double gamma_s, deltaGamma, gamma_l, gamma_h, deltaM, Phi_s, Azero_sq, Aperp_sq, Apara_sq, delta_1, delta_2, lambda;
	double _CC, _DD, _SS;
	AccCorr_v3* accep1;
	AngularAcceptance* angAccPETE;
	// Mistag class
	MistagCalib3fb* mCalib;
    private:
	Bs2PhiPhi_v4& operator=( const Bs2PhiPhi_v4& );
        void MakePrototypes();
	
	//Integraed TD values
	mutable double IntExpGammas, IntExpGammal, IntExpGammah;
	mutable double intcos,intsin,intexpsin,intexpcos;
	mutable double expSin, expCos, expL, expH;

	// sWave parameters
	double Ass_sq, As_sq, deltaS, deltaSS;
	double As_sq_in,Ass_sq_in;
	bool _useSWAVE;
	inline bool useSWAVE() const { return _useSWAVE ; }
        ObservableRef Ass_sqName;
        ObservableRef As_sqName;
        ObservableRef deltaSSName;
        ObservableRef deltaSName;
	// sWave coupling
        ObservableRef CspName;
	double Csp;

        //Cached values
        double cach_v1, cach_v2;
        double cach_Azero, cach_Apara, cach_Aperp;
        double cach_sinPhis, cach_cosPhis;
        double cach_sinDelta_[2];
        double cach_cosDelta_[2];
        double cach_cosDeltaSubtr;
        bool normalisationCacheValid, evaluationCacheValid;

        // These contain the ObservableRefs that correspond
        // to the physics parameter names that will be
        // used in the minimiser.
        ObservableRef gamma_sName;       // gamma_s
        ObservableRef deltaGammaName;
	ObservableRef deltaMName;        // delta mass
        ObservableRef Phi_sName;         // what we want to measure!
        ObservableRef Azero_sqName;      // amplitude
        ObservableRef Apara_sqName;      // amplitude
        ObservableRef Aperp_sqName;      // amplitude
        ObservableRef delta_1Name;       // strong phase
        ObservableRef delta_2Name;       // strong phase
        ObservableRef lambdaName;       // DCPV

        // These contain the ObservableRefs that correspond
        // to the observable names that are used in the
        // PDF. 
        ObservableRef timeName;    // proper time
        ObservableRef ctheta_1Name;  // angle between Phi1 p vector and its K+ daughter's p vector
        ObservableRef ctheta_2Name;  // angle between Phi2 p vector and its K+ daughter's p vector
        ObservableRef phiName;     // angle between the decay planes of the 2 Phi resonances
        /*
	ObservableRef tagName;     // B tag
        ObservableRef mistagName;  // B mistag
	ObservableRef mistagP0Name;
	ObservableRef mistagP1Name;
	ObservableRef mistagSetPointName;
	ObservableRef mistagDeltaP0Name;
	ObservableRef mistagDeltaP1Name;
	ObservableRef mistagDeltaSetPointName;
	*/
	// Time resolution
	ObservableRef resScaleName;			// Scale to multiply all Gaussians with 
	ObservableRef resOffsetName;			// time resolution offset 
	ObservableRef res1Name;				// time resolution narrow
        ObservableRef perEvTimeName;    // decay time error

	double resolutionScale ;
	double resolutionOffset ;
	double resolution1 ;
	double resolution ;
	double perEvTime ;
	
	// Observables as member variables
	double ctheta_1, ctheta_2, phi, time; 
	//double _mistag;
	//int tag;
	double tlo, thi ;
	double tlow, thigh;

	//double _mistagP0, _mistagP1, _mistagSetPoint;
	//double _mistagDeltaP0, _mistagDeltaP1, _mistagDeltaSetPoint;
	// Functions required for helicity*******************
	inline double cHt1sq() const { return (ctheta_1*ctheta_1) ; }
	inline double cHt2sq() const { return (ctheta_2*ctheta_2) ; }
	inline double sHt1sq() const { return (sin(acos(ctheta_1))*sin(acos(ctheta_1))) ; }
	inline double sHt2sq() const { return (sin(acos(ctheta_2))*sin(acos(ctheta_2))) ; }
	inline double c2Ht1() const { return (cos(2.*acos(ctheta_1))) ; }
	inline double c2Ht2() const { return (cos(2.*acos(ctheta_2))) ; }
	inline double s2Ht1() const { return (sin(2.*acos(ctheta_1))) ; }
	inline double s2Ht2() const { return (sin(2.*acos(ctheta_2))) ; }
	inline double cHphi() const { return (cos(phi)); }
	inline double sHphi() const { return (sin(phi)); }
	inline double c2Hphi() const { return (cos(2.*phi)); }
	inline double s2Hphi() const { return (sin(2.*phi)); }
	// **************************************************
	
	// TIME ACCEPTANCE *************************************
	inline bool useTimeAcceptance() const { return _useTimeAcceptance ; }			
	// TIME ERROR *************************************
	inline bool usePerEvTime() const { return _usePerEvTime ; }			
	// ACCEPTANCE WEIGHTS *************************************
	inline bool useAccepWeights() const { return _useAccepWeights ; }		
/*	
	// Using tagging information
	inline double q() const { return (double)tag ;}
	
	// Extra functions taken from the J/psi phi v6 PDF ********************************************************
		inline double mistagB() const {
			double returnValue = -1000.;

			if( (fabs(q()) < 0.5) || (fabs(q()) > 1.) ) {
				returnValue = 0.5;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0+(_mistagDeltaP0*0.5) + (_mistagP1+(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint+(_mistagDeltaSetPoint*0.5)) );
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagB() : _mistag < 0 so deltaMistag set to 0 also " << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagB() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0.5;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagB() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue;
		}

		inline double mistagBbar() const {
			double returnValue = -1000.;

			if( fabs(q()) < 0.5 ) {
				returnValue = 0.5 ;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0-(_mistagDeltaP0*0.5) + (_mistagP1-(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint-(_mistagDeltaSetPoint*0.5)) );
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagBbar() : _mistag < 0 so deltaMistag set to 0 also " << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagBbar() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0.5;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2PhiPhi_v4::mistagBbar() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue;
		}
		
		inline double D1() const {  return 1.0 - q()*(mistagB()-mistagBbar()); }
		inline double D2() const {  return q()*( 1.0 - mistagB() -mistagBbar() ); }
*/
		inline double D1() const {  return mCalib->D1(); }
		inline double D2() const {  return mCalib->D2(); }
		//.....................
		// C, D, S
		inline double cosphis() const { return _DD ; } //  _cosphis ; }
		inline double sinphis() const { return _SS ; } //  _sinphis ; }
		inline double CC() const { return _CC ; } //  _sinphis ; }
	// *******************************************************************************************************************
/*	
	inline double mistag( ) const {
		double returnValue;
		//Normal case
		returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
		//returnValue=_mistag;
		if( returnValue < 0. )  returnValue = 0. ;
		if( returnValue > 0.5) returnValue = 0.5 ; 
		return returnValue; }
*/
	inline double GF( ) const { return 9./32./TMath::Pi() ; }
	// Angle factor fot three angle PDFs in HELICITY basis ***
	// .................
	inline double HangleFactorA0A0( ) const { return 4. * cHt1sq() * cHt2sq() * GF(); }

	//..................
	inline double HangleFactorAPAP( ) const { return sHt1sq() * sHt2sq() * (1.+c2Hphi()) * GF(); }

	//..................
	inline double HangleFactorATAT( ) const { return sHt1sq() * sHt2sq() * (1.-c2Hphi()) * GF(); }

	//..................
	inline double HangleFactorImAPAT( ) const { return -2. * sHt1sq() * sHt2sq() * s2Hphi() * GF(); }

	//..................
	inline double HangleFactorReA0AP( )  const { return sqrt(2.) * s2Ht1() * s2Ht2() * cHphi() * GF(); }
		
	//..................
	inline double HangleFactorImA0AT( )  const { return -1. * sqrt(2.) * s2Ht1() * s2Ht2() * sHphi() * GF(); }
	// *******************************************************

	// sWave terms *******************************************
	//.................. f0-f0 (CP even)
	inline double f7( ) const { return 4./9. * GF(); }

	//.................. phi-f0 (CP odd)
	//inline double f8( ) const { return 4./3. *( (ctheta_1*ctheta_1) + (ctheta_2*ctheta_2) + 2.0*Csp*Csp*ctheta_1*ctheta_2 ) * GF() ; } // v1
	//inline double f8( ) const { return 4./3. * (ctheta_1 + ctheta_2) * (ctheta_1 + ctheta_2) * GF() ; } // Original
	inline double f8( ) const { return 4./3. * ( (ctheta_1 * ctheta_1) + (ctheta_2 * ctheta_2) ) * GF() ; }
	inline double f8a( ) const { return 4./3. * 2. * Csp*Csp * ctheta_1 * ctheta_2 * GF() ; }

	//.................. f0f0-phif0 interference
	inline double f9( ) const { return 8./(3.*sqrt(3.)) * (ctheta_1 + ctheta_2) * GF() ; }

	//.................. f0f0-phiphi interferences ...........
	// K10 - Re(A_0 A_SS*)
	inline double f10( ) const { return 8./3. * ctheta_1 * ctheta_2 * GF() ; }
	// K11 - Re(A_para A_SS*)
	inline double f11( ) const { return 4.*sqrt(2.)/3. * sin(acos(ctheta_1)) * sin(acos(ctheta_2)) * cHphi() * GF() ; }
	// K12 - Re(A_perp A_SS*)
	inline double f12( ) const { return -1.*4.*sqrt(2.)/3. * sin(acos(ctheta_1)) * sin(acos(ctheta_2)) * sHphi() * GF() ; }
	//........................................................

	//.................. phif0-phiphi interferences ...........
	// K13 - Re(A_0 A_S*)
	inline double f13( ) const { return 8./sqrt(3.) * ctheta_1 * ctheta_2 * (ctheta_1 + ctheta_2) * GF() ; }
	// K14 - Re(A_para A_S*)
	inline double f14( ) const { return 4.*sqrt(2.)/sqrt(3.) * sin(acos(ctheta_1)) * sin(acos(ctheta_2)) * (ctheta_1 + ctheta_2) * cHphi() * GF() ; }
	// K15 - Re(A_perp A_S*)
	inline double f15( ) const { return -1.*4.*sqrt(2.)/sqrt(3.) * sin(acos(ctheta_1)) * sin(acos(ctheta_2)) * (ctheta_1 + ctheta_2) * sHphi() * GF() ; }
	//........................................................
	//********************************************************
	
	// Time factors ******************************************
	// K1
	inline double timeA0A0( ) const { return 0.5 * Azero_sq * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K2
	inline double timeAPAP( ) const { return 0.5 * Apara_sq * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K3
	inline double timeATAT( ) const { return 0.5 * Aperp_sq * (
			
			D1() * ( (1.-cosphis())*expL + (1.+cosphis())*expH ) 
			+ D2() * ( -2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K4
	inline double timeAPAT( ) const { return sqrt(Apara_sq) * sqrt(Aperp_sq) * (
			
			D2() * ( sin(delta_1)*expCos - cos(delta_1)*expSin*cosphis() ) 
			+ D1() * ( -0.5*(expH-expL)*cos(delta_1)*sinphis()  + 0.5*(expH+expL)*sin(delta_1)*CC() )

			); }
	// K5
	inline double timeA0AP( ) const { return 0.5 * sqrt(Azero_sq) * sqrt(Apara_sq) * cos(delta_2-delta_1) * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K6
	inline double timeA0AT( ) const { return sqrt(Azero_sq) * sqrt(Aperp_sq) * (
			
			D2() * ( sin(delta_2)*expCos - cos(delta_2)*expSin*cosphis() ) 
			+ D1() * ( -0.5*(expH-expL)*cos(delta_2)*sinphis()  + 0.5*(expH+expL)*sin(delta_2)*CC() )

			); }

	// sWave time factors ************************************
	// K7
	inline double timeASSASS( ) const { return 0.5 * Ass_sq * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K8
	inline double timeASAS( ) const { return 0.5 * As_sq * (
			
			D1() * ( (1.-cosphis())*expL + (1.+cosphis())*expH ) 
			+ D2() * ( -2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K9
	inline double timeASASS( ) const { return Csp * sqrt(As_sq)*sqrt(Ass_sq) * (
			
			D2() * ( cos(deltaSS-deltaS)*expCos - sin(deltaSS-deltaS)*expSin*cosphis() ) 
			+ D1() * ( -0.5*(expH-expL)*sin(deltaSS-deltaS)*sinphis()  + 0.5*(expH+expL)*cos(deltaSS-deltaS)*CC() )

			); }
	// K10
	inline double timeA0ASS( ) const { return Csp * Csp * 0.5 * sqrt(Azero_sq) * sqrt(Ass_sq) * cos(deltaSS) * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K11
	inline double timeAPASS( ) const { return Csp * Csp * 0.5 * sqrt(Apara_sq) * sqrt(Ass_sq) * cos(delta_2-delta_1-deltaSS) * (
			
			D1() * ( (1.+cosphis())*expL + (1.-cosphis())*expH ) 
			+ D2() * ( 2.*expSin*sinphis() + 2.*expCos*CC() )

			); }
	// K12
	inline double timeATASS( ) const { return Csp * Csp * sqrt(Aperp_sq) * sqrt(Ass_sq) * (
			
			D2() * ( sin(delta_2-deltaSS)*expCos - cos(delta_2-deltaSS)*expSin*cosphis() ) 
			+ D1() * ( -0.5*(expH-expL)*cos(delta_2-deltaSS)*sinphis()  + 0.5*(expH+expL)*sin(delta_2-deltaSS)*CC() )
			
			); }
	// K13
	inline double timeA0AS( ) const { return Csp * sqrt(Azero_sq) * sqrt(As_sq) * (
			
			D2() * ( cos(deltaS)*expCos + sin(deltaS)*expSin*cosphis() ) 
			+ D1() * ( 0.5*(expH-expL)*sin(deltaS)*sinphis()  + 0.5*(expH+expL)*cos(deltaS)*CC() )

			); }
	// K14
	inline double timeAPAS( ) const { return Csp * sqrt(As_sq) * sqrt(Apara_sq) * (
			
			D2() * ( cos(delta_2-delta_1-deltaS)*expCos - sin(delta_2-delta_1-deltaS)*expSin*cosphis() ) 
			+ D1() * ( -0.5*(expH-expL)*sin(delta_2-delta_1-deltaS)*sinphis()  + 0.5*(expH+expL)*cos(delta_2-delta_1-deltaS)*CC() )

			); }
	// K15
	inline double timeATAS( ) const { return Csp * 0.5 * sqrt(Aperp_sq) * sqrt(As_sq) * sin(delta_2-deltaS) * (
			
			D1() * ( (1.-cosphis())*expL + (1.+cosphis())*expH ) 
			+ D2() * ( -2.*expSin*sinphis() + 2.*expCos*CC() )

			); }


	// Analytic integrals of time factors
	// K1
	inline double Int_timeA0A0( ) const { return 0.5 * Azero_sq * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

					); }
	// K2
	inline double Int_timeAPAP( ) const { return 0.5 * Apara_sq * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

					); }
	// K3
	inline double Int_timeATAT( ) const { return 0.5 * Aperp_sq * (
			
			D1() * ( (1.-cosphis())*IntExpGammal + (1.+cosphis())*IntExpGammah ) 
			+ D2() * ( -2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

					); }
	// K4
	inline double Int_timeAPAT( ) const { return sqrt(Apara_sq)*sqrt(Aperp_sq) * (
			
			D2() * ( sin(delta_1)*intexpcos - cos(delta_1)*intexpsin*cosphis() ) 
			+ D1() * ( -0.5*(IntExpGammah-IntExpGammal)*cos(delta_1)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*sin(delta_1)*CC() )

			); }
	// K5
	inline double Int_timeA0AP( ) const { return 0.5 * sqrt(Azero_sq) * sqrt(Apara_sq) * cos(delta_2-delta_1) * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

			); }
	// K6
	inline double Int_timeA0AT( ) const { return sqrt(Azero_sq)*sqrt(Aperp_sq) * (
			
			D2() * ( sin(delta_2)*intexpcos - cos(delta_2)*intexpsin*cosphis() ) 
			+ D1() * ( -0.5*(IntExpGammah-IntExpGammal)*cos(delta_2)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*sin(delta_2)*CC() )
			
			); }
	
	// sWave time integrals ****************************
	// K7
	inline double Int_timeASSASS( ) const { return 0.5 * Ass_sq * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

			); }
	// K8
	inline double Int_timeASAS( ) const { return 0.5 * As_sq * (
			
			D1() * ( (1.-cosphis())*IntExpGammal + (1.+cosphis())*IntExpGammah ) 
			+ D2() * ( -2.*intexpsin*sinphis() + 2.*intexpcos*CC() )
			
			); }
	// K9
	inline double Int_timeASASS( ) const { return Csp * sqrt(As_sq) * sqrt(Ass_sq) * (
			
			D2() * ( cos(deltaSS-deltaS)*intexpcos - sin(deltaSS-deltaS)*intexpsin*cosphis() ) 
			+ D1() * ( -0.5*(IntExpGammah-IntExpGammal)*sin(deltaSS-deltaS)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*cos(deltaSS-deltaS)*CC() )

			); }
	// K10
	inline double Int_timeA0ASS( ) const { return Csp * Csp * 0.5 * sqrt(Azero_sq) * sqrt(Ass_sq) * cos(deltaSS) * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

			); }
	// K11
	inline double Int_timeAPASS( ) const { return Csp * Csp * 0.5 * sqrt(Apara_sq) * sqrt(Ass_sq) * cos(delta_2-delta_1-deltaSS) * (
			
			D1() * ( (1.+cosphis())*IntExpGammal + (1.-cosphis())*IntExpGammah ) 
			+ D2() * ( 2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

			); }
	// K12
	inline double Int_timeATASS( ) const { return Csp * Csp * sqrt(Aperp_sq) * sqrt(Ass_sq) * (
			
			D2() * ( sin(delta_2-deltaSS)*intexpcos - cos(delta_2-deltaSS)*intexpsin*cosphis() ) 
			+ D1() * ( -0.5*(IntExpGammah-IntExpGammal)*cos(delta_2-deltaSS)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*sin(delta_2-deltaSS)*CC() )
			
			); }
	// K13
	inline double Int_timeA0AS( ) const { return Csp * sqrt(Azero_sq) * sqrt(As_sq) * (
			
			D2() * ( cos(deltaS)*intexpcos + sin(deltaS)*intexpsin*cosphis() ) 
			+ D1() * ( 0.5*(IntExpGammah-IntExpGammal)*sin(deltaS)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*cos(deltaS)*CC() )

			); }
	// K14
	inline double Int_timeAPAS( ) const { return Csp * sqrt(Apara_sq) * sqrt(As_sq) * (
			
			D2() * ( cos(delta_2-delta_1-deltaS)*intexpcos - sin(delta_2-delta_1-deltaS)*intexpsin*cosphis() ) 
			+ D1() * ( -0.5*(IntExpGammah-IntExpGammal)*sin(delta_2-delta_1-deltaS)*sinphis()  + 0.5*(IntExpGammah+IntExpGammal)*cos(delta_2-delta_1-deltaS)*CC() )

			); }
	// K15
	inline double Int_timeATAS( ) const { return Csp * 0.5 * sqrt(Aperp_sq) * sqrt(As_sq) * sin(delta_2-deltaS) * (
			
			D1() * ( (1.-cosphis())*IntExpGammal + (1.+cosphis())*IntExpGammah ) 
			+ D2() * ( -2.*intexpsin*sinphis() + 2.*intexpcos*CC() )

			); }
        
	void setupTimeAmplitudeIntegrals(PhaseSpaceBoundary*);
	
	//inline double CT1A( ) const { return 1.0 +1.31925e-02*ctheta_1 -7.43694e-02*ctheta_1*ctheta_1 -2.58816e-02*ctheta_1*ctheta_1*ctheta_1 ; }
	//inline double CT2A( ) const { return 1.0 +1.31925e-02*ctheta_2 -7.43694e-02*ctheta_2*ctheta_2 -2.58816e-02*ctheta_2*ctheta_2*ctheta_2 ; }
	// 2011
	//inline double CT1A( ) const { return 1.0 -2.18131e-03*ctheta_1 -7.63397e-02*ctheta_1*ctheta_1; }
	//inline double CT2A( ) const { return 1.0 -2.18131e-03*ctheta_2 -7.63397e-02*ctheta_2*ctheta_2; }
	//inline double PHIA( ) const { return 1.0 -3.30944e-03*phi; }
	// Use 2012 as this is most of our data
	inline double CT1A( ) const { return 1.0 -1.31493e-02*ctheta_1 -1.06010e-01*ctheta_1*ctheta_1; }
	inline double CT2A( ) const { return 1.0 -1.31493e-02*ctheta_2 -1.06010e-01*ctheta_2*ctheta_2; }
	inline double PHIA( ) const { return 1.0 +1.74833e-03*phi; }

	//Time acceptance 
	SlicedAcceptance * timeAcc ;
	
	//Configurationparameters
	bool _useTimeAcceptance ;	
	bool _useAccepWeights ;	
	vector<double> weightsVect;
	inline bool useAccepInEval() const { return _useAccepInEval ; }
	bool _useAccepInEval;
	bool _usePerEvTime;

        // The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	double w1n,w2n;	
	double returnValueInit;
};

#endif
