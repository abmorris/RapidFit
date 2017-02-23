#ifndef Bs2PhiKKHelpers_H
#define Bs2PhiKKHelpers_H
#include <array>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include "PDFConfigurator.h"
#include "ParameterSet.h"
#include "DPMassShape.hh"

class Bs2PhiKK
{
	public:
		Bs2PhiKK(PDFConfigurator*);
		Bs2PhiKK(const Bs2PhiKK&);
		~Bs2PhiKK() {}
		static double mBs;
		static double mK;
		static double mpi;
		typedef std::array<double,4> datapoint_t;
		typedef std::array<std::complex<double>,2> amplitude_t;
		static std::vector<std::string> LineShapeParameterNames(std::string name, std::string shape);
		static bool IsPhysicalDataPoint(const Bs2PhiKK::datapoint_t&); // Return whether or not this datapoint makes sense
		// Simplify the case where a value and a name correspond 1:1
		struct PhysPar
		{
			// Construct this however you want
			PhysPar() {}
			PhysPar(ObservableRef _name) : name(_name), value(0) {}
			PhysPar(ObservableRef _name, double _value) : name(_name), value(_value) {}
			PhysPar(PDFConfigurator* config, std::string _name) : name(config->getName(_name)), value(0) {}
			PhysPar(PDFConfigurator* config, std::string _name, double _value) : name(config->getName(_name)), value(_value) {}
			PhysPar(const PhysPar& other) : value(other.value), name(other.name) {}
			void Update(const ParameterSet* pars);
			double value;
			ObservableRef name;
		};
		static void UpdateLineshape(const std::string&, DPMassShape&, const std::vector<PhysPar>&); // Update the parameters of a resonance line shape
	protected:
		ObservableRef mKKName, ctheta_1Name, ctheta_2Name, phiName; // Datapoint stuff: K+K− mass and helicity angles
		Bs2PhiKK::datapoint_t ReadDataPoint(DataPoint*) const; // Retrieve an array of doubles from a RapidFit Datapoint object
};

#endif

