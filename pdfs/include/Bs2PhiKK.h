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
		// Fixed mass values
		static constexpr double mBs  = 5.36689;
		static constexpr double mphi = 1.019460;
		static constexpr double mK   = 0.493677;
		static constexpr double mpi  = 0.13957061;
		typedef std::array<double,5> datapoint_t; // Datapoint type
		typedef std::array<std::complex<double>,2> amplitude_t; // Two complex amplitudes (B and B̅)
		static std::vector<std::string> LineShapeParameterNames(std::string name, std::string shape); // Return the necessary parameter names given a lineshape name
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
		void MakePrototypeDataPoint(std::vector<std::string>&); // Create a prototype datapoint to use in MakePrototypes() in the PDFs
		std::vector<ObservableRef> ObservableNames;// Datapoint stuff: K+K− mass and helicity angles
		Bs2PhiKK::datapoint_t ReadDataPoint(DataPoint*) const; // Retrieve an array of doubles from a RapidFit Datapoint object
};

#endif

