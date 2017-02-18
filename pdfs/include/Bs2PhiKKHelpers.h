#ifndef Bs2PhiKKHelpers_H
#define Bs2PhiKKHelpers_H
#include <array>
#include <complex>
#include <iostream>
#include <string>
#include "PDFConfigurator.h"
#include "ParameterSet.h"

class Bs2PhiKK
{
	public:
		static double mBs;
		static double mK;
		static double mpi;
		typedef std::array<double,4> datapoint_t;
		typedef std::array<std::complex<double>,2> amplitude_t;
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
};

#endif

