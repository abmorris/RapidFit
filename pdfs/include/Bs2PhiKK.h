#ifndef Bs2PhiKK_H
#define Bs2PhiKK_H
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
		static constexpr double mBs = 5.36677;
		static constexpr double mphi = 1.019461;
		static constexpr double mK = 0.493677;
		static constexpr double mpi = 0.13957018;
		// Datapoint type
		enum dim { _mKK_ = 0, _phi_, _ctheta_1_, _ctheta_2_, _trigger_};
		typedef std::array<double,5> datapoint_t;
		typedef std::array<std::complex<double>,2> amplitude_t; // Two complex amplitudes (B and B̅)
		static std::vector<std::string> LineShapeParameterNames(std::string name, std::string shape); // Return the necessary parameter names given a lineshape name
		static datapoint_t Parity(const datapoint_t&); // Return the datapoint with a parity operation applied
		// Simplify the case where a value and a name correspond 1:1
		struct PhysPar
		{
			PhysPar() {}
			PhysPar(PDFConfigurator* config, std::string _name) : name(config->getName(_name)), value(0) {}
			PhysPar(const PhysPar& other) : value(other.value), name(other.name) {}
			void Update(const ParameterSet* pars);
			double value;
			ObservableRef name;
		};
		static void UpdateLineshape(const std::string&, DPMassShape&, const std::vector<PhysPar>&); // Update the parameters of a resonance line shape
	protected:
		void MakePrototypeDataPoint(std::vector<std::string>&); // Create a prototype datapoint to use in MakePrototypes() in the PDFs
		std::map<dim, ObservableRef> ObservableNames; // Store observable references to retrieve numbers from the datapoints
		std::vector<bool> UseObservable; // Cache whether or not we use each dimension. Saves a lot of cycles compared to repeatedly calling ObservableNames.count()
		Bs2PhiKK::datapoint_t ReadDataPoint(DataPoint*) const; // Retrieve an array of doubles from a RapidFit Datapoint object
};
#endif

