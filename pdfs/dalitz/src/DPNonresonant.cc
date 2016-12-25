#include "DPNonresonant.hh"

#include <iostream>

DPNonresonant::DPNonresonant()
{
}

DPNonresonant::DPNonresonant(const DPNonresonant& other) : DPMassShape(other)
{
}

DPNonresonant::~DPNonresonant()
{
}

std::complex<double> DPNonresonant::massShape(double m)
{
  (void)m;
  std::complex<double> result(1,0);
  return result;
}

