#include "DPWignerFunctionJ1.hh"
#include <cmath>

double DPWignerFunctionJ1::function(double cosTheta, int m)
{

  switch (m)
  {
    case 0: return d00(cosTheta);
            break;
    case 1: return dp10(cosTheta);
            break;
    case -1: return -dp10(cosTheta);
            break;
    default: return -1000; // Give crazy number, alternatively we can exit or throw exception
  }
}

double DPWignerFunctionJ1::d00(double cosTheta)
{
  return cosTheta;
}

double DPWignerFunctionJ1::dp10(double cosTheta)
{
  double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
  return -sinTheta / std::sqrt(2.);
}

