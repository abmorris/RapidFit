#include "DPWignerFunctionJ2.hh"
#include <cmath>

double DPWignerFunctionJ2::function(double cosTheta, int m)
{
  switch (m)
  {
    case 0: return d00(cosTheta);
            break;
    case 1: return dp10(cosTheta);
            break;
    case -1: return -dp10(cosTheta);
            break;
    case 2: return dp20(cosTheta);
            break;
    case -2: return -dp20(cosTheta);
            break;
    default: return -1000; // Give crazy number, alternatively we can exit or throw exception
  }
}

double DPWignerFunctionJ2::d00(double cosTheta)
{
  return 1.5 * cosTheta * cosTheta - 0.5;
}

double DPWignerFunctionJ2::dp10(double cosTheta)
{
  double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
  return -std::sqrt(1.5) * sinTheta * cosTheta;
}
double DPWignerFunctionJ2::dp20(double cosTheta)
{
  double sinTheta_sq = 1 - cosTheta * cosTheta;

  return (std::sqrt(6.)/4.)*sinTheta_sq;
}
