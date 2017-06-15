#include "DPWignerFunctionJ1.hh"
#include <cmath>

double DPWignerFunctionJ1::function(double cosTheta, double mm, double nn)
{

  int m=(int)mm;
  int n=(int)nn;
  if ( m == 0 )
  {
    switch (n)
    {
      case 0: return d00(cosTheta);
              break;
      case 1: return -dp10(cosTheta);
              break;
      case -1: return dp10(cosTheta);
              break;
    }
  }
  else if ( m == 1 )
  {
    switch (n)
    {
      case 0: return dp10(cosTheta);
              break;
      case 1: return dp1p1(cosTheta);
              break;
      case -1: return dp1m1(cosTheta);
              break;
    }
  }
  else if ( m == -1 )
  {
    switch (n)
    {
      case 0: return -dp10(cosTheta);
              break;
      case 1: return dp1m1(cosTheta);
              break;
      case -1: return dp1p1(cosTheta);
              break;
    }
  }

  return -1000; // Give crazy number, alternatively we can exit or throw exception
}

double DPWignerFunctionJ1::d00(double cosTheta)
{
  return cosTheta;
}

double DPWignerFunctionJ1::dp10(double cosTheta)
{
  double sinTheta=std::sqrt(1-cosTheta*cosTheta);
  return -sinTheta/std::sqrt(2.);
}

double DPWignerFunctionJ1::dp1p1(double cosTheta)
{
  return 0.5*(1+cosTheta);
}

double DPWignerFunctionJ1::dp1m1(double cosTheta)
{
  return 0.5*(1-cosTheta);
}

