#ifndef DP_WIGNER_FUNCTION_J1
#define DP_WIGNER_FUNCTION_J1

class DPWignerFunctionJ1
{
  public:
    static double function(double cosTheta, double m, double n);
  private:
    static double d00(double cosTheta);
    static double dp10(double cosTheta);
    static double dp1p1(double cosTheta);
    static double dp1m1(double cosTheta);

};

#endif
