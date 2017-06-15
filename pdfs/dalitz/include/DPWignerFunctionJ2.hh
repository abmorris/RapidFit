#ifndef DP_WIGNER_FUNCTION_J2
#define DP_WIGNER_FUNCTION_J2

class DPWignerFunctionJ2
{
  public:
    static double function(double cosTheta, double m, double n);
  private:
    static double d00(double cosTheta);
    static double dp10(double cosTheta);
    static double dp1p1(double cosTheta);
    static double dp1m1(double cosTheta);
    static double dp2p2(double cosTheta);
    static double dp2p1(double cosTheta);
    static double dp20(double cosTheta);
    static double dp2m1(double cosTheta);
    static double dp2m2(double cosTheta);
};

#endif
