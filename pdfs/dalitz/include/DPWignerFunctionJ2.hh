#ifndef DP_WIGNER_FUNCTION_J2
#define DP_WIGNER_FUNCTION_J2

class DPWignerFunctionJ2
{
  public:
    static double function(double cosTheta, int m);
  private:
    static double d00(double cosTheta);
    static double dp10(double cosTheta);
    static double dp20(double cosTheta);
};

#endif
