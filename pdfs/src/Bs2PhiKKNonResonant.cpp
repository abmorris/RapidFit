#include "TMath.h"
#include "DPHelpers.hh"
#include "Bs2PhiKKNonResonant.h"
double Bs2PhiKKNonResonant::Evaluate(double mKK)
{
  const double mK   = Bs2PhiKKComponent::mK;
  const double mBs  = Bs2PhiKKComponent::mBs;
  const double mPhi = Bs2PhiKKComponent::mphi;
  double p1_st = DPHelpers::daughterMomentum(mKK, mK, mK);
  double p3    = DPHelpers::daughterMomentum(mBs,mKK,mPhi);
  double eval = p1_st*p3/1e6;
  return TMath::IsNaN(eval) ? 0 : eval;
}
