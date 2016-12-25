#include "DPBarrierL0.hh"

DPBarrierL0::DPBarrierL0(double _radius) : DPBarrierFactor(_radius)
{
 //std::cout << "DPBarrierL0 const " << std::endl;
}

double DPBarrierL0::barrier(double p0, double p)
{
  (void)p0;
  (void)p;
  return 1;
}
