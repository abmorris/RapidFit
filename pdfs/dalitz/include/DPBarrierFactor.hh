#ifndef DP_BARRIER_FUNCTION
#define DP_BARRIER_FUNCTION

class DPBarrierFactor
{
	public:
		DPBarrierFactor();
		DPBarrierFactor(unsigned spin, double radius);
		~DPBarrierFactor() {}
		double barrier(const double p0, const double p) const;
		static double FunctionL0(const double z);
		static double FunctionL1(const double z);
		static double FunctionL2(const double z);
		static double FunctionL3(const double z);
	private:
		double (*function)(const double z);
		double radius;  // Blatt-Weiskopf radius
};

#endif
