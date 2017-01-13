#ifndef DP_BARRIER_FUNCTION
#define DP_BARRIER_FUNCTION
#include <memory>
class DPBarrierFactor
{
	public:
		DPBarrierFactor();
		DPBarrierFactor(unsigned spin, double radius);
		DPBarrierFactor(const DPBarrierFactor& other) : function(other.function), radius(other.radius) {}
		DPBarrierFactor(DPBarrierFactor&& other) : function(std::move(other.function)), radius(std::move(other.radius)) {}
		DPBarrierFactor& operator=(const DPBarrierFactor& other) {radius = other.radius; function = other.function; return *this;}
		DPBarrierFactor& operator=(DPBarrierFactor&& other) {radius = std::move(other.radius); function = std::move(other.function); return *this;}
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
