#include "Domain_d.cuh"

namespace SPH{

__global__ inline void CalcAccelKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors);

__global__ inline void CalcDensIncKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors);

  
__global__ inline void CalcRateTensorsKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors);


}; //SPH