#include "Domain.h"

namespace SPH{

void __device__ inline Domain_d::CalcIntEnergy(){
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if (i<particle_count){ 
    double dint_energy_dt = m[i]/rho[i] * (
          shearstress[0]*strrate[0] + 
          shearstress[1]*strrate[1] +
          shearstress[2]*strrate[2] +
          2.0*(shearstress[3]*strrate[3] +
          shearstress[4]*strrate[4] +
          shearstress[5]*strrate[5] )
          );
    int_energy[i] += dint_energy_dt * deltat;
    int_energy_sum += int_energy[i];
  }      
}

__global__ inline void CalcIntEnergyKernel(Domain_d *dom){
  dom->CalcIntEnergy();
}


__device__ inline void Domain_d::CalcKinEnergy(
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors,
  int KernelType
	)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];

	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];
		int j = neighbors[writeOffset + k];
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		
		if (!IsFree[i]) {
			di = DensitySolid(PresEq[i], Cs[j], P0[j],p[i], rho_0[j]);
			mi = FPMassC[i] * m[j];
		} else {
			di = rho[i];
			mi = m[i];
		}
		if (!IsFree[j]) {
			dj = DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
			mj = FPMassC[j] * m[i];
		} else {
			dj = rho[j];
			mj = m[j];
		}

		double3 vij	= v[i] - v[j];
		double h_ = (h[i] + h[j])/2.0;
			
		double GK	= GradKernel(3, KernelType, rij/h_, h_);
		double K	= Kernel(3, 0, rij/h_, h_);
    
  } //nb count
  
  }
}
  
  
__global__ inline void CalcKinEnergyKernel(Domain_d *dom) {
  
}

};