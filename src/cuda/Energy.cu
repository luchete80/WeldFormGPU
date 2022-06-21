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
    int_energy_sum[i] += dint_energy_dt * deltat;
  }      
}

__global__ inline void CalcIntEnergyKernel(Domain_d *dom){
  dom->CalcIntEnergy();
}

};