#include "Domain_d.cuh"
#include "Functions.cuh"
// #include <iostream>

// #include <chrono>
// //#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
// #include <ctime> //Clock
#include "tensor3.cu" //INLINE
// #include "Interaction.cu"

// #include "InteractionAlt.cuh"

// #include "Geometry.cu"
// #include "Contact.cu"

// #include "cuNSearch.h"
// #include "Material.cuh"

namespace SPH {
  
__global__ inline void CalcGradCorrMatrixKernel (Domain_d *dom,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors,
  int KernelType){
  dom->CalcGradCorrMatrix(particlenbcount,neighborWriteOffsets,neighbors,KernelType);
}
//New, for Bonet gradient correction
__device__ void Domain_d::CalcGradCorrMatrix (
  const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors,
  int KernelType) {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;

  int i = threadIdx.x + blockDim.x*blockIdx.x;
  
  tensor3 m_, mt[2], temp;

  if (i<particle_count){
    int neibcount = particlenbcount[i];
    const uint writeOffset = neighborWriteOffsets[i];
	
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];

		int j = neighbors[writeOffset + k];

		double3 xij	= x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
    
		double h_ = (h[i] + h[j])/2.0;
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, KernelType, rij/h_, h_);

    di = rho[i];    dj = rho[j];
    mi = m[i];      mj = m[j];
    m_ = Dyad (GK*xij,xij); 
    temp = temp  - mj/dj * m_;

  } // NEIBCOUNT
    // if (i==0){
    // printf ("corr matrx\n");
    // print (temp);
    // printf("inv temp\n");
    // print (Inv(temp));
    // }  
  //Is Symmetric??
  ToFlatPtr(Inv(temp),gradCorrM,9*i);
  //ToFlatSymPtr(Inv(temp),gradCorrM,9*i);
  
  }//i < particle_count 
}

}; //SPH