#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "Domain_d.cuh"
#include "Functions.cuh"
#include "tensor.cuh"
#include "tensor3.cu" //INLINE

namespace SPH {
  
__global__ inline void CalculateSurfaceKernel(Domain_d *dom_d,	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
																	/*const int &id, */const double &totmass){
	dom_d->CalculateSurface(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
	/*id, */dom_d->totmass);

}
                                  
// Calculate Free Surface (for contact and heat convection)
void __device__ inline Domain_d::CalculateSurface(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors,
                                                  /*const int &id, */const double &totmass){
	//id_free_surf = id;

	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
    normal[i] = make_double3(0.,0.,0.);

    int neibcount = particlenbcount[i];
    const uint writeOffset = neighborWriteOffsets[i];
    
    
      for (int k=0;k < neibcount; k++) { //Or size
        //if fixed size i = part * NB + k
        //int j = neib[i][k];
        int j = neighbors[writeOffset + k];
        //double h	= partdata->h[i]+P2->h)/2;
        double3 xij = x[i] - x[j];
        normal[i] += m[j] * xij; 

      }//
      normal[i]*= ((double)particle_count/(totmass *(double)neibcount)); //Attention parenthesis, if not it crashes

      if ( length(normal[i]) >= 0.25 * h[i] && neibcount <= 46) {//3-114 Fraser {
        // //if (!Particles[i]->not_write_surf_ID)
        // printf("I: %d\n",i);
        ID[i] = id_free_surf; //THIS CRASH IS ASSIGNED BY PARAMETER
        // //surf_part++;
      }
  
  }//i < particle_count
  
}

}; //SPH

#endif