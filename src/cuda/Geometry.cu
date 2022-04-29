#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "Domain_d.cuh"
#include "Functions.cuh"
#include "tensor.cuh"
#include "tensor3.cu" //INLINE

namespace SPH {
  
__global__ void CalculateSurfaceKernel(Domain_d *dom_d,	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
																	const int &id, const double &totmass){
	dom_d->CalculateSurface(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
	id, dom_d->totmass);

}
                                  
// Calculate Free Surface (for contact and heat convection)
void __device__ inline Domain_d::CalculateSurface(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors,
                                                  const int &id, const double &totmass){
	//id_free_surf = id;

	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
    int Dimension = 3; //TODO, put in another 
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
        normal[i] += m[j] * xij; 

      }//

      normal[i]*= 1./totmass;
      
      if ( length(normal[i]) >= 0.25 * h[i] && neibcount <= 46) {//3-114 Fraser {
        //if (!Particles[i]->not_write_surf_ID)
        ID[i] = id;
        //surf_part++;
      }
    
  
  }//i < particle_count
  
}

}; //SPH

#endif