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
	
	if ( i < first_fem_particle_idx ) { //In Contact Surface Particles, normal are updated different way
    normal[i] = make_double3(0.,0.,0.);

    int neibcount = particlenbcount[i];
    int nbcount_corr = 0; //WITHOUT CONTACT SURFACE!!
    const uint writeOffset = neighborWriteOffsets[i];
    
    
      for (int k=0;k < neibcount; k++) { //Or size
        //if fixed size i = part * NB + k
        //int j = neib[i][k];
        int j = neighbors[writeOffset + k];
        //double h	= partdata->h[i]+P2->h)/2;
        double3 xij = x[i] - x[j];
        if (ID[j]!=contact_surf_id){  //EXCLUDE RIGID PAIRS!
          normal[i] += m[j] * xij; 
          //printf("particle %d Nb %d xij: %f %f %f \n", i, j, xij.x, xij.y, xij.z);
          nbcount_corr++;
        }

      }//
      //printf("particle %d normal : %f %f %f , nb %d\n", i, normal[i].x, normal[i].y, normal[i].z, nbcount_corr);
      //normal[i]*= ((double)particle_count/(totmass *(double)neibcount)); //Attention parenthesis, if not it crashes
      normal[i]*= ((double)particle_count/(totmass *(double)nbcount_corr)); //Attention parenthesis, if not it crashes
      if ( length(normal[i]) >= 0.25 * h[i] && nbcount_corr <= 46) {//3-114 Fraser {
        // //if (!Particles[i]->not_write_surf_ID)
        // printf("I: %d\n",i);
        ID[i] = id_free_surf; //THIS CRASH IS ASSIGNED BY PARAMETER
        // //surf_part++;
      }
  
  }//i < particle_count
  
}

}; //SPH

#endif