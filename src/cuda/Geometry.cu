#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "Domain_d.cuh"
#include "Functions.cuh"
#include "tensor.cuh"
#include "tensor3.cu" //INLINE

__global__ void CalculateSurface(	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
																	Domain_d *dom_d, const int &id){
	dom_d->CalculateSurface(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
	id);

}
                                  
// Calculate Free Surface (for contact and heat convection)
void Domain_d::CalculateSurface(const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
                                  const int &id){
	id_free_surf = id;

	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
	int Dimension = 3; //TODO, put in another 
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];
  
  double totmass=0.;
  
    for (int k=0;k < neibcount; k++) { //Or size
      //if fixed size i = part * NB + k
      //int j = neib[i][k];
      int j = neighbors[writeOffset + k];
      //double h	= partdata->h[i]+P2->h)/2;
      double3 xij = x[i] - x[j];
      double rij = length(xij);
      double di=0.0,dj=0.0,mi=0.0,mj=0.0;

    }//
  }
  

	// Particle *P1,*P2;
	// Vec3_t xij;
	// //TODO: SAVE THIS AT THE BEGINING AND PARALLELIZE
	// double totmass=0.;
	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		// totmass += Particles[i]->Mass;
		
	// totmass /= Particles.Size();;
	// //cout << "Totmass" <<	totmass <<endl;
	
	// int maxid;
	// if (contact)
		// maxid = first_fem_particle_idx;
	// else 
		// first_fem_particle_idx = Particles.Size();
	

	// for (size_t i=0; i < maxid; i++)	{//Like in Domain::Move
		// Particles[i] -> normal = 0.;
		// Particles[i] -> ID = Particles [i] -> ID_orig;
	// }
	
	// #pragma omp parallel for schedule (static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t k=0; k<Nproc;k++) 
	// #else
	// for (int k=0; k<Nproc;k++) 
	// #endif	
	// {
		// for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			// P1	= Particles[SMPairs[k][a].first];
			// P2	= Particles[SMPairs[k][a].second];
			// xij	= P1->x - P2->x;
						
			// mi = P1->Mass;
			// mj = P2->Mass;	
			// //Eqn 3-112 Fraser Thesis
			// P1->normal += mj * xij; 
			// P2->normal -= mi * xij;					
		// } //Nproc //Pairs

	// }//Nproc
	
	// //Calculate Particle Neighbours
	
	// //TODO: Parallelize with lock

	// int surf_part =0;
	// for (size_t i=0; i < maxid; i++)	{//Like in Domain::Move
	
		// Particles[i]->normal *= 1./totmass;
		
		// if ( norm(Particles[i]->normal) >= 0.25 * Particles[i]->h && Particles[i]->Nb <= 46) {//3-114 Fraser {
			// if (!Particles[i]->not_write_surf_ID)
      // Particles[i]->ID = id;
			// surf_part++;
		// }
	// }
	// //cout << "Surface particles: " << surf_part<<endl;
  
  } //i < particle count
}

#endif