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

  
  
	// std::vector < Mat3_t> temp(Particles.Size());
	// Mat3_t m,mt[2];
	
	// //cout << "Applying grad corr"<<endl;
	// //#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	// for ( size_t k = 0; k < Nproc ; k++) {
		// Particle *P1,*P2;
		// Vec3_t xij;
		// double h,GK;
		// //cout << "SMPairs[k].Size()"<<SMPairs[k].Size()<<endl;
		// //TODO: DO THE LOCK PARALLEL THING
		// for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			// P1	= Particles[SMPairs[k][a].first];
			// P2	= Particles[SMPairs[k][a].second];
			// xij	= P1->x - P2->x;
			// h	= (P1->h+P2->h)/2.0;
			// GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			// di = P1->Density; mi = P1->Mass;
			// dj = P2->Density; mj = P2->Mass;
		
			// Dyad (Vec3_t(GK*xij),xij,m);
			// mt[0] = mj/dj * m;
			// mt[1] = mi/di * m;
			// //cout << "mt"<<mt[0]<<endl;
			// //omp_set_lock(&P1->my_lock);
			// //SIGN IS NEGATIVE (IF POSITIVE, GRADIENT SIGN IS OPPOSITE)
			
			// temp[SMPairs[k][a].first]  = temp[SMPairs[k][a].first]  - mt[0];  
			// temp[SMPairs[k][a].second] = temp[SMPairs[k][a].second] - mt[1];
		// }
	// }//Nproc


	// //cout << "Inverting"<<endl;
	// //#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	// //cout << "Inverting"<<endl;
	// int max_id = Particles.Size();
	// if (contact)
		// max_id = first_fem_particle_idx[0];
		
	// for (int i=0; i<max_id; i++){
		// // cout << "part "<<i<<endl;
		// //cout << "x: "<<Particles[i]->x<<endl;
		// // cout << "nb: "<<Particles[i]->Nb<<endl;
		// // if (!Particles[i]->IsFree) cout << "Fixed"<<endl;
		// //cout << "temp "<<temp[i]<<endl;
		// if (Dimension == 2)
			// temp[i](2,2) = 1;
		// /** Inverse.*/
		// //inline void Inv (Mat3_t const & M, Mat3_t & Mi, double Tol=1.0e-10)}	
		// if (Particles[i]->IsFree){
			// if (!InvLog(temp[i],m)) {//TODO: CHANGE THIS
        // m = I;
        // nulldetcount++;
      // }

			// Particles[i] ->gradCorrM = m;
			// //cout << "Corr Matrix: " << m <<endl;
		// } else {
			// Particles[i] ->gradCorrM = I;
		// }
	// }
  // if (nulldetcount>0) cout << nulldetcount << " particles with correction matrix determinant."<<endl;	
  } // NEIBCOUNT
  
  ToFlatSymPtr(Inv(temp),gradCorrM,9*i);
  
  }//i < particle_count 
}

}; //SPH