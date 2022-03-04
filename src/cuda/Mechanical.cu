#include "Domain_d.cuh"
namespace SPH{

//Called by __global__
//void __device__ Domain_d::WholeVelocity() {
    // //Apply a constant velocity to all particles in the initial time step
    // if (norm(BC.allv)>0.0 || BC.allDensity>0.0) {
    	// Vec3_t vel = 0.0;
    	// double den = 0.0;

    	// for (int i=0 ; i<Particles.Size() ; i++) {
		// AllCon(Particles[i]->x,vel,den,BC);
    		// if (Particles[i]->IsFree && norm(BC.allv)>0.0) {
			// Particles[i]->v		= vel;
 		// }
    		// if (IsFree[i] && BC.allDensity>0.0) {
			// rho[i]	= den;
			// p[i]		= EOS(Particles[i]->PresEq[i], Particles[i]->Cs[i], Particles[i]->P0,Particles[i]->rho[i], RefDensity);
    		// }
    	// }
    // }
//}


// __device__ PartData_d::PrimaryComputeAcceleration(){
	// // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
	
	
	// // Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
// }

// void __global__ MechSolveKernel (double dt, PartData_d *partdata) {
	// int i = threadIdx.x+blockDim.x*blockIdx.x;
	// dTdt[i] = 0.;
	
	// int neibcount;
	// #ifdef FIXED_NBSIZE
	// neibcount = neib_offs[i];
	// #else
	// neibcount =	neib_offs[i+1] - neib_offs[i];
	// #endif
	// printf("Solving\n");
	// for (int k=0;k < neibcount;k++) { //Or size
		// //if fixed size i = part * NB + k
		// //int j = neib[i][k];
		// int j = NEIB(i,k);
		// printf("i,j\n",i,j);
		// double3 xij; 
		// xij = x[i] - x[j];
		// double h_ = (h[i] + h[j])/2.0;
		// double nxij = length(xij);
		
		// double GK	= GradKernel(3, 0, nxij/h_, h_);
		// //		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
		// //   mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
		// dTdt[i] += m[j]/rho[j]*( 4.0*k_T[i]*k_T[j]/(k_T[i]+k_T[j]) * (T[i] - T[j])) * dot( xij , GK*xij )/(nxij*nxij);
	// }
	// dTdt[i] *=1/(rho[i]*cp[i]);
// }

__device__ __forceinline__ void Domain_d::LastComputeAcceleration(){
	//partdata->();
}

void Domain_d::MechSolve(const double &tf){

	int N = particle_count;
	int threadsPerBlock = 256; //Or BlockSize
	int blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
  Time =0.;
	
	isfirst_step =true;
	
	//This was in Original LastCompAcceleration
	CalcForcesKernel	<<<blocksPerGrid,threadsPerBlock >>>(this);
	cudaDeviceSynchronize(); //REQUIRED!!!!
	
	//TODO: Pass toPartData
	//CalcForcesMember	<<<blocksPerGrid,threadsPerBlock >>>(partdata);
	//MechSolveKernel<<< >>>();
	
	
}

};//SPH