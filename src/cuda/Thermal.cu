#include "Domain_d.cuh"
#include <iostream>
#include "Functions.cuh"
#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

using namespace std;

namespace SPH{
//This is the future solver, ir order to not pass so much data
void __global__ ThermalSolveKernel (double dt, PartData_d *partdata){
	partdata->ThermalSolveKernel(dt);
}

__global__ void TempCalcLeapfrogFirst(double *T, double *Ta, double *Tb, //output
																			double *dTdt, double dt, int count){//input
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	if ( i < count ) {
	Ta[i] = T[i] - dt/2.0*dTdt[i];
	}
}
__global__ void TempCalcLeapfrog     (double *T, double *Ta, double *Tb,
																		double *dTdt, double dt, int count){
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	if ( i < count ) {
	Tb[i]  = Ta[i];
	Ta[i] += dTdt[i] * dt;
	T [i] = ( Ta[i] + Tb[i] ) / 2.;
	}
}

////////////////////////////////////////////////////
//// THIS IS NOT NECESSARY TO PASS THIS IN DOMAIN //
////////////////////////////////////////////////////
void __global__ CalcConvHeatKernel (double *dTdt,
																		double *m, double *rho, double *cp_T,
																		double *T, double T_inf,
																		int *BC_T,
																		double h_conv, int count) {
	
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	if ( i < count ) {	
		if ( BC_T[i] == 1) {
			double dS2 = pow(m[i]/rho[i],0.666666666);
			//printf("dS2 %f\n",dS2);
			//cout << Particles[i]->Density<<endl;
			//Fraser Eq 3.121
			//double q_conv = rho[i] * h_conv * dS2 * (T_inf[i] - T[i])/m[i]; //Fraser Eq - 3.121, J/(m3s)

			//dTdt[i] += 1./(rho[i] * cp_T[i] ) * q_conv; //J /*+ Particles[i]->q_source + Particles[i]->q_plheat);	*/
			
			
			double q_conv =  /*rho[i]  * */h_conv * dS2 * (T_inf - T[i])/m[i];
			dTdt[i] += 1./(/*rho[i] * */cp_T[i] ) * q_conv;
		
			// if (Particles[i]->q_conv>max){
				// max= Particles[i]->q_conv;
				// imax=i;
			// }
			//cout << "Particle  "<<Particles[i]->Mass<<endl;
		}
	}
}

void __global__ ThermalSolveKernel (double *dTdt,
																		double3 *x, double *h,
																		double *m, double *rho,
																		double *T, double *k_T, double *cp, 
																		int *neib_part, int *neib_offs,/*orcount*/
																		int count) {

//	printf("searching nb..\n");	
	
	int i = threadIdx.x+blockDim.x*blockIdx.x;

	if ( i < count ) {
		dTdt[i] = 0.;

		int neibcount;
		#ifdef FIXED_NBSIZE
		neibcount = neib_offs[i];
		#else
		neibcount =	neib_offs[i+1] - neib_offs[i];
		#endif
		// printf("neibcount %d\n",neibcount);
		// printf("Nb indexed,i:%d\n",i);
		for (int k=0;k < neibcount;k++) { //Or size
			// //if fixed size i = part * NB + k
			// //int j = neib[i][k];
			int j = NEIB(i,k);
			//printf("i,j: %d,%d\n",i,j);
			double3 xij; 
			xij = x[i] - x[j];
			//printf("xij: %f,%f,%f:\n",x[i].x,x[i].y,x[i].z);
			double h_ = (h[i] + h[j])/2.0;
			double nxij = length(xij);
			
			double GK	= GradKernel(3, 0, nxij/h_, h_);
			//printf("i, j, rho, GK, nxij,h: %d, %d, %f, %f, %f, %f\n",i, j, rho[j], GK,nxij,h_);
			//		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
			//   mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
			dTdt[i] += m[j]/rho[j]*( 4.0*k_T[i]*k_T[j]/(k_T[i]+k_T[j]) * (T[i] - T[j])) * dot( xij , GK*xij )/(nxij*nxij); //Fraser, Eqn 3.99
		}
		dTdt[i] *= 1./(rho[i]*cp[i]);
	}
}

////////////////////////////////////////////////////
//// ALTERNATIVE KERNEL TO "ONLY" PASS PARTICLE DATA
void __device__ PartData_d::ThermalSolveKernel (double dt) {

//	printf("searching nb..\n");	
	
	int i = threadIdx.x+blockDim.x*blockIdx.x;

	if ( i < particle_count ) {
		dTdt[i] = 0.;

		int neibcount;
		#ifdef FIXED_NBSIZE
		neibcount = neib_offs[i];
		#else
		neibcount =	neib_offs[i+1] - neib_offs[i];
		#endif
		// printf("neibcount %d\n",neibcount);
		// printf("Nb indexed,i:%d\n",i);
		for (int k=0;k < neibcount;k++) { //Or size
			// //if fixed size i = part * NB + k
			// //int j = neib[i][k];
			int j = NEIB(i,k);
			//printf("i,j: %d,%d\n",i,j);
			double3 xij; 
			xij = x[i] - x[j];
			//printf("xij: %f,%f,%f:\n",x[i].x,x[i].y,x[i].z);
			double h_ = (h[i] + h[j])/2.0;
			double nxij = length(xij);
			
			double GK	= GradKernel(3, 0, nxij/h_, h_);
			//printf("i, j, rho, GK, nxij,h: %d, %d, %f, %f, %f, %f\n",i, j, rho[j], GK,nxij,h_);
			//		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
			//   mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
			dTdt[i] += m[j]/rho[j]*( 4.0*k_T[i]*k_T[j]/(k_T[i]+k_T[j]) * (T[i] - T[j])) * dot( xij , GK*xij )/(nxij*nxij); //Fraser, Eqn 3.99
		}
		dTdt[i] *= 1./(rho[i]*cp_T[i]);
	}
}



//Originally in Particle::TempCalcLeapfrog
//Host function
__host__ void Domain_d::ThermalSolve(const double &tf){
	int N = particle_count;
	int threadsPerBlock = 256; //Or BlockSize
	int blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
  Time =0.;
	
	isfirst_step =true;
	
	int step = 0;
	
	clock_t clock_beg;
	double time_spent;
	
	clock_beg = clock();
	double t_inf = 500.;
	double h_conv = 100.;
	
	double t_out,dt_out;
	t_out = dt_out = 0.1;
	while (Time<tf) {
	// cout << "Callign Kernel"<<endl;
	// cout << "blocksPerGrid (Blocksize)"<<blocksPerGrid<<endl;
	// cout << "threads per block (grid size)"<<threadsPerBlock<<endl;
	
		ThermalSolveKernel<<<blocksPerGrid,threadsPerBlock>>>(dTdt,	
																		x, h, //Vector has some problems
																		m, rho, 
																		T, k_T, cp_T,
																		neib_part, neib_offs,
																		particle_count);
		cudaDeviceSynchronize(); //REQUIRED!!!!
	
		CalcConvHeatKernel <<< blocksPerGrid,threadsPerBlock >>> (dTdt,
												m, rho, cp_T,
												T, t_inf,
												BC_T,
												h_conv,
												particle_count);
		cudaDeviceSynchronize();
												
		//cout << "Kernel called"<<endl;
		 if (isfirst_step) {
			TempCalcLeapfrogFirst<<< blocksPerGrid,threadsPerBlock >>>(T, Ta, Tb,
																			 dTdt, deltat, particle_count);	
			cudaDeviceSynchronize(); //After ANY COMMAND!!!!
			isfirst_step = false;
		} else {
			TempCalcLeapfrog <<< blocksPerGrid,threadsPerBlock >>>(T, Ta, Tb,
																			 dTdt, deltat, particle_count);		
			cudaDeviceSynchronize();
		}
		
		Time += deltat;
		if (Time >= t_out) {
	//		cout << "Copying to host"<<endl;
			cudaMemcpy(T_h, T, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);	
			double max=0;
			for (int i=0;i<particle_count;i++){
				if (T_h[i]>max) max = T_h[i];
			}
			t_out += dt_out;
			cout << "dTdt max\n"<<max<<endl;
		}

		step ++;
	}//main time while
	
	time_spent = (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
	
	printf("Total steps: %d, time spent %f\n",step, time_spent);
	
}//Thermal Solve

};