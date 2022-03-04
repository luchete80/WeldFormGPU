#include "Domain_d.cuh"
#include <iostream>
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