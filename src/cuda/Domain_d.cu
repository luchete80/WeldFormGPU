#include "Domain_d.cuh"
#include "Functions.cuh"
#include "Domain.h"

#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

//Else (offset)
//Allocating from host
namespace SPH {
void Domain_d::SetDimension(const int &particle_count){
	this->particle_count = particle_count;
	//Allocae arrays (as Structure of arryays, SOA)

	cudaMalloc((void **)&x, particle_count * sizeof (double3));
	cudaMalloc((void **)&v, particle_count * sizeof (Vector));
	cudaMalloc((void **)&a, particle_count * sizeof (Vector));

	cudaMalloc((void **)&h, 	particle_count * sizeof (double));
	cudaMalloc((void **)&m, 	particle_count * sizeof (double));
	cudaMalloc((void **)&rho, particle_count * sizeof (double));
	
	//THERMAL
	cudaMalloc((void **)&k_T, particle_count * sizeof (Vector));
	cudaMalloc((void **)&cp_T, particle_count * sizeof (Vector));
		
	cudaMalloc((void **)&T		, particle_count * sizeof (double));
	cudaMalloc((void **)&Ta		, particle_count * sizeof (double));
	cudaMalloc((void **)&Tb		, particle_count * sizeof (double));
	
	cudaMalloc((void **)&BC_T, particle_count * sizeof (int));
	
	//Host things
	T_h =  new double [particle_count];
	
	cudaMalloc((void **)&dTdt	, particle_count * sizeof (double));
	//printf("Size of dTdt: %d, particle count %d\n",sizeof(dTdt)/sizeof (double),particle_count);

	//Nb data
	cudaMalloc((void **)&neib_offs	, (particle_count + 1) * sizeof (int));
	cudaMalloc((void **)&neib_part	, (particle_count * 400) * sizeof (int));
	
	//cudaMalloc((void **)&partdata, sizeof(PartData_d));
	
	//cudaMalloc((void **)&partdata->dTdt,particle_count * sizeof (double)); //TODO, pass to PartData
	
	// cudaMalloc((void**)&ppArray_a, 10 * sizeof(int*));
	// for(int i=0; i<10; i++) {
		// cudaMalloc(&someHostArray[i], 100*sizeof(int)); /* Replace 100 with the dimension that u want */
	// }

	
	//To allocate Neighbours, it is best to use a equal sized double array in order to be allocated once
}

void Domain_d::CheckData(){
	printf("dTdt partdta: %d",sizeof(this->dTdt)/sizeof(double));
	printf("dTdt[200] %f",dTdt[200]);
	printf("neibpart %f",neib_part[300000]);
	//dom->CheckData();
}

__global__ void CheckData(Domain_d *dom){
	//printf("dTdt partdta: %d",sizeof(dom->partdata->dTdt)/sizeof(double));
	dom->CheckData();
}

void Domain_d::Set_h(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->h, k_, size, cudaMemcpyHostToDevice);
	delete k_;
}

void Domain_d::SetConductivity(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->k_T, k_, size, cudaMemcpyHostToDevice);
	delete k_;
}

void Domain_d::SetDensity(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->rho, k_, size, cudaMemcpyHostToDevice);
	delete k_;
}

void Domain_d::SetHeatCap(const double &cp){
	double *cp_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		cp_[i] = cp;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->cp_T, cp_, size, cudaMemcpyHostToDevice);
	delete cp_;
}

// // Templatize data type, and host and device vars (of this type)
// template <typename T> copydata (const Domain &d, T *var_h, T *var_d){
	// T *var_h =  (Vector *)malloc(dom.Particles.size());
	// for (int i=0;i<dom.Particles.size();i++){
		// var_h[i] = dom.Particles[i]->T;
	// }
	// int size = dom.Particles.size() * sizeof(Vector);
	// cudaMemcpy(this->T, T, size, cudaMemcpyHostToDevice);
// }

//TEMPORARY, UNTIL EVERYTHING WILL BE CREATED ON DEVICE
void __host__ Domain_d::CopyData(const Domain& dom){
	
	//TODO TEMPLATIZE THIS!!
	double *T =  (double *)malloc(dom.Particles.size());
	for (int i=0;i<dom.Particles.size();i++){
		T[i] = dom.Particles[i]->T;
	}
	int size = dom.Particles.size() * sizeof(double);
	cudaMemcpy(this->T, T, size, cudaMemcpyHostToDevice);

	// for (int i=0;i<dom.Particles.size();i++){
		// T[i] = dom.Particles[i]->cp_T;
	// }
	// int size = dom.Particles.size() * sizeof(double);
	// cudaMemcpy(this->cp_T, T, size, cudaMemcpyHostToDevice);
	
}

void __device__ Domain_d::CalcThermalTimeStep(){
	deltat = 0.3*h[0]*h[0]*rho[0]*cp_T[0]/k_T[0];
	printf("Time Step: %f\n",deltat);
}

void __global__ CalcConvHeatKernel (double *dTdt,
																		double *m, double *rho, double *cp_T,
																		double *T, double &T_inf,
																		int *BC_T,
																		double &h_conv, int count) {
	
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	if ( i < count ) {	
		if ( BC_T[i] == 1) {
			double dS2 = pow(m[i]/rho[i],0.666666666);
			//printf("dS2 %f\n",dS2);
			//cout << Particles[i]->Density<<endl;
			//Fraser Eq 3.121
			//double q_conv = rho[i] * h_conv * dS2 * (T_inf[i] - T[i])/m[i]; //Fraser Eq - 3.121, J/(m3s)

			//dTdt[i] += 1./(rho[i] * cp_T[i] ) * q_conv; //J /*+ Particles[i]->q_source + Particles[i]->q_plheat);	*/
			
			
			double q_conv = h_conv * dS2 * (T_inf - T[i])/cp_T[i];
			//dTdt[i] = dTdt[i] + q_conv;
		
			// if (Particles[i]->q_conv>max){
				// max= Particles[i]->q_conv;
				// imax=i;
			// }
			//cout << "Particle  "<<Particles[i]->Mass<<endl;
		}
	}
}
//Thread per particle
//dTdt+=1/cp* (mass/dens^2)*4(k)
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

//This is the future solver, ir order to not pass so much data
void __global__ ThermalSolveKernel (double dt, PartData_d *partdata){
	
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
void Domain_d::ThermalSolve(const double &tf){
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
	double T_inf = 500.;
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
												T, T_inf,
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


//NEXT SOLVER
// void Domain_d::ThermalSolve(const double &tf){
	
	
// }

Domain_d::~Domain_d(){
	
		cudaFree(a);		
		cudaFree(v);

		cudaFree(h);		
		cudaFree(m);
		cudaFree(rho);

		cudaFree(neib_offs);
		cudaFree(neib_part);		
}

};//SPH
    // // Create host pointer to array-like storage of device pointers
    // Obj** h_d_obj = (Obj**)malloc(sizeof(Obj*) * 3); //    <--------- SEE QUESTION 1
    // for (int i = 0; i < 3; i++) {
        // // Allocate space for an Obj and assign
        // cudaMalloc((void**)&h_d_obj[i], sizeof(Obj));
        // // Copy the object to the device (only has single scalar field to keep it simple)
        // cudaMemcpy(h_d_obj[i], &(h_obj[i]), sizeof(Obj), cudaMemcpyHostToDevice);
    // }

    // /**************************************************/
    // /* CREATE DEVICE ARRAY TO PASS POINTERS TO KERNEL */
    // /**************************************************/

    // // Create a pointer which will point to device memory
    // Obj** d_d_obj = NULL;
    // // Allocate space for 3 pointers on device at above location
    // cudaMalloc((void**)&d_d_obj, sizeof(Obj*) * 3);
    // // Copy the pointers from the host memory to the device array
    // cudaMemcpy(d_d_obj, h_d_obj, sizeof(Obj*) * 3, cudaMemcpyHostToDevice);
