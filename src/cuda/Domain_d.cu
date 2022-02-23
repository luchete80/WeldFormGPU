#include "Domain_d.cuh"
#include "Functions.cuh"
#include "Domain.h"

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
	cudaMalloc((void **)&dTdt	, particle_count * sizeof (double));
	printf("Size of dTdt: %d, particle count %d\n",sizeof(dTdt)/sizeof (double),particle_count);

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

void Domain_d::SetConductivity(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->k_T, k_, size, cudaMemcpyHostToDevice);
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

//Thread per particle
//dTdt+=1/cp* (mass/dens^2)*4(k)
void __global__ ThermalSolveKernel (double *dTdt,
																		double3 *x, double *h,
																		double *m, double *rho,
																		double *T, double *k_T, double *cp, 
																		int *neib_part, int *neib_offs/*orcount*/) {

	//printf("searching nb..\n");	
	//int i = threadIdx.x+blockDim.x*blockIdx.x;
	// dTdt[i] = 0.;

		// int neibcount;
		// #ifdef FIXED_NBSIZE
		// neibcount = neib_offs[i];
		// #else
		// neibcount =	neib_offs[i+1] - neib_offs[i];
		// #endif
		// printf("neibcount %d\n",neibcount);
		//printf("Nb indexed,i:%d\n",i);
	//for (int k=0;k < neibcount;k++) { //Or size
		// // //if fixed size i = part * NB + k
		// // //int j = neib[i][k];
		//int j = NEIB(i,k);
		//printf("i,j\n",i,j);
		// // double3 xij; 
		// // xij = x[i] - x[j];
		// // double h_ = (h[i] + h[j])/2.0;
		// // double nxij = length(xij);
		
		// // double GK	= GradKernel(3, 0, nxij/h_, h_);
		// // //		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
		// // //   mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
		// // dTdt[i] += m[j]/rho[j]*( 4.0*k_T[i]*k_T[j]/(k_T[i]+k_T[j]) * (T[i] - T[j])) * dot( xij , GK*xij )/(nxij*nxij);
	//}
	//dTdt[i] *=1/(rho[i]*cp[i]);
	//printf("dT: %f\n",dTdt[i]);
}

__global__ void TempCalcLeapfrogFirst(double *T, double *Ta, double *Tb, //output
																			double *dTdt, double dt){//input
	int i = threadIdx.x+blockDim.x*blockIdx.x;

	Ta[i] = T[i] - dt/2.0*dTdt[i];

}
__global__ void TempCalcLeapfrog     (double *T, double *Ta, double *Tb,
																		double *dTdt, double dt){
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	
	Tb[i]  = Ta[i];
	Ta[i] += dTdt[i] * dt;
	T [i] = ( Ta[i] + Tb[i] ) / 2.;
}
//Originally in Particle::TempCalcLeapfrog
//Host function
void Domain_d::ThermalSolve(const double &tf){
	int N = particle_count;
	int threadsPerBlock = 256;
	int blocksPerGrid =
	(N + threadsPerBlock - 1) / threadsPerBlock;
  Time =0.;
	//	while (Time<tf) {
	cout << "Callign Kernel"<<endl;
	
		ThermalSolveKernel<<<1,4>>>(dTdt,	
																		x, h, //Vector has some problems
																		m, rho, 
																		T, k_T, cp_T,
																		neib_part, neib_offs);
		cudaDeviceSynchronize(); //REQUIRED!!!!
		//cout << "Kernel called"<<endl;
		// if (isfirst_step) {
			// TempCalcLeapfrogFirst<<< 1,1 >>>(T, Ta, Tb,
																			 // dTdt, deltat);		
		// } else {
			// TempCalcLeapfrog <<< 1,1 >>>(T, Ta, Tb,
																			 // dTdt, deltat);				
		// }
		Time += deltat;
	//}//main time while
}//Thermal Solve

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
