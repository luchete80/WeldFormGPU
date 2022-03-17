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
	
	cudaMalloc((void **)&h, 	particle_count * sizeof (double));
	cudaMalloc((void **)&m, 	particle_count * sizeof (double));
	cudaMalloc((void **)&rho, particle_count * sizeof (double));
	
	///////////
	//THERMAL //
	cudaMalloc((void **)&k_T, 	particle_count * sizeof (double));
	cudaMalloc((void **)&cp_T, 	particle_count * sizeof (double));
		
	cudaMalloc((void **)&T		, particle_count * sizeof (double));
	cudaMalloc((void **)&Ta		, particle_count * sizeof (double));
	cudaMalloc((void **)&Tb		, particle_count * sizeof (double));
	
	cudaMalloc((void **)&BC_T, particle_count * sizeof (int));
	
	//Host things
	T_h =  new double  [particle_count];
	x_h =  new double3 [particle_count];
	v_h =  new double3 [particle_count];
	u_h =  new double3 [particle_count];
	a_h =  new double3 [particle_count];
	
	cudaMalloc((void **)&dTdt	, particle_count * sizeof (double));
	//printf("Size of dTdt: %d, particle count %d\n",sizeof(dTdt)/sizeof (double),particle_count);

	//Nb data
	cudaMalloc((void **)&neib_offs	, (particle_count + 1) * sizeof (int));
	cudaMalloc((void **)&neib_part	, (particle_count * 400) * sizeof (int));
	
	cudaMalloc((void **)&partdata, sizeof(PartData_d));
	
	////////////////////////////
	/////// SPH ////////////////
	cudaMalloc((void **)&SumKernel, 	particle_count  * sizeof (double));	
	
	////////////////////////////
	//// MECHANICAL DATA ///////
	////////////////////////////
	cudaMalloc((void **)&v, particle_count * sizeof (double3));
	cudaMalloc((void **)&a, particle_count * sizeof (double3));
	cudaMalloc((void **)&u, particle_count * sizeof (double3));
	
	cudaMalloc((void **)&va, particle_count * sizeof (double3));
	cudaMalloc((void **)&vb, particle_count * sizeof (double3));
	
	cudaMalloc((void **)&p, 			particle_count * sizeof (double));	
	/// DensitySolid ///
	//DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
	cudaMalloc((void **)&PresEq, 	particle_count  * sizeof (double));	
	cudaMalloc((void **)&Cs, 			particle_count  * sizeof (double));		
	cudaMalloc((void **)&P0, 			particle_count  * sizeof (double));		
	cudaMalloc((void **)&FPMassC, particle_count  * sizeof (double));	
	cudaMalloc((void **)&rho_0, 	particle_count  * sizeof (double));	

	cudaMalloc((void **)&G, 	particle_count  * sizeof (double));	
	
	cudaMalloc((void **)&rhoa, 	particle_count  * sizeof (double));	
	cudaMalloc((void **)&rhob, 	particle_count  * sizeof (double));	
	cudaMalloc((void **)&drho, 	particle_count  * sizeof (double));		
	
	// STRESS AND STRAIN TENSORS - FLATTENED ARRAY!!!!
	cudaMalloc((void **)&sigma		, particle_count  * 6 * sizeof (double));		
	cudaMalloc((void **)&strrate	, particle_count  * 6 * sizeof (double));			
	cudaMalloc((void **)&rotrate	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
	
	cudaMalloc((void **)&shearstress	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
	cudaMalloc((void **)&shearstressa	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
	cudaMalloc((void **)&shearstressb	, particle_count  * 6 * sizeof (double));		 //ANTISYMM

	cudaMalloc((void **)&strain		, particle_count  * 6 * sizeof (double));		
	cudaMalloc((void **)&straina	, particle_count  * 6 * sizeof (double));		
	cudaMalloc((void **)&strainb	, particle_count  * 6 * sizeof (double));		
	
	// BOUNDARY CONDITIONS
	cudaMalloc((void **)&IsFree	, particle_count  * sizeof (bool));	
	cudaMalloc((void **)&NoSlip	, particle_count  * sizeof (bool));
	cudaMalloc((void **)&NSv, 		particle_count  * sizeof (double3));	
	cudaMalloc((void **)&ID, 			particle_count  * sizeof (int));	
	
	//////////////////////////
	/////// TENSILE INST /////
	cudaMalloc((void **)&TI, 					particle_count  * sizeof (double));	
	cudaMalloc((void **)&TIn, 				particle_count  * sizeof (double));		
	cudaMalloc((void **)&TIInitDist, 	particle_count  * sizeof (double));		
	cudaMalloc((void **)&TIR, 				6 * particle_count  * sizeof (double));	
	
	//////////////////////////
	/// CORRECTIONS /////////
	cudaMalloc((void **)&VXSPH, 	particle_count  * sizeof (double3));		
	//cudaMalloc((void **)&partdata->dTdt,particle_count * sizeof (double)); //TODO, pass to PartData
	
	// cudaMalloc((void**)&ppArray_a, 10 * sizeof(int*));
	// for(int i=0; i<10; i++) {
		// cudaMalloc(&someHostArray[i], 100*sizeof(int)); /* Replace 100 with the dimension that u want */
	// }
	
	Alpha= 1.;
	auto_ts = true;
	
	//To allocate Neighbours, it is best to use a equal sized double array in order to be allocated once
}


__host__ void Domain_d::SetFreePart(const Domain &dom){
	bool *k_ =  new bool[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = dom.Particles[i]->IsFree;
	}
	int size = particle_count * sizeof(bool);
	cudaMemcpy(this->IsFree, k_, size, cudaMemcpyHostToDevice);
	delete k_;	
}

__host__ void Domain_d::SetID(const Domain &dom){
	int *k_ =  new int[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = dom.Particles[i]->ID;
	}
	int size = particle_count * sizeof(int);
	cudaMemcpy(this->ID, k_, size, cudaMemcpyHostToDevice);
	delete k_;	
}

__host__ void Domain_d::SetCs(const Domain &dom){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = dom.Particles[i]->Cs;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->Cs, k_, size, cudaMemcpyHostToDevice);
	delete k_;	
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

void Domain_d::SetShearModulus(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->G, k_, size, cudaMemcpyHostToDevice);
	delete k_;
}


void Domain_d::SetDensity(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->rho, k_, size, cudaMemcpyHostToDevice);
	cudaMemcpy(this->rho_0, k_, size, cudaMemcpyHostToDevice);
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


#include <cstdio>

void Domain_d::WriteCSV(char const * FileKey){
	FILE *f = fopen(FileKey,"w");;
	
	fprintf(f, "X, Y, Z\n");

	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move

	for (int i=0; i<particle_count; i++) {
		fprintf(f,"%f, %f, %f\n",a_h[i].x,a_h[i].y,a_h[i].z);
		//Particles[i]->CalculateEquivalentStress();		//If XML output is active this is calculated twice
		//oss << Particles[i]->Sigma_eq<< ", "<< Particles[i]->pl_strain <<endl;
	}

 fclose(f);
}

//TODO: CHANGE TO MECH
__device__ inline void Domain_d::AdaptiveTimeStep(){
		// //Min time step check based on the acceleration
		// double test	= 0.0;
		// deltatmin	= deltatint;

		// if (IsFree[i]) {
			// test = sqrt(h[i]/norm(Particles[i]->a));
			// if (deltatmin > (sqrt_h_a*test))
			// {
				// omp_set_lock(&dom_lock);
					// deltatmin = sqrt_h_a*test;
				// omp_unset_lock(&dom_lock);
			// }
		// }
		
}

};//SPH