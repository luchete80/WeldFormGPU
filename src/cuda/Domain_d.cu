#include "Domain_d.cuh"
#include "Functions.cuh"
#include "Domain.h"

#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

#include "Mesh.h"
#include "cudautils.cuh"
#include "Mesh.cuh"

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
  if (thermal_solver){
    cudaMalloc((void **)&k_T, 	particle_count * sizeof (double));
    cudaMalloc((void **)&cp_T, 	particle_count * sizeof (double));
      
    cudaMalloc((void **)&T		, particle_count * sizeof (double));
    cudaMalloc((void **)&Ta		, particle_count * sizeof (double));
    cudaMalloc((void **)&Tb		, particle_count * sizeof (double));
    
    cudaMalloc((void **)&BC_T, particle_count * sizeof (int));


    // cudaMalloc((void **)&T_inf , particle_count * sizeof (double));
    // cudaMalloc((void **)&q_conv, particle_count * sizeof (double));
    // cudaMalloc((void **)&h_conv, particle_count * sizeof (double));	
  }  
	//Host things
	T_h =  new double  [particle_count];
	x_h =  new double3 [particle_count];
	v_h =  new double3 [particle_count];
	u_h =  new double3 [particle_count];
	a_h =  new double3 [particle_count];
  

  
  normal_h = new double3 [particle_count];
  nb_h = new unsigned int [particle_count];
	
  ID_h =  new int    [particle_count];
	
	sigma_eq_h =  new double [particle_count];
	pl_strain_h = new double [particle_count];
	
	p_h =  new double [particle_count];
	rho_h =  new double [particle_count];
	
	max_deltat_h = new double [particle_count];
	
	cudaMalloc((void **)&dTdt	, particle_count * sizeof (double));
	//printf("Size of dTdt: %d, particle count %d\n",sizeof(dTdt)/sizeof (double),particle_count);
  
  //Contact 
  contneib_count_h = new int [particle_count];

	//Nb data
	cudaMalloc((void **)&neib_offs	, (particle_count + 1) * sizeof (int));
	cudaMalloc((void **)&neib_part	, (particle_count * 100) * sizeof (int));
	
	cudaMalloc((void **)&partdata, sizeof(PartData_d));
  
  cudaMalloc((void **)&test, sizeof(particle_count)* sizeof (double));
	
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
  
  report_gpu_mem_();
  
	/// DensitySolid ///
	//DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
  cout << "Allocating  state vars "<<endl;
	cudaMalloc((void **)&PresEq, 	particle_count  * sizeof (double));	
	cudaMalloc((void **)&Cs, 			particle_count  * sizeof (double));		
	cudaMalloc((void **)&P0, 			particle_count  * sizeof (double));		
	cudaMalloc((void **)&FPMassC, particle_count  * sizeof (double));	
	cudaMalloc((void **)&rho_0, 	particle_count  * sizeof (double));	

	cudaMalloc((void **)&G, 	particle_count  * sizeof (double));	
	
	// cudaMalloc((void **)&rhoa, 	particle_count  * sizeof (double));	
	// cudaMalloc((void **)&rhob, 	particle_count  * sizeof (double));	
	cudaMalloc((void **)&drho, 	particle_count  * sizeof (double));		
	
	// STRESS AND STRAIN TENSORS - FLATTENED ARRAY!!!!
	cudaMalloc((void **)&sigma		, particle_count  * 6 * sizeof (double));		
	cudaMalloc((void **)&strrate	, particle_count  * 6 * sizeof (double));			
	cudaMalloc((void **)&rotrate	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
  
	
	cudaMalloc((void **)&shearstress	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
	// cudaMalloc((void **)&shearstressa	, particle_count  * 6 * sizeof (double));		 //ANTISYMM
	// cudaMalloc((void **)&shearstressb	, particle_count  * 6 * sizeof (double));		 //ANTISYMM

	cudaMalloc((void **)&sigma_eq		, particle_count  * sizeof (double));		
	cudaMalloc((void **)&pl_strain	, particle_count  * sizeof (double));	
	cudaMalloc((void **)&sigma_y		, particle_count  * sizeof (double));		
	cudaMalloc((void **)&Et		      , particle_count  * sizeof (double));		
	
	cudaMalloc((void **)&strain		, particle_count  * 6 * sizeof (double));		
	// cudaMalloc((void **)&straina	, particle_count  * 6 * sizeof (double));		
	// cudaMalloc((void **)&strainb	, particle_count  * 6 * sizeof (double));		
  
  cudaMalloc((void **)&eff_strain_rate	, particle_count  * sizeof (double));		
	
	// BOUNDARY CONDITIONS
	cudaMalloc((void **)&IsFree	, particle_count  * sizeof (bool));	
	cudaMalloc((void **)&NoSlip	, particle_count  * sizeof (bool));
	cudaMalloc((void **)&NSv, 		particle_count  * sizeof (double3));	
	cudaMalloc((void **)&ID, 			particle_count  * sizeof (int));	
	
  
  

  cudaMalloc((void**)&mat, particle_count * sizeof(Material_ *));

	//////////////////////////
	/////// TENSILE INST /////
	cudaMalloc((void **)&TI, 					particle_count  * sizeof (double));	
	cudaMalloc((void **)&TIn, 				particle_count  * sizeof (double));		
	cudaMalloc((void **)&TIInitDist, 	particle_count  * sizeof (double));		
	cudaMalloc((void **)&TIR, 				6 * particle_count  * sizeof (double));	

	cudaMalloc((void **)&max_deltat, 	particle_count  * sizeof (double));		
  
  
  /////////////CONTACT////////////////
	cudaMalloc((void **)&normal,          particle_count * sizeof (double3));	
  
  cudaMalloc((void **)&contneib_count,  particle_count * sizeof (int));
  cudaMalloc((void **)&contneib_part,   particle_count * 20 * sizeof (int));
  cudaMalloc((void **)&contneib_offs,   particle_count * sizeof (int));
  cudaMalloc((void **)&contforce,       particle_count * sizeof (double3));
  
  contforce_h = new double3 [particle_count];
  
  //////////////////////////
  ///// ENERGY /////////////
  cudaMalloc((void **)&int_energy,   particle_count * sizeof (double));
  cudaMalloc((void **)&kin_energy,   particle_count * sizeof (double));
	
	//////////////////////////
	/// CORRECTIONS /////////
	cudaMalloc((void **)&VXSPH, 	particle_count  * sizeof (double3));		
  
  
  //pplane_h =  new double [particle_count-solid_part_count]; ////TEST

  
  trimesh = NULL;
	//cudaMalloc((void **)&partdata->dTdt,particle_count * sizeof (double)); //TODO, pass to PartData
	
	// cudaMalloc((void**)&ppArray_a, 10 * sizeof(int*));
	// for(int i=0; i<10; i++) {
		// cudaMalloc(&someHostArray[i], 100*sizeof(int)); /* Replace 100 with the dimension that u want */
	// }
	
	Alpha = 1.;
  Beta  = 0.;
	auto_ts = true;
	
	deltat	= 0.0;
	deltatint	= 0.0;
	deltatmin	= 0.0;
	sqrt_h_a = 0.0025;	
	
	//Initiate pl_strain, is it necessary
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = 0.;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->pl_strain, k_, size, cudaMemcpyHostToDevice);
	delete k_;
  
  //Change to first particle
  cudaMalloc((void **)&not_write_surf_ID, 	particle_count  * sizeof (bool));		
  bool *surf = new bool[particle_count];
  for (int i=0;i<particle_count;i++){
    surf[i] = false;
  }
  cudaMemcpy(not_write_surf_ID, surf,  particle_count * sizeof(bool), cudaMemcpyHostToDevice);

	
	// CONTACT THINGS
	contact_force_factor =1.;
	PFAC =0.8;
	DFAC =0.2;
  contact = false;
  
  //SetDouble(this->int_energy,0.);
	double *ki_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		ki_[i] = 0.0;
	}

	cudaMemcpy(int_energy, ki_, particle_count * sizeof(double), cudaMemcpyHostToDevice);
	delete ki_;	
  
  solid_part_count = particle_count;
	
	//To allocate Neighbours, it is best to use a equal sized double array in order to be allocated once
  friction_sta = friction_dyn = 0.;
  cout << "Domain allocation done. "<<endl;
  report_gpu_mem_();
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
  cout << "Copying "<<particle_count<<" particles id"<<endl;
	int size = particle_count * sizeof(int);
	cudaMemcpy(this->ID, k_, size, cudaMemcpyHostToDevice);
	delete k_;	
}

// NOW ALLOWED ACCORDING TO CUDA 12
// __global__
// void init_hollomon_mat_kernel(Material_ ** p)
// {
  // // if (threadIdx.x == 0 && blockIdx.x == 0){
  // // Elastic_ el;
  // // (*p) = new Hollomon(el,1.0,1.0,1.0);
  // // }
// }

__global__
void init_hollomon_mat_kernel(Domain_d *dom)
{
  if (threadIdx.x == 0 && blockIdx.x == 0){
  dom->init_hollomon_material();
  }
  // (*p) = new Hollomon(el,1.0,1.0,1.0);
  // }
}

// void   __device__ Domain_d::init_hollomon_material(){
  // materials = new Hollomon();
// }

void   __device__ Domain_d::init_hollomon_material(){
  *materials_ptr = new Hollomon();
  materials = new Hollomon();
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

//TODO: MOVE TO TEMPLATE
__host__ void Domain_d::SetDouble(double *arr, double val){
	int *k_ =  new int[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = val;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(arr, k_, size, cudaMemcpyHostToDevice);
	delete k_;	
}

// template <typename T>
// __host__ void Domain_d::SetType(T* arr, double val){
	// T *k_ =  new T[particle_count];
	// for (int i=0;i<particle_count;i++){
		// k_[i] = val;
	// }
	// int size = particle_count * sizeof(T);
	// cudaMemcpy(arr, k_, size, cudaMemcpyHostToDevice);
	// delete k_;	
// }

// template __host__ void Domain_d::SetType<double>(double *arr, double val);

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

__device__ void Domain_d::AssignMatAddress(int i){
  if (i< solid_part_count)
    mat[i] = &materials[0];
  
}

__global__ void AssignMatAddressKernel(Domain_d *dom){
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  dom->AssignMatAddress(i);
}

// __global__ void InitMatHollomonKernel(Domain_d *dom, const Material_ &mat_h){
  // int i = threadIdx.x + blockDim.x*blockIdx.x;
  // if (i< dom->solid_part_count)
    // dom->InitMatHollomon(mat_h);
// }

// __device__ void Domain_d::InitMatHollomon(const Material_ &mat){

  // materials[0].K = 0.0;
// }


void Domain_d::Set_h(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->h, k_, size, cudaMemcpyHostToDevice);
	h_glob = k;
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

void Domain_d::SetSigmay(const double &k){
	double *k_ =  new double[particle_count];
	for (int i=0;i<particle_count;i++){
		k_[i] = k;
	}
	int size = particle_count * sizeof(double);
	cudaMemcpy(this->sigma_y, k_, size, cudaMemcpyHostToDevice);
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
	
	fprintf(f, "ID, X, Y, Z, posX,posY,posZ,Ux, Uy, Uz, Vx, Vy, Vz, Ax, Ay, Az, rho, p, SigmaEq, Pl_Strain, Nb, ContNb, CFx, CFy, CFz, Nx, Ny, Nz\n");

	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
  double ppl;
	for (int i=0; i<particle_count; i++) {
		fprintf(f,"%d, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %d, %d, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e\n", 
							ID_h[i],
              x_h[i].x,x_h[i].y,x_h[i].z, 
              x_h[i].x,x_h[i].y,x_h[i].z, 
              u_h[i].x,u_h[i].y,u_h[i].z,
							v_h[i].x,v_h[i].y,v_h[i].z, 
							a_h[i].x,a_h[i].y,a_h[i].z,
							//ID[i],
						rho_h[i],
						p_h[i],
						sigma_eq_h[i],
						pl_strain_h[i],
            nb_h[i],
            contneib_count_h[i],
            contforce_h[i].x,contforce_h[i].y,contforce_h[i].z,
            normal_h[i].x,normal_h[i].y,normal_h[i].z
            );
		//Particles[i]->CalculateEquivalentStress();		//If XML output is active this is calculated twice
		//oss << Particles[i]->Sigma_eq<< ", "<< Particles[i]->pl_strain <<endl;
	}

 fclose(f);
}

__global__ void CalcMinTimeStepKernel(Domain_d *dom){
	
	dom->CalcMinTimeStep();	//Stablish deltatmin based on acceleration
	
}

//TODO: CHANGE TO MECH
__device__ void Domain_d::CalcMinTimeStep(){
		int i = threadIdx.x + blockDim.x*blockIdx.x;
		//THIS WAS IN LASTCOMPUTEACCELERATION original code
		// //Min time step check based on the acceleration
		if (i<particle_count){
			double test	= 0.0;
			deltatmin	= deltatint;
			double sqrt_h_a = 0.0025;
			//Appears to be safe
			//https://stackoverflow.com/questions/8416374/several-threads-writing-the-same-value-in-the-same-global-memory-location
	
			//ORIGINAL: test = sqrt(h[i]/length(a[i]));
			if (IsFree[i]) {
				
				//test = h[i]/(Cs[i]+length(v[i]));
				//max_deltat[i]=0.6*test;
				//if (deltatmin > (0.6*test)) 
					//deltatmin = 0.6*test;
        
				test = sqrt(h[i]/length(a[i]));
				max_deltat[i]=sqrt_h_a*test;				
				if (deltatmin > (sqrt_h_a*test)) 
						deltatmin = sqrt_h_a*test;


						
						//printf("particle i: %d Min time step %f\n",i,deltatmin);
				
			}
		}		
}

__host__ void Domain_d::AdaptiveTimeStep(){
		if (deltatint>deltatmin) {
		if (deltat<deltatmin)
			deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		else
			deltat		= deltatmin;
	} else {
		if (deltatint!=deltat)
			deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		else
			deltat		= deltatint;
	}
}

//If this is called, all particle has to be reallocated
__host__ void Domain_d::AddTrimeshParticles(TriMesh_d &mesh, const double &hfac, const int &id){
	
	first_fem_particle_idx = particle_count;
	double Density =0.;
	double h;
	bool Fixed = false;	//Always are fixed ...
	contact_surf_id = id;
	trimesh = &mesh;
	//TODO: MAKE SEVERAL ONE PER MESH
  //cudaMalloc((void **)&element,       mesh.element.size() * sizeof (int));
  
	//for ( int e = 0; e < mesh.element.size(); e++ ){
		// Vector pos = mesh.element[e]->centroid;
		// h = hfac * mesh.element[e]->radius;
		// Particles.Push(new Particle(id,pos,Vector(0,0,0),0.0,Density,h,Fixed));
		// Particles[first_fem_particle_idx + e] -> normal  = mesh.element[e] -> normal;
		// Particles[first_fem_particle_idx + e] -> element = e; 
	//}
	// cout << Particles.size() - first_fem_particle_idx << "particles added with ID " << contact_surf_id <<endl;
	// cout << first_fem_particle_idx << " is the first solid particle index."<<endl;
}


// THIS SHOULD BE DONE
	// if (deltatint>deltatmin)
	// {
		// if (deltat<deltatmin)
			// deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		// else
			// deltat		= deltatmin;
	// }
	// else
	// {
		// if (deltatint!=deltat)
			// deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		// else
			// deltat		= deltatint;
	// }
	
	// if (contact){
		// if (min_force_ts < deltat)
		// //cout << "Step size changed minimum Contact Forcess time: " << 	min_force_ts<<endl;
		// deltat = min_force_ts;
	// }

	// if (deltat<(deltatint/1.0e5))
		// //cout << "WARNING: Too small time step, please choose a smaller time step initially to make the simulation more stable"<<endl;
		// throw new Fatal("Too small time step, please choose a smaller time step initially to make the simulation more stable");
// }


};//SPH