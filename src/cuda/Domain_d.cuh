/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *
*             and soils) using Smoothed Particle Hydrodynamics method              *
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#ifndef SPH_Domain_d_CUH
#define SPH_Domain_d_CUH

// #include <stdio.h>    // for NULL
// #include <algorithm>  // for min,max


// #include <omp.h>

// #include "Particle_d.cuh"
// #include "Functions.h"
#include "tensor.cuh"
//#include "Boundary_Condition.h"

// //#ifdef _WIN32 /* __unix__ is usually defined by compilers targeting Unix systems */
// #include <sstream>
// //#endif
// #include <sstream>
// #include <string>
// #include <cmath>
// #include "tensor.h"
#include "Vector.h"
#include "vector_math.h"
#include "PartData.cuh"

#include "cuNSearch.h"
#include "cuNSearchDeviceData.h"


//C++ Enum used for easiness of coding in the input files

//enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
//enum Gradient_Type { Squared_density=0, Multiplied_density=1 };

//#include <cuNSearch.h>

#define DEBUG_MODE

namespace SPH {

class TriMesh;
class TriMesh_d;
class Boundary;

/******************************************/
/* CELL STRUCT LEADING TO ARRAY OF STRUCT */
/******************************************/
// struct cellAoS {

    // unsigned int    x1;
    // unsigned int    x2;
    // unsigned int    code;
    // bool            done;

// };

// /*******************************************/
// /* CELL STRUCT LEADING TO STRUCT OF ARRAYS */
// /*******************************************/
// struct cellSoA {

    // unsigned int    *x1;
    // unsigned int    *x2;
    // unsigned int    *code;
    // bool            *done;

// };

class Domain;

// LUCIANO
// DOMAIN_D is basically a SOA.
// In the future it is of interest to compare passing this in kernel vs
// passing less arguments (oir by grouping arrays per Mech, Therm, Stress, etc.)
// in order to check if computational speed is improved
class Domain_d
{
	public:
	///////////////////////////////////
	///////////// NEIGHBOR THING //////
	
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neib_part;	//1D array, faster
	int *neib_offs;	//Offset or count
	
	int *neib_count; //Optional
	int *neibcount;	//Useful??
	int particle_count;
	
	
	//cuNSearch::PointSet part_pointset; 	// IN THE FUTURE IS GOOD TO HAVE NEIGHBOUR DEVICE DATA WHICH IS IN DEVICE	
														// INSTEAD OF HOST
														
	cuNSearch::NeighborhoodSearch nb_search;
														
	cuNSearch::cuNSearchDeviceData nb_device_data;
	
	//SPH
	double *h;
	double *SumKernel;
	
	double h_glob;	//Initial h
	
	double3* x; //Vector is double
	double3* v;
	double3* a;
	double3* u;
	
	double3 *a_h, *x_h, *v_h, *u_h;
	
	PartData_d *partdata;
  	
	//Time things
	bool isfirst_step;
	int step = 0;

	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step
	double 					sqrt_h_a;
	
	bool						auto_ts;
	
	double 					*max_deltat, *max_deltat_h;		//According to each particle


	
	/// TODO: PASS THIS TO PARTICLE DATA
	double *rho, *m;	//Mass and density
	//THERMAL
	double *T, *Ta, *Tb, *dTdt;
	double *T_h;	//host (for finding max, writing)
	
	double *k_T, *cp_T,*h_conv;
	
	int *BC_T;	//0 nothing, 1 convection
	//double *T_inf;

	//Mechanical
	double *p,*PresEq;
	double *Cs, *P0;
	
	//Material params elastic
	double *G;
	
	double *rho_0;	///< Reference Density of Particle	
	
	double *rhoa,*rhob,*drho;
	double3 *va,*vb;
	//
	
	//STRESS AND STRAIN TENSORS
	double *sigma; //To convert after to tensor;
	
	double *strrate,*rotrate;//all flattened, six component (rotation rate is upped matrix component, since it is antisymm)
	double *shearstress,*shearstressa,*shearstressb;
	double *strain,*straina,*strainb;
	
	double *sigma_eq, *sigma_eq_h;
	double *pl_strain,*pl_strain_h;
	double *p_h,*rho_h;
	
	double *sigma_y,*sigma_y_h;
	//Be in another class
	double  *FPMassC;        ///< Mass coefficient for fixed particles to avoid leaving particles
		
	//ARTIFFICIAL VISC	
	double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
	double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
	
	//TENSILE INSTABILITY
	double	*TI;		///< Tensile instability factor
	double	*TIR;		///< Tensile Instability stress tensor R
	double	*TIn;		///< Tensile instability power
	double 	*TIInitDist;	///< Initial distance of particles for calculation of tensile instability
	
	double3 *VXSPH;
	
	//BOUNDARY
	bool 			*IsFree, *NoSlip;
	double3 	*NSv;	///< Velocity of the fixed particle for no-slip BC
	
	int 			*ID, *ID_h;
	
	
	////////////////////////////////////
	/////// CONTACT THINGS /////////////
	////////////////////////////////////
  int id_free_surf, contact_surf_id;
  bool contact;
  double totmass;
	double 	max_contact_force;
	double3 *normal; 
	int 		  **contneib;	//array of lists, NOT IN USE
	
  int 		  *contneib_part;	//1D array, faster, THESE ARE THE NEIGHBOURS ITSELF (OF THE SURFACE CONTACT)
	int 		  *contneib_offs;	//Offset or count
  int       *contneib_count; //REDUNDANCE WITH contneib_offs
  
	double 	  *cont_stiff;
	double3   *contforce;	//SOA
  int       first_fem_particle_idx;
	double    contact_force_factor;
	double    PFAC,DFAC;
	double    fritcion_sta,fritcion_dyn;
  int       *contneib_count_h;
  double min_force_ts;
  double friction_dyn;

  

	// TODO, EACH RIGID PARTICLE SHOULD 
  int 		*element; //ELEMENT OF TRIMESH FROM "RIGID" PARTICLE, ALL FIRST PARTICLES ARE ZERO
  TriMesh_d *trimesh;
	
	/////////////////////////////////////////
	///////// MEMBER FUNCTIONS /////////////
	Domain_d(){isfirst_step=true;};
	Domain_d(const int &particle_count);
	__host__ void SetDimension(const int &particle_count);//Called from kernel to assign with CUDA_MALLOC
	__host__ void Set_h(const double &);
	__host__ void ThermalSolve(const double &tf);

	__host__ void MechSolve(const double &tf, const double &dt_out);
	
	//General
	__host__ void SetDensity (const double &k);
	//__host__ void SetDensity0(const double &k);
	__host__ void SetConductivity(const double &k);
	__host__ void SetSigmay(const double &k);
	__host__ void SetHeatCap(const double &);
	
	//Boundary
	__host__ void SetFreePart(const Domain &dom);
	__host__ void SetShearModulus(const double &);
	__host__ void SetID(const Domain &dom);
	
	//Mechanical
	__host__ void SetCs(const Domain &dom);
	~Domain_d();
	
	__host__ void Domain_d::CopyData(const Domain &dom);
	__device__ void CheckData();
	__device__ void CalcThermalTimeStep();

	//MAIN MECHANICAL FUNCTIONS ARE THESE THREE (Start Acc & Whole Velocity are not)
	__device__ __forceinline__ void PrimaryComputeAcceleration ();	
	__device__ __forceinline__ void LastComputeAcceleration();
	__device__ /*__forceinline__*/inline void CalcForce2233(const uint *particlenbcount,
																													const uint *neighborWriteOffsets,
																													const uint *neighbors,
																													int KernelType, float XSPH);

	__device__ /*__forceinline__*/inline void CalculateSurface(const uint *particlenbcount,
																													const uint *neighborWriteOffsets,
																													const uint *neighbors,

																													/*const int &id, */const double &totmass);
  
  //CREATE NB LIST FROM ORIGINAL CONTACT NB LIST
  __device__ inline void CalcContactNb(const uint *particlenbcount,
                            const uint *neighborWriteOffsets,
                            const uint *neighbors);

	__device__ inline void CalcContactForces(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors);


	__device__ void StressStrain(int i);
	
	__device__ void ApplyBCVel(int bcid, 
														double3 bcv);
	__host__ void WriteCSV(char const * FileKey);
	__device__ void CalcMinTimeStep();
	__host__ void AdaptiveTimeStep();
 ////////////////////////
	


	__device__ void WholeVelocity();
  
  //Contact things 
  void __device__ inline CalcContactForces(/* int i*/);
  __host__ void AddTrimeshParticles(TriMesh_d &mesh, const float &hfac, const int &id);
    
}; 



__global__ void CheckData(Domain_d *dom);

//Called by Solve host function
//TODO: REMOVE; PASS ENTIRE CLASS
__global__ void ThermalSolveKernel(double *dTdt, 
																		double3 *x, double *h,
																		double *mass, double *rho, 
																		double *T, double *k_T, double *cp_T, 
																		int *neib_part, int *neib_offs,
																		int count); //Idea is to pass minimum data as possible

//NEXT SOLVER
void __global__ ThermalSolveKernel (double dt, PartData_d *partdata);

// void __global__ CalcConvHeatKernel (double *dTdt,
																		// double *m, double *rho, double *cp_T,
																		// double *T, double T_inf,
																		// int *BC_T,
																		// double &h_conv);

__global__ void TempCalcLeapfrogFirst(double *T,double *Ta, double *Tb, 
																			double *dTdt, double dt,
																			int count);
__global__ void TempCalcLeapfrog     (double *T, double *Ta, double *Tb, 
																			double *dTdt, double dt,
																			int count);
	

void __global__ MechSolveKernel (double dt, PartData_d *partdata);

__device__ void CalcForcesExt(PartData_d *partdata);

__global__ inline void CalcForcesKernelMember(PartData_d *partdata);

                                  
__global__ void WholeVelocityKernel(Domain_d *dom_d);

void __global__ CalcConvHeatKernel (double *dTdt,
																		double *m, double *rho, double *cp_T,
																		double *T, double T_inf,
																		int *BC_T,
																		double h_conv, int count);
																		
void __global__ MoveKernelExt(double3 *v, double3 *va, double3 *vb,
													double *rho, double *rhoa, double *rhob, double *drho,
													double3 *x, double3 *a,
													double3 *u, /*Mat3_t I, */double dt,
													bool , int count);
													
void __global__ MoveKernelDom(Domain_d *dom);

__global__ void StressStrainExtKernel(double *sigma,	//OUTPUT
																								double *strain,double *straina,double *strainb, //OUTPUT
																								//INPUT
																								double *p, double *rotrate, 
																								double *shearstress,double *shearstressa, double *shearstressb,
																								
																								double dt, int particle_count);

__global__ void PressureKernelExt(double *p,
																	double *PresEq, double *Cs, double *P0,double *Density, double *RefDensity, int particle_count);																		

__global__ void ApplyBCVelExtKernel(double *v,
																double *va,
																int *ID, 	//Input
																int bcid, 
																double bcv,
																double Time,
																int particle_count);

__global__ void ApplyBCVelKernel (Domain_d *dom, int bcid, double3 bcv);

__global__ void StressStrainKernel(Domain_d *dom);

__global__ void CalcMinTimeStepKernel(Domain_d *dom);

__global__ void TimestepCheckKernel(const double &CFL,
																double *h,
																double *Cs);

	/* const double &Dimension*/

// /*inline*/ __host__ void StartAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void StartAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void PrimaryComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void PrimaryComputeAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void LastComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void LastComputeAccelerationKernel(Domain_d &sd);

}; // namespace SPH

#include "Mesh.h"

#endif // SPH_DOMAIN_H
