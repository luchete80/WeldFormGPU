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
//C++ Enum used for easiness of coding in the input files

//enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
//enum Gradient_Type { Squared_density=0, Multiplied_density=1 };

//#include <cuNSearch.h>



namespace SPH {

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
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neib_part;	//1D array, faster
	int *neib_offs;	//Offset or count
	
	int *neib_count; //Optional
	int *neibcount;	//Useful??
	int particle_count;
	
	//SPH
	double *h;
	double *SumKernel;
	
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
	
	bool						auto_ts;

	
	
	
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
	
	int 			*ID;

		/// FUNCTIONS
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
	__device__ /*__forceinline__*/inline void CalcForce2233(int KernelType, float XSPH);
	__device__ void StressStrain(int i);
	
	__device__ void ApplyBCVel(int bcid, 
														double3 bcv);
	__host__ void WriteCSV(char const * FileKey);
	__device__ inline void AdaptiveTimeStep();
 ////////////////////////
	


	__device__ void WholeVelocity();

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

__global__ void CalcForcesKernelMember(PartData_d *partdata);

void __global__ /*inline*/ CalcForcesKernel(
																		double3 *a, double *drho,				//OUTPUT
																		double3 *x, double *h, double3* v,
																		double *m, double *rho, double *FPMassC,
																		double *Cs, double *P0,double *p, double *rho_0,
																		int *neib_part, int *neib_offs,
																		double *sigma,
																		double *strrate, double *rotrate,
																		int particle_count);
																		
__global__ void CalcForcesKernel(Domain_d *dom_d);

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

	/* const double &Dimension*/

// /*inline*/ __host__ void StartAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void StartAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void PrimaryComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void PrimaryComputeAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void LastComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void LastComputeAccelerationKernel(Domain_d &sd);

}; // namespace SPH

#endif // SPH_DOMAIN_H
