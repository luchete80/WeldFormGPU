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

//C++ Enum used for easiness of coding in the input files

//enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
//enum Gradient_Type { Squared_density=0, Multiplied_density=1 };

//#include <cuNSearch.h>

//#define neibs (i,j)	neibs[i*neibcount[i-1]+j]

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

//#ifdef FIXED_NBSIZE //fixed nb per part
#define MAXNB_PPART
//#define NEIB(i, k) neib [ MAXNB_PPART * i + k]  
//#else
//#define NEIB(i, k) neib [ neibcount[i-1]  i + k]  
//#endif

class Domain;

class Domain_d
{
	public:
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neib_idx;	//1D array, faster
	int *neib_part;	//1D array, faster
	
	int *neib_count; //Optional
	int *neibcount;	//Useful??
	
	double *h;
	double3* x; //Vector is double
	Vector* v;
	Vector* a;
	
	//Time things
	bool isfirst_step;

	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step

	
	double *rho, *m;	//Mass and density
	//THERMAL
	double *T, *Ta, *Tb, *dTdt;
	double *k_T, *cp_T,*h_conv, *T_inf;
			
	Domain_d(){};
	Domain_d(const int &particle_count);
	__host__ void SetDimension(const int &particle_count);//Called from kernel to assign with CUDA_MALLOC
	__host__ void ThermalSolve(const double &tf);
	~Domain_d();
	
	__host__ void Domain_d::CopyData(const Domain &dom);

};

PartData_d{
	
	public:
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neibs;	//1D array, faster
	int *neib_offs;	//1D array, faster
	int *neibcount;	//Useful??
	
	//General
	double *rho, *m;	//Mass and density
	//THERMAL
	double *T, *Ta, *Tb, *dTdt;
	double *k_T, *cp_T,*h_conv, *T_inf;
	
};

//Called by Solve host function
__global__ void ThermalSolveKernel(double *dTdt, 
																		double3 *x, double *h,
																		double *mass, double *rho, 
																		double *T, double *k_T, double *cp_T, 
																		int *neib, int *neibcount); //Idea is to pass minimum data as possible


__global__ void TempCalcLeapfrogFirst(double *T,double *Ta, double *Tb, 
																			double *dTdt, double dt);
__global__ void TempCalcLeapfrog     (double *T, double *Ta, double *Tb, 
																			double *dTdt, double dt);
	
	

// /*inline*/ __host__ void StartAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void StartAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void PrimaryComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void PrimaryComputeAccelerationKernel(Domain_d &sd);

// /*inline*/ __host__ void LastComputeAcceleration(Domain_d &sd); // This is the buffer function which calls the kernel
// __global__ void LastComputeAccelerationKernel(Domain_d &sd);

}; // namespace SPH

#endif // SPH_DOMAIN_H
