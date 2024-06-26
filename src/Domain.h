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

#ifndef SPH_DOMAIN_H
#define SPH_DOMAIN_H

// #include <stdio.h>    // for NULL
// #include <algorithm>  // for min,max
#include <vector>
#include <iostream>

//#include <omp.h>

#include "Particle.h"
//#include "cuda/Functions.cuh"
#include "Boundary_Condition.h"
#include "cuda_runtime.h"
//#include "cuda/Domain_d.cuh"

// //#ifdef _WIN32 /* __unix__ is usually defined by compilers targeting Unix systems */
// //#include <sstream>
// //#endif
// //#include <sstream>
//#include <string>
// #include <cmath>

#include "Vector.h"
#define NONLOCK_SUM
#define MAX_NB_PER_PART 100

enum Kernels_Type { Qubic_Spline=0, Quintic=1, Quintic_Spline=2 ,Hyperbolic_Spline=3};
//C++ Enum used for easiness of coding in the input files
enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
enum Gradient_Type { Squared_density=0, Multiplied_density=1 };


namespace SPH {

class TriMesh;


int ComputeCylinderParticles( double Rxy, double Lz, double r);

class Domain
{

	//cuNSearch::NeighborhoodSearch neib;
	public:
	typedef void (*PtVel) (Vector & position, Vector & Vel, double & Den, Boundary & bdry);
	typedef void (*PtOut) (Particle * Particles, double & Prop1, double & Prop2,  double & Prop3);
	typedef void (*PtDom) (Domain & dom);
	// Constructor
	Domain();

	// Destructor
	~Domain();

  void AddCylUniformLength(int tag, double Rxy, double Lz, 
																				double r, double Density, double h, 
                                        double ang = 2.0*M_PI, int rows = 1, double r_i = 0.0);

	//Cylinder Slice 
	void AddXYSymCylinderLength(int tag, double Rxy, double Lz, 
								double r, double Density, double h, bool Fixed, bool symlength = false);
                
	// Domain Part
	void __host__ __device__ AddSingleParticle	(int tag, Vector const & x, double Mass, double Density, double h, bool Fixed);		//Add one particle
	void __host__ /*__device__ */ AddBoxLength				(int tag, Vector const &V, double Lx, double Ly, double Lz,double r, double Density,
																double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined dimensions

	void __host__ /*__device__*/ AddCylinderLength(int tag, Vector const & V, double Rxy, double Lz, 
								double r, double Density, double h, bool Fixed);

	void __host__ __device__ AddTractionProbeLength(int tag, Vector const & V, double Rxy, double Lz_side,
										double Lz_neckmin,double Lz_necktot,double Rxy_center,
										double r, double Density, double h, bool Fixed);
										
	void __host__ __device__ Calculate3DMass(double Density);
	void __host__ __device__ Add3DCubicBoxParticles(int tag, Vector const & V, double Lx, double Ly, double Lz, 
								double r, double Density, double h);


	void __host__ __device__ AddBoxNo						(int tag, Vector const &V, size_t nx, size_t ny, size_t nz,double r, double Density,
																double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined numbers
	void __host__ __device__ DelParticles				(int const & Tags);					//Delete particles by tag
	void __host__ __device__ CheckParticleLeave	();													//Check if any particles leave the domain, they will be deleted

	void YZPlaneCellsNeighbourSearch(int q1);						//Create pairs of particles in cells of XZ plan
	inline void MainNeighbourSearch				();									//Create pairs of particles in the whole domain
	inline void AllocateNbPair(const int &temp1, const int &temp2, const int &T);
	
	inline bool CheckRadius(Particle* P1, Particle *P2);
	
	void Move						(double dt);										//Move particles
	
	//This should be host mode
	void Solve					(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function

	//void Step(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);

	inline void CellInitiate		();															//Find the size of the domain as a cube, make cells and HOCs
	inline void ListGenerate		();															//Generate linked-list
	void CellReset			();															//Reset HOCs and particles' LL to initial value of -1

	void ClearNbData();	

	void WriteXDMF			(char const * FileKey);					//Save a XDMF file for the visualization
	__host__ __device__ inline void WriteCSV(char const * FileKey);

//	void InFlowBCLeave	();
//	void InFlowBCFresh	();
	void WholeVelocity	();

	void Kernel_Set									(Kernels_Type const & KT);
	void Viscosity_Eq_Set						(Viscosity_Eq_Type const & VQ);
	void Gradient_Approach_Set			(Gradient_Type const & GT);

	//Thermal Solver
	void CalcTempInc 	(); 		//LUCIANO: Temperature increment
	inline void CalcConvHeat ();
	inline void CalcPlasticWorkHeat();
	inline void CalcGradCorrMatrix();	//BONET GRADIENT CORRECTION
  
  inline void CalcPairPosList();
  void CheckParticlePairs();
  void InitReductionArraysOnce();
  
  void AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id);
  int AssignZone(Vector &start, Vector &end, int &id);


	// Data
	std::vector<Particle*>			Particles; 	///< Array of particles
	//Particle**			Particles; 	///< Array of particles
	double					R;		///< Particle Radius in addrandombox

	double					sqrt_h_a;				//Coefficient for determining Time Step based on acceleration (can be defined by user)

	int 					Dimension;    	///< Dimension of the problem

	double					MuMax;		///< Max Dynamic viscosity for calculating the timestep
	double					CsMax;		///< Max speed of sound for calculating the timestep
	double 				Vol;		///LUCIANO

	Vector					Gravity;       	///< Gravity acceleration


	Vector                 			TRPR;		///< Top right-hand point at rear of the domain as a cube
	Vector                  			BLPF;           ///< Bottom left-hand point at front of the domain as a cube
	Vector          CellSize;      	///< Calculated cell size according to (cell size >= 2h)
	int		          CellNo[3];      ///< No. of cells for linked list
	double 					hmax;						///< Max of h for the cell size  determination
	Vector          DomSize;				///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well

	double					rhomax;

	int						*** HOC;					///< Array of "Head of Chain" for each cell

	// BONET KERNEL CORRECTION
	bool 					gradKernelCorr;	

	double 					XSPH;		///< Velocity correction factor
	double 					InitialDist;	///< Initial distance of particles for Inflow BC

	double					AvgVelocity;	///< Average velocity of the last two column for x periodic constant velocity
	double 					getCellfac(){return Cellfac;}

	#ifdef __GNUC__
	size_t					Nproc;		///< No of threads which are going to use in parallel calculation
	#else
	int						Nproc;
	#endif
	//omp_lock_t 					dom_lock;	///< Open MP lock to lock Interactions array

	Boundary					BC;
	PtOut					UserOutput;
	PtVel 					InCon;
	PtVel 					OutCon;
	PtVel 					AllCon;
	Vector					DomMax;
	Vector					DomMin;
	PtDom					GeneralBefore;	///< Pointer to a function: to modify particles properties before CalcForce function
	PtDom					GeneralAfter;	///< Pointer to a function: to modify particles properties after CalcForce function
	size_t					Scheme;		///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog

	// Array<Array<std::pair<size_t,size_t> > >	SMPairs;
	std::vector< std:: vector<std::pair<int,int> > >  SMPairs;
	std::vector< std:: vector<std::pair<int,int> > > FSMPairs;
	std::vector< std:: vector<std::pair<int,int> > > NSMPairs;
	
	// Array<Array<std::pair<size_t,size_t> > >	NSMPairs;
	// Array<Array<std::pair<size_t,size_t> > >	FSMPairs;
	// Array< int > 				FixedParticles;
	
	 std::vector< int > 				FixedParticles;

	double 	& getTime (){return Time;}		//LUCIANO

	//Array<std::pair<size_t,size_t> >		Initial;
	std::vector< std::pair<int,int> >			Initial;

	//String					OutputName[3];
	double T_inf;			//LUCIANO: IN CASE OF ONLY ONE CONVECTION TEMPERAURE

	//iKernel m_kernel;
	bool					m_isNbDataCleared;
	bool						auto_ts;				//LUCIANO: Auto Time Stepping
  
  // CONTACT THINGS (ONLY FOR MESH GENERATION)
  int first_fem_particle_idx;

    //NEW: For parallel sum/reduction
    std::vector<size_t> first_pair_perproc;                   // Almost like pair count        
    int pair_count;                                                   //var names as stated as Nishimura (2015) ipl is njgi 
    std::vector < std::pair<int,int> >     pair_test,pair_ord;                    //OLY FOR TESTING
    //Array< Array <size_t> >               ilist_SM,jlist_SM;          // Size [Pairs] i and j particles of pair list [l], already flattened
    std::vector < size_t >                ipair_SM,jpair_SM;          //[Particles]// This is nb count for each particle i<j and j>i (called njgi) FLATTENED
    std::vector < size_t >                ipl_SM;                     // [Particles] position of link/pair (nb sum), called s_jgi in 1991 paper
    std::vector < std::vector <size_t>  > Aref;                      // Entry[Particle ,nb], indicates link
    std::vector < std::vector <size_t>  > Anei;                      //[Particles][MAX_NB_PER_PART] neighbiour list for j > i
    //std::vector <Vec3_t>                  pair_force;
    std::vector <double>                  temp_force;
    std::vector <double>                  pair_densinc;
    //std::vector <Mat3_t>                  pair_StrainRate;
    //std::vector <Mat3_t>                  pair_RotRate;    

private:
	void Periodic_X_Correction	(Vector & x, double const & h, Particle * P1, Particle * P2);		//Corrects xij for the periodic boundary condition
	void AdaptiveTimeStep				();		//Uses the minimum time step to smoothly vary the time step

	//void StartAcceleration					(Vector const & a = Vector(0.0,0.0,0.0));	//Add a fixed acceleration such as the Gravity
	//void PrimaryComputeAcceleration	();									//Compute the solid boundary properties
	//void LastComputeAcceleration		();									//Compute the acceleration due to the other particles

	inline __device__ void StartAcceleration					();
	inline __device__ void PrimaryComputeAcceleration	();									//Compute the solid boundary properties
    inline __device__ void LastComputeAcceleration		();									//Compute the acceleration due to the other particles
    inline __device__ void CalcForce2233	(Particle * P1, Particle * P2);		//Calculates the contact force between soil-soil/solid-solid particles

		
	void PrintInput			(char const * FileKey);		//Print out some initial parameters as a file
	void InitialChecks	();		//Checks some parameter before proceeding to the solution
	void TimestepCheck	();		//Checks the user time step with CFL approach

	size_t					VisEq;					//Choose viscosity Eq based on different SPH discretisation
	size_t					KernelType;			//Choose a kernel
	size_t					GradientType;		//Choose a Gradient approach 1/Rho i^2 + 1/Rho j^2 or 1/(Rho i * Rho j)
	double 					Cellfac;				//Define the compact support of a kernel

	double					Time;    				//Current time of simulation at each solving step
	double					deltat;					//Time Step
	double					deltatmin;			//Minimum Time Step
	double					deltatint;			//Initial Time Step

	std::vector<std::pair<size_t,size_t> > GhostPairs;	//If used  


};

__global__ void WriteCSV_kernel (Domain *d);
/*{
	d->WriteCSV("test");
}*/

}; // namespace SPH

// #include "Interaction.cpp"
//#include "Domain.cpp"
#include "NbSearch.cpp"
//#include "Output.cpp"
//#include "InOutFlow.cpp"
// #include "Thermal.cpp"

#endif // SPH_DOMAIN_H
