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

#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H

#include "Vector.h"
//#include "Functions.h"

#define TH_BC_NONE			0
#define TH_BC_CONVECTION	1
#define TH_BC_CONDUCTION	2

enum Ghost_Type {Symmetric = 0, Periodic = 1, Mirror_XYZ = 2 };

namespace SPH {

	class Particle
	{
	public:
		// Shepard density correction
		bool   	Shepard;	///< Shepard Filter for the density
		size_t	ShepardCounter;	///< Count number of contributing particles
		size_t	ShepardStep;	///< Cycle number for shepard correction
		double	ZWab;		///< Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
		double	SumDen;		///< Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)

		bool   	IsFree;		///< Check the particle if it is free to move or not
		size_t	InOut;		///< Check the particle if it is in-flow or out-flow or not
		bool   	IsSat;		///< Check the particle if it is Saturated or not
		bool   	SatCheck;	///< Check the particle Saturation at each time step
		bool   	NoSlip;		///< No-Slip BC

		int    	ID;		///< an Integer value to identify the particle set
		int 	Thermal_BC;
		int    	Material;	///< an Integer value to identify the particle material type: 1 = Fluid, 2 = Solid, 3 = Soil

		Vector	x;		///< Position of the particle n
		Vector	vb;		///< Velocity of the particle n-1 (Modified Verlet)
		Vector	va;		///< Velocity of the particle n+1/2 (Leapfrog)
		Vector	v;		///< Velocity of the particle n+1
		Vector	NSv;		///< Velocity of the fixed particle for no-slip BC
		Vector	VXSPH;		///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
		Vector	a;		///< Acceleration of the particle n

    bool not_write_surf_ID;

		size_t	PresEq;		///< Selecting variable to choose an equation of state
		double	Cs;		///< Speed of sound
		double	P0;		///< background pressure for equation of state
		double 	Pressure;	///< Pressure of the particle n+1

		double	Density;	///< Density of the particle n+1
		double 	Densitya;	///< Density of the particle n+1/2 (Leapfrog)
		double 	Densityb;	///< Density of the particle n-1 (Modified Verlet)
		double 	dDensity;	///< Rate of density change in time based on state equations n
		double 	RefDensity;	///< Reference Density of Particle
		double 	FPMassC;	///< Mass coefficient for fixed particles to avoid leaving particles
		double 	Mass;		///< Mass of the particle
		Vector	Displacement;	///< Density of the particle n+1

		//Mat3_t	StrainRate;	///< Global shear Strain rate tensor n
		//Mat3_t	RotationRate;	///< Global rotation tensor n
		double	ShearRate;	///< Global shear rate for fluids
		double	SBar;		///< shear component for LES

		//Mat3_t	ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1
		//Mat3_t	ShearStressa;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1/2 (Leapfrog)
		//Mat3_t	ShearStressb;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-1 (Modified Verlet)

		//Mat3_t	Sigma;		///< Cauchy stress tensor (Total Stress) n+1

		//Mat3_t	Sigmaa;		///< Cauchy stress tensor (Total Stress) n+1/2 (Leapfrog)
		//Mat3_t	Sigmab;		///< Cauchy stress tensor (Total Stress) n-1 (Modified Verlet)
		
		double Sigma_eq;	//Von Mises
		
		////////////////// PLASTIC THINGS
		//Mat3_t	Strain;		///< Total Strain n+1
		//Mat3_t	Straina;	///< Total Strain n+1/2 (Leapfrog)
		//Mat3_t	Strainb;	///< Total Strain n-1 (Modified Verlet)
		//Mat3_t  Strain_pl;	//// Plastic Strain
		
		
		double 	pl_strain,delta_pl_strain;	//Accum and incremental Effective (Von Mises) plastic strain 
		
		// BONET GRADIENT CORRECTION MATRIX
		//Mat3_t	gradCorrM;

		//Mat3_t	TIR;		///< Tensile Instability stress tensor R
		double	TI;		///< Tensile instability factor
		double	TIn;		///< Tensile instability power
		double 	TIInitDist;	///< Initial distance of particles for calculation of tensile instability

		double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Mu;		///< Dynamic viscosity coefficient of the fluid particle
		double 	MuRef;		///< Reference Dynamic viscosity coefficient
		double 	T0;		///< Yield stress for Bingham fluids
		double 	m;		///< Normalization value for Bingham fluids
		size_t	VisM;		///< Non-Newtonian viscosity method
		bool		LES;		///< Large eddy simulation using sub-particle scale
		double	CSmag;		///< Coefficient of Smagorinsky-Lilly model

		double 	G;		///< Shear modulus //TODO: Move To Material
		double 	K;		///< Bulk modulus  //TODO: Move To Material
		double	Sigmay;		///< Tensile yield stress
		size_t	Fail;		///< Failure criteria
		
		double 	Ep;		//TODO: Move To Material

		double	V;		///< Volume of a particle

		double 	h;		///< Smoothing length of the particle
		int    	LL;		///< Linked-List variable to show the next particle in the list of a cell
		int    	CC[3];		///< Current cell No for the particle (linked-list)
		int		ct;		///< Correction step for the Modified Verlet Algorithm
		double	SumKernel;	///< Summation of the kernel value for neighbour particles
		bool	FirstStep;	///< to initialize the integration scheme
		
		//LUCIANO: THERMAL PROPERTIES
		double T,k_T,cp_T,dTdt;			// Temperature, avoid permeability	
		double Ta,Tb;						//Temperature (t-1) for leapfrog
		double q_source;
		double q_conv,T_inf,h_conv;				//Different heat source terms
		double q_plheat;				//Plastic Work Heat generation
		int 	Nb;
		
		//omp_lock_t my_lock;		///< Open MP lock
    
        //SYMMETRY AND GHOST
    int     inner_mirr_part;       //NEW: SYMMETRIC PARTICLE 
    bool    is_ghost;      
    int     ghost_plane_axis;      //Assuming is cartesian
		bool    impose_vel;
    Ghost_Type ghost_type;

		// Constructor
		Particle						(int Tag, Vector const & x0, Vector const & v0, double Mass0, double Density0, double h0, bool Fixed=false);
		

	};
}; // namespace SPH

#include "Particle.cpp"

#endif //SPH_PARTICLE_H
