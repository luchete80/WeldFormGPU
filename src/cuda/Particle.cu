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

#include "Particle.h"
#include "vector_math.h"
#ifndef __CUDACC__
#include <cmath>
#endif

#include "Functions.h"

static __forceinline__  __host__ __device__ void operator+=(float3 &a, const float &b) {
	a.x=a.y=a.z=b;
}


namespace SPH {

inline Particle::Particle(int Tag, float3 const & x0, float3 const & v0, double Mass0, double Density0, double h0,bool Fixed)
{
	ct = 0;

    x = x0;

    Cs		= 0.0;
    P0		= 0.0;
    PresEq	= 0;
    Alpha	= 0.0;
    Beta	= 0.0;

    v = v0;

    TI		= 0.0;
    TIn		= 4.0;
    TIInitDist  = 0.0;

    Densitya = 0.0;
    Densityb = 0.0;
    Density = Density0;
    RefDensity = Density0;

    Mass = Mass0;
    FPMassC = 1.0;
    IsFree = !Fixed;
    h = h0;
    Pressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    LL=0;
    ZWab = 0.0;
    SumDen = 0.0;
    dDensity=0.0;
    ShearRate = 0.0;
    MuRef = Mu = 0.0;
		VisM = 0;
    T0 = 0.0;
    m = 300.0;
    SumKernel = 0.0;
    G = 0.0;
    K = 0.0;
	Ep = 0.0;
	
    Material = 0;
    Fail = 0;
    
    Sigmay = 0.0;
    NoSlip = false;
    Shepard = false;
    InOut = 0;
    FirstStep = true;
    V = Mass/RefDensity;
    // RhoF = 0.0;
    IsSat = false;
    SatCheck = false;
    ShepardStep = 40;
    ShepardCounter = 0;

	LES = false;
	SBar = 0.0;
	CSmag = 0.17;
	
	///////// THERMAL ///////
	Thermal_BC = TH_BC_NONE;
	pl_strain=0;
	q_source =0;
	
	Nb=0;
	


    // set_to_zero(Strainb);
    // set_to_zero(Strain);
    // set_to_zero(Sigmab);
    // set_to_zero(Sigma);

    // set_to_zero(Sigmaa);
    // set_to_zero(ShearStress);
    // set_to_zero(ShearStressb);
    // set_to_zero(TIR);
    // set_to_zero(StrainRate);
    // set_to_zero(RotationRate);
    

}

inline void Particle::Move(double dt, float3 Domainsize, float3 domainmax, float3 domainmin, size_t Scheme, symtensor3 I)
{
//	if (Scheme == 0)
//		Move_MVerlet(I, dt);
//	else if (Scheme == 1)
		Move_Leapfrog(I, dt);
//	else if (Scheme == 2)
//		Move_Verlet(I, dt);
//	else if (Scheme == 3)
//		Move_Euler(I, dt);

	//Periodic BC particle position update
	// if (Domainsize(0)>0.0)
	// {
		// (x.x>(domainmax.x)) ? x.x -= Domainsize(0) : x.x;
		// (x.x<(domainmin(0))) ? x.x += Domainsize(0) : x.x;
	// }
	// if (Domainsize.y>0.0)
	// {
		// (x.y>(domainmax.y)) ? x.y -= Domainsize.y : x.y;
		// (x.y<(domainmin.y)) ? x.y += Domainsize.y : x.y;
	// }
	// if (Domainsize.z>0.0)
	// {
		// (x.z>(domainmax.z)) ? x.z -= Domainsize.z : x.z;
		// (x.z<(domainmin.z)) ? x.z += Domainsize.z : x.z;
	// }

}

inline void Particle::Move_Leapfrog(symtensor3 I, double dt) {
	if (FirstStep) {
		Densitya = Density - dt/2.0*dDensity;
		va = v - dt/2.0*a;
	}
	Densityb = Densitya;
	Densitya += dt*dDensity;
	Density = (Densitya+Densityb)/2.0;
	vb = va;
	va += dt*a;
	v = (va + vb)/2.0;
	x += dt*va;
	
	Displacement += dt*va;

    Mat2Leapfrog(dt);
	if (FirstStep) FirstStep = false;

}

inline void Particle::CalculateEquivalentStress () {
	// Sigma_eq	= sqrt ( Sigma(0,0)*Sigma(0,0) + Sigma(1,1)*Sigma(1,1) + Sigma(2,2)*Sigma(2,2) -
						// ( Sigma(0,0)*Sigma(1,1) + Sigma(1,1)*Sigma(2,2) + Sigma(0,0)*Sigma(2,2) ) + 
					// 3.0*(Sigma(0,1)*Sigma(0,1) + Sigma(1,2)*Sigma(1,2) + Sigma(0,2)*Sigma(0,2)));

	double J2	= 0.5*(ShearStress[0]*ShearStress[0] + 2.0*ShearStress[1]*ShearStress[3] +
						2.0*ShearStress[2]*ShearStress[6] + ShearStress[4]*ShearStress[4] +
						2.0*ShearStress[5]*ShearStress[7] + ShearStress[8]*ShearStress[8]);
	
	Sigma_eq = sqrt(3.0*J2);	
}

inline void Particle::Mat2Leapfrog(double dt) {
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	// Jaumann rate terms
	float *RotationRateT,*SRT,*RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	if (FirstStep)
		ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

	ShearStressb	= ShearStressa;
	ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;

	if (Fail == 1) {
		double J2	= 0.5*(ShearStressa(0,0)*ShearStressa(0,0) + 2.0*ShearStressa(0,1)*ShearStressa(1,0) +
						2.0*ShearStressa(0,2)*ShearStressa(2,0) + ShearStressa(1,1)*ShearStressa(1,1) +
						2.0*ShearStressa(1,2)*ShearStressa(2,1) + ShearStressa(2,2)*ShearStressa(2,2));
		//Scale back, Fraser Eqn 3-53
		ShearStressa= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;
		
		double sig_trial = sqrt(3.0*J2);
		if ( sig_trial > Sigmay) {
			double dep=( sig_trial - Sigmay)/ (3.*G + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			pl_strain += dep;
			Sigmay += dep*Ep;

		}
	}
	ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
	
	Sigma = -Pressure * OrthoSys::I + ShearStress;	//Fraser, eq 3.32

	if (FirstStep)
		Straina	= -dt/2.0*StrainRate + Strain;
	
	Strainb	= Straina;
	Straina	= dt*StrainRate + Straina;
	Strain	= 1.0/2.0*(Straina+Strainb);


	if (Fail > 1){
		//std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}



inline void Particle::translate(double dt, float3 Domainsize, float3 domainmax, float3 domainmin)
{
	x = x + dt*v + 0.5*dt*dt*a;

	// Evolve velocity
	float3 temp;
	temp = v;
	v = vb + 2*dt*a;
	vb = temp;

	//Periodic BC particle position update
	if (Domainsize.x > 0.0)
	{
		(x.x>(domainmax.x)) ? x.x -= Domainsize.x : x.x;
		(x.x<(domainmin.x)) ? x.x += Domainsize.x : x.x;
	}
	if (Domainsize.y > 0.0)
	{
		(x.y>(domainmax.y)) ? x.y -= Domainsize.y : x.y;
		(x.y<(domainmin.y)) ? x.y += Domainsize.y : x.y;
	}
	if (Domainsize.z > 0.0)
	{
		(x.z>(domainmax.z)) ? x.z -= Domainsize.z : x.z;
		(x.z<(domainmin.z)) ? x.z += Domainsize.z : x.z;
	}
}

//inline void Particle::CalcPlasticWorkHeat(){
//	
//	q_plheat 	= 	0.5*(
//					Sigma(0,0)*StrainRate(0,0) + 
//					2.0*Sigma(0,1)*StrainRate(1,0) + 2.0*Sigma(0,2)*StrainRate(2,0) + 
//					Sigma(1,1)*StrainRate(1,1) +
//					2.0*Sigma(1,2)*StrainRate(2,1) + 
//					Sigma(2,2)*StrainRate(2,2)
//					);
//}

}; // namespace SPH
