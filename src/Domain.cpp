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


#include "Domain.h"
#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

#include <vector>

#define MIN_PS_FOR_NBSEARCH		0.001	//TODO: MOVE TO CLASS MEMBER


#include <thrust/device_vector.h>

#include "cuda/Functions.cuh"	//For EOS at Whole Vel

#include "Mesh.h" //TO DELETE

using namespace std;

int iDivUp(int a, int b) // Round a / b to nearest higher integer value

	{ return (a % b != 0) ? (a / b + 1) : (a / b); }
	

namespace SPH {

////// GENERIC



// void General(Domain & dom)
// {
// }

// void OutPut(Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
// {
	// Prop1 = 0.0;
	// Prop2 = 0.0;
	// Prop3 = 0.0;
// }

// void InFlowCon(Vector & position, Vector & Vel, double & Den, Boundary & bdry)
// {
	// Vel = bdry.inv;
	// Den = bdry.inDensity;
// }

// void OutFlowCon(Vector & position, Vector & Vel, double & Den, Boundary & bdry)
// {
	// Vel = bdry.outv;
	// Den = bdry.outDensity;
// }

// void AllFlowCon(Vector & position, Vector & Vel, double & Den, Boundary & bdry)
// {
	// Vel = bdry.allv;
	// Den = bdry.allDensity;
// }

// Constructor
Domain::Domain ()
{
    // OutputName[0] = "Property1";
    // OutputName[1] = "Property2";
    // OutputName[2] = "Property3";
    Time    = 0.0;

    Dimension = 2;
    DomSize	= 0.0;

    Gravity	= 0.0,0.0,0.0;

    Cellfac = 2.0;

    KernelType	= 0;
    VisEq	= 0;
    Scheme	= 0;
		GradientType = 0;

    XSPH	= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    //omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;
    deltatint	= 0.0;
    deltatmin	= 0.0;
    sqrt_h_a = 0.0025;

    TRPR = 0.0;
    BLPF = 0.0;

    // InCon = & InFlowCon;
    // OutCon = & OutFlowCon;
    // AllCon = & AllFlowCon;
    // GeneralBefore = & General;
    // GeneralAfter = & General;
    // UserOutput = & OutPut;

    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
    //I = OrthoSys::I;
	
	Vol=0.;
		auto_ts = true;
		
	
	gradKernelCorr = false;
}

Domain::~Domain ()
{
	size_t Max = Particles.size();
	for (size_t i=1; i<=Max; i++)  {//Particles.DelItem(Max-i);
	//Particles.
	
	}
	Particles.clear();
}

// inline void Domain::Periodic_X_Correction(Vector & x, double const & h, Particle * P1, Particle * P2)
// {
	// if (Domsize(0)>0.0) {if (x(0)>2*Cellfac*h || x(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? x(0) -= Domsize(0) : x(0) += Domsize(0);}}
	// if (Domsize(1)>0.0) {if (x(1)>2*Cellfac*h || x(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? x(1) -= Domsize(1) : x(1) += Domsize(1);}}
	// if (Domsize(2)>0.0) {if (x(2)>2*Cellfac*h || x(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? x(2) -= Domsize(2) : x(2) += Domsize(2);}}
// }

/*inline*/ void Domain::Kernel_Set(Kernels_Type const & KT)
{
 KernelType = KT;
 if (KernelType==2) Cellfac = 3.0; else Cellfac = 2.0;
}

// inline void Domain::Viscosity_Eq_Set(Viscosity_Eq_Type const & VQ)
// {
	// VisEq = VQ;
// }

// inline void Domain::Gradient_Approach_Set(Gradient_Type const & GT)
// {
	// GradientType = GT;
// }

 inline void Domain::AdaptiveTimeStep()
 {
	 if (deltatint>deltatmin)
	 {
		 if (deltat<deltatmin)
			 deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		 else
			 deltat		= deltatmin;
	 }
	 else
	 {
		 if (deltatint!=deltat)
			 deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		 else
			 deltat		= deltatint;
	 }

//	 if (deltat<(deltatint/1.0e5))
//		 throw new Fatal("Too small time step, please choose a smaller time step initially to make the simulation more stable");
 }

// inline void Domain::AddSingleParticle(int tag, Vector const & x, double Mass, double Density, double h, bool Fixed)
// {
   	// Particles.push_back(new Particle(tag,x,Vector(0,0,0),Mass,Density,h,Fixed));
// }

void Domain::AddBoxLength(int tag, Vector const & V, double Lx, double Ly, double Lz, 
									double r, double Density, double h, int type, int rotation, bool random, bool Fixed) {
    if ( !(type == 0 || type == 1) ) {
	   	//std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		//std::cout << "0 => Hexagonal Close Packing" << std::endl;
		//std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90)) {
	   	//std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		//std::cout << "0 => " << std::endl;
		//std::cout << "0 0 0 0" << std::endl;
		//std::cout << " 0 0 0 0" << std::endl;
		//std::cout << "0 0 0 0" << std::endl;
		//std::cout << " 0 0 0 0" << std::endl;
		//std::cout << std::endl;
		//std::cout << "90 => Cubic Close Packing" << std::endl;
		//std::cout << "  0   0" << std::endl;
		//std::cout << "0 0 0 0" << std::endl;
		//std::cout << "0 0 0 0" << std::endl;
		//std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    //std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.size();

    double x,y,xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	printf("Building domain...\n");
	
    if (Dimension==3) {
    	if (type==0) {
    		//Hexagonal close packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while ( zp <= (V(2)+Lz-r) ) {
				
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r)) {
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
						y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
						z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
						if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vector(0,0,0),0.0,Density,h,Fixed));
						else    	Particles.push_back(new Particle(tag,Vector(x,y,z),Vector(0,0,0),0.0,Density,h,Fixed));
						i++;
						if ((k%2!=0) && (j%2!=0)) xp = V(0) + (2*i+(j%2)+(k%2)-1)*r; else xp = V(0) + (2*i+(j%2)+(k%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				}
				k++;
				zp = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
				//cout << "Z: "<<z<<endl;
			}
    	}
    	else {
    		//Cubic packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r)) {
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					//cout << "Y: "<<yp<<endl;
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2.0*i+1)*r;
						y = V(1) + (2.0*j+1)*r;
						z = V(2) + (2.0*k+1)*r;
						if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vector(0,0,0),0.0,Density,h,Fixed));
							else    Particles.push_back(new Particle(tag,Vector(x,y,z),Vector(0,0,0),0.0,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+1)*r; //COMMENTED BY LUCIANO
						//printf("X: %f\n",xp);
					}
					j++;
					yp = V(1) + (2.0*j+1)*r;//COMMENTED BY LUCIANO
				}
				k++;
				zp = V(2) + (2.0*k+1)*r;//COMMENTED BY LUCIANO
				//cout << "Z: "<<z<<endl;
			}
    	}

        //Calculate particles' mass in 3D
        Vector temp, Max=V;
		for (size_t i=PrePS; i<Particles.size(); i++) {
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		//cout << "BoxDimensions: "<<temp(0)<<", "<<temp(1)<<", "<<temp(2)<<", "<<endl;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.size()-PrePS);
		
		cout << "Particle mass: " << Mass <<endl;
		
		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}
    } else if (Dimension==2) {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
				j = 0;
				yp = V(1);

				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2*i+(j%2)+1)*r;
						y = V(1) + (sqrt(3.0)*j+1)*r;
						if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vector(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.push_back(new Particle(tag,Vector(x,y,0.0),Vector(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+(j%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*j+1)*r;
				}
			}
    		else
    		{
				i = 0;
				xp = V(0);

				while (xp <= (V(0)+Lx-r))
				{
					j = 0;
					yp = V(1);
					while (yp <= (V(1)+Ly-r))
					{
						x = V(0) + (sqrt(3.0)*i+1)*r;
						y = V(1) + (2*j+(i%2)+1)*r;
						if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vector(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.push_back(new Particle(tag,Vector(x,y,0.0),Vector(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						j++;
						yp = V(1) + (2*j+(i%2)+1)*r;
					}
					i++;
					xp = V(0) + (sqrt(3.0)*i+1)*r;
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
    		j = 0;
			yp = V(1);

			while (yp <= (V(1)+Ly-r))
			{
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2*i+1)*r;
					y = V(1) + (2*j+1)*r;
					if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vector(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.push_back(new Particle(tag,Vector(x,y,0.0),Vector(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
					i++;
					xp = V(0) + (2*i+1)*r;
				}
				j++;
				yp = V(1) + (2*j+1)*r;
			}

    	}
    }
		
		//cout << "Particle Count: "<<Particles.size()<< endl;

	R = r;
}

inline void Domain::Add3DCubicBoxParticles(int tag, Vector const & V, double Lx, double Ly, double Lz, 
									double r, double Density, double h) {
//	Util::Stopwatch stopwatch;
    //std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    double x,y,xp,yp;
    size_t i,j;
	size_t PrePS = Particles.size();
    
	if (Dimension==3) {
		//Cubic packing
		double z,zp;
		size_t k=0;
		zp = V(2);

		while (zp <= (V(2)+Lz-r)) {
			
			j = 0;
			yp = V(1);
			while (yp <= (V(1)+Ly-r))
			{
				//cout << "Y: "<<yp<<endl;
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2.0*i+1)*r;
					y = V(1) + (2.0*j+1)*r;
					z = V(2) + (2.0*k+1)*r;
					Particles.push_back(new Particle(tag,Vector(x,y,z),Vector(0,0,0),0.0,Density,h,false));
					i++;
					xp = V(0) + (2*i+1)*r; //COMMENTED BY LUCIANO
					//cout << "X: "<<xp<<endl;
				}
				j++;
				yp = V(1) + (2.0*j+1)*r;//COMMENTED BY LUCIANO
			}
			k++;
			cout << "Z: "<<z<<endl;
			zp = V(2) + (2.0*k+1)*r;//COMMENTED BY LUCIANO
		}
    	//Vol+=(Lx*Ly*Lz);
        Vector temp, Max=V;
		for (size_t i=PrePS; i<Particles.size(); i++) {
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		//cout << "BoxDimensions: "<<temp(0)<<", "<<temp(1)<<", "<<temp(2)<<", "<<endl;
		Vol+=temp(0)*temp(1)*temp(2);
		//double Mass = temp(0)*temp(1)*temp(2)
    }//Dimension
}

// Calculate Mass for 3D particles
inline void Domain::Calculate3DMass(double Density){
	double Mass = Vol*Density/Particles.size();
	//cout << "Particle Mass: "<<Mass<<endl;
	#pragma omp parallel for num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
	#endif
	{
		Particles[i]->Mass = Mass;
	}

}
	
//////Return half (on the quadrant) particle count from a single position in an axis
int calcHalfPartCount(const double &r, const double &R, const int xinc){
	int ypartcount = -1;
	if ( xinc > 0 ){
		ypartcount = 1;
		double yp = r;
		double xp = r + (double)(xinc - 1 ) *2.*r; 
		double rad = sqrt(yp*yp + xp*xp);
		while( rad <= R -r ){
			yp += 2.*r;
			rad = sqrt(yp*yp + xp*xp);
			ypartcount++;
		}
		ypartcount-=1;
	}
	return ypartcount;
}

void Domain::AddCylinderLength(int tag, Vector const & V, double Rxy, double Lz, 
									double r, double Density, double h, bool Fixed) {

//	Util::Stopwatch stopwatch;
    //std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy=1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	cout << "Start "<<V(0)<<", "<<V(1)<<", "<<V(2)<<endl;
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	
    if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2)+r;

		while (zp <= (V(2)+Lz-r)) {
			j = 0;
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vector(0,0,0),0.0,Density,h,Fixed));
					//	else   
					Particles.push_back(new Particle(tag,Vector(xp,yp,zp),Vector(0,0,0),0.0,Density,h,Fixed));
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
			zp += 2.*r;
		}
    
    		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.size()-PrePS);
		double Mass = Vol * Density /Particles.size();
		
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;
}

inline void Domain::AddTractionProbeLength(int tag, Vector const & V, double Rxy, double Lz_side,
											double Lz_neckmin,double Lz_necktot,double Rxy_center,
											double r, double Density, double h, bool Fixed) {

//	Util::Stopwatch stopwatch;
    //std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly, Lz;
	
	Lz = 2. * Lz_side + Lz_necktot;
	
	//Particles are tried to be aligned 
	int numpartxy;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	
	
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	double z_radiusreduction = (Lz_necktot-Lz_neckmin)/2.;
	double tan = (Rxy - Rxy_center)/(z_radiusreduction);
	double R;
	double z1 = V(2) + Lz_side - r;
	double z2 = V(2) + Lz_side + z_radiusreduction - r;
	double z3 = V(2) + Lz_side + z_radiusreduction + Lz_neckmin - r;
	double z4 = V(2) + Lz_side + Lz_necktot - r;
	
    if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2);

		while (zp <= ( V(2) + Lz - r )) {
			if 		( zp <= z1 || zp >  z4)		R = Rxy;
			else if ( zp > 	z1 && zp <= z2 )	R = Rxy - (zp - z1) * tan;
			else if ( zp >= z2 && zp < z3 )		R = Rxy_center;
			else if ( zp >= z3 && zp < z4 )		R = Rxy_center + (zp - z3) * tan;
							
			
			numpartxy = calcHalfPartCount(r, R, 1);
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, R, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//if (random) Particles.push_back(new Particle(tag,Vector((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vector(0,0,0),0.0,Density,h,Fixed));
					//	else    
					Particles.push_back(new Particle(tag,Vector(xp,yp,zp),Vector(0,0,0),0.0,Density,h,Fixed));
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
			zp = V(2) + (2.0*k+1)*r;
		}

//
// double Rxy, double Lz_side,
											// double Lz_neckmin,double Lz_necktot,double Rxy_center,
											// double r, double Density, double h, bool Fixed) {
												
		//Calculate particles' mass in 3D
		// Vector temp, Max=V;
		// for (size_t i=PrePS; i<Particles.size(); i++) {
			// if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			// if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			// if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		// }
		// Max +=r;
		// temp = Max-V;
		// double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.size()-PrePS);

		double L_cone = ( Lz_necktot - Lz_neckmin )/2.;
		double Vol = 2.0 * Lz_side * M_PI * Rxy * Rxy + 
								       Lz_neckmin * M_PI * Rxy_center * Rxy_center +
							   2.0 * L_cone * M_PI * (Rxy_center * Rxy_center + Rxy * Rxy + Rxy * Rxy_center ) / 3.0 ; //Cones
		double Mass = Vol * Density / (Particles.size()-PrePS);
		
		//cout << "Particle mass: " << Mass <<endl;
		
		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;
}


inline void Domain::Move (double dt) {
	// #pragma omp parallel for schedule (static) num_threads(Nproc)
	// for (int i=0; i<Particles.size(); i++)
		// if (Particles[i]->IsFree) {
			// if (Particles[i]->InOut>0) {
				// Particles[i]->a = 0.0;
				// if (Particles[i]->InOut == 1) {
					// Particles[i]->dDensity = 0.0;
					// Particles[i]->ZWab = 0.0;
				// } else {
					// if (BC.outDensity>0.0) {
						// Particles[i]->dDensity = 0.0;
						// Particles[i]->ZWab = 0.0;
					// }
				// }
			// }
		// Particles[i]->Move(dt,Domsize,TRPR,BLPF,Scheme,I);
		// }
}

inline void Domain::WholeVelocity() {
    //Apply a constant velocity to all particles in the initial time step
    // if (BC.allv.norm()>0.0 || BC.allDensity>0.0) {
    	// Vector vel = 0.0;
    	// double den = 0.0;


    // for (int i=0 ; i<Particles.size() ; i++) {
		// AllCon(Particles[i]->x,vel,den,BC);
    		// if (Particles[i]->IsFree && BC.allv.norm()>0.0) {
			// Particles[i]->v		= vel;
 		// }
    		// if (Particles[i]->IsFree && BC.allDensity>0.0) {
			// Particles[i]->Density	= den;
			// Particles[i]->Pressure	= EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
    		// }
    	// }
    // }
}

inline void Domain::InitialChecks() {
	//initializing identity matrix
	// if (Dimension == 2) I(2,2) = 0;

	// if (Dimension<=1 || Dimension>3) {
		// //std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		// abort();
	// }

	// if (BC.InOutFlow>0 && BC.Periodic[0])
		// throw new Fatal("Periodic BC in the X direction cannot be used with In/Out-Flow BC simultaneously");


	// //#pragma omp parallel for schedule (static) num_threads(Nproc)
	// for (int i=0; i<Particles.size(); i++) //Initializing pressure of solid and fluid particles
			// Particles[i]->Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
}

inline void Domain::TimestepCheck ()
{
	// Check the time step
	double t1,t2;
	t1 = 0.25*hmax/(CsMax);
	if (MuMax>0.0) t2 = 0.125*hmax*hmax*rhomax/MuMax; else t2 =1000000.0;

	//std::cout << "Max allowable time step using CFL = "<< std::min(t1,t2) << " S" << std::endl;
	//std::cout << "User Time Step = "<< deltatint  << " S" << std::endl;

	// if (deltatint > std::min(t1,t2))
	// throw new Fatal("Please decrease the time step to the allowable range");
}

inline void Domain::ClearNbData(){
	
	for (int i=0 ; i<Nproc ; i++) { //In the original version this was calculated after
		SMPairs[i].clear();
		FSMPairs[i].clear();
		NSMPairs[i].clear();
	}
	CellReset();
	ListGenerate();
	m_isNbDataCleared = true;
}

////In Original this function is in contact file, but here contac t
//// WeldForm 
// void Domain::AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id){
	
	// first_fem_particle_idx = Particles.size();
	// double Density =0.;
	// double h;
	// bool Fixed = false;	//Always are fixed ...
	// contact_surf_id = id;
	// trimesh = &mesh;
	
	// for ( int e = 0; e < mesh.element.size(); e++ ){
		// Vector pos = mesh.element[e]->centroid;
		// h = hfac * mesh.element[e]->radius;
		// Particles.Push(new Particle(id,pos,Vector(0,0,0),0.0,Density,h,Fixed));
		// Particles[first_fem_particle_idx + e] -> normal  = mesh.element[e] -> normal;
		// Particles[first_fem_particle_idx + e] -> element = e; 
	// }
	// cout << Particles.size() - first_fem_particle_idx << "particles added with ID " << contact_surf_id <<endl;
	// cout << first_fem_particle_idx << " is the first solid particle index."<<endl;
// }


// inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	
	// double *a;
	// calc_interact_kernel <<<1,1>>>(a);
	// // //std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	// // size_t idx_out = 1;
	// // double tout = Time;

	// // //Initializing adaptive time step variables
	// // deltat = deltatint = deltatmin	= dt;
	
	// // auto start_whole = std::chrono::steady_clock::now();

	// InitialChecks();
	// // CellInitiate(); //In NbSearch
	// // ListGenerate();
	// // PrintInput(TheFileKey);
	// // TimestepCheck();
	// // WholeVelocity();
	
	// // std::chrono::duration<double> total_time,neighbour_time;
	
	// // clock_t clock_beg;
	// // double clock_time_spent,pr_acc_time_spent,acc_time_spent;
	// // double neigbour_time_spent_per_interval=0.;
	
	// // clock_time_spent=pr_acc_time_spent=acc_time_spent=0.;


	// // //Initial model output
	// // if (TheFileKey!=NULL) {
		// // //String fn;
		// // //fn.Printf    ("%s_Initial", TheFileKey);
		// // //WriteXDMF    (fn.CStr());
		// // ////std::cout << "\nInitial Condition has been generated\n" << std::endl;
	// // }
	

	// // unsigned long steps=0;
	// // unsigned int first_step;
	
	// // int ts_nb_inc=5;	// Always > 0
	// // int ts_i=0;

	// // bool isfirst = true;
	// // bool isyielding = false;

	// // const int N = 2048 * 2048 * 4;

	// // //TimingGPU timerGPU;
  // // #define BLOCKSIZE   1024
  
  // // //SubDomain sd; //TEST
    	
	// // while (Time<=tf && idx_out<=maxidx) {
		// // StartAcceleration(sd); //TODO KERNEL
		// // //if (BC.InOutFlow>0) InFlowBCFresh();
		// // auto start_task = std::chrono::system_clock::now();
		

		// // double max = 0;
		// // int imax;
		// // #pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
		// // for (int i=0; i<Particles.size(); i++){
			// // if (Particles[i]->pl_strain>max){
				// // max= Particles[i]->pl_strain;
				// // imax=i;
			// // }
		// // }
		
		// // if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			// // ClearNbData();
			// // MainNeighbourSearch();
			// // isyielding  = true ;
		// // }
		// // if ( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE
			// // if ( ts_i == 0 ){
				// // clock_beg = clock();
				// // if (m_isNbDataCleared)
					// // MainNeighbourSearch();
				// // neigbour_time_spent_per_interval += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
			// // }
			// // isfirst = false;
		// // }
		
		// // // NEIGHBOUR SETS
		// // std::vector <int> nb(Particles.size());
		// // std::vector <int> nbcount(Particles.size());
		// // for ( size_t k = 0; k < Nproc ; k++) {
			// // for (size_t a=0; a<SMPairs[k].size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// // //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
				// // nb[SMPairs[k][a].first ]+=1;
				// // nb[SMPairs[k][a].second]+=1;
				
			// // }
		// // }	
			// // for (int p=0;p<Particles.size();p++){
			// // Particles[p]->Nb=nb[p];
		// // }

		// // auto end_task = std::chrono::system_clock::now();
		 // // neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		// // ////std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		// // GeneralBefore(*this);
		// // clock_beg = clock();
		// // PrimaryComputeAcceleration(sd);
		// // pr_acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		// // clock_beg = clock();
		// // LastComputeAcceleration(sd);
		// // acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		// // GeneralAfter(*this);
		// // steps++;
		// // //cout << "steps: "<<steps<<", time "<< Time<<", tout"<<tout<<endl;
		// // // output
		// // if (Time>=tout){
			// // // if (TheFileKey!=NULL) {
				// // // String fn;
				// // // fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				// // // WriteXDMF    (fn.CStr());

			// // // }
			// // idx_out++;
			// // tout += dtOut;
			// // total_time = std::chrono::steady_clock::now() - start_whole;
			// // //std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			// // //std::cout << "Current Time Step = " <<deltat<<std::endl;
			
			// // clock_time_spent += neigbour_time_spent_per_interval;
			// // //std::cout << "Total CPU time: "<<total_time.count() << ", Neigbour search time: " << clock_time_spent << ", Pr Accel Calc time: " <<
			// // pr_acc_time_spent << "Las Acel Calc Time" << acc_time_spent<<
			// // std::endl;
						
			// // cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
			
			// // //std::cout << "Steps count in this interval: "<<steps-first_step<<"Total Step count"<<steps<<endl;
			// // cout << "Total Neighbour search time in this interval: " << neigbour_time_spent_per_interval;
			// // cout << "Average Neighbour search time in this interval: " << neigbour_time_spent_per_interval/(float)(steps-first_step);
			// // first_step=steps;
			// // neigbour_time_spent_per_interval=0.;
		// // }
		
		// // if (auto_ts)
			// // AdaptiveTimeStep();
		// // Move(deltat);
		// // Time += deltat;
		// // //if (BC.InOutFlow>0) InFlowBCLeave(); else CheckParticleLeave ();
		
		// // if (max>MIN_PS_FOR_NBSEARCH){	//TODO: CHANGE TO FIND NEIGHBOURS
			// // if ( ts_i == (ts_nb_inc - 1) ){
				// // ClearNbData();
			// // }

			// // ts_i ++;
			// // if ( ts_i > (ts_nb_inc - 1) ) 
				// // ts_i = 0;
		
		// // }
		
	
	// // }
	

	// // //std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

// }

void Domain::AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id){
	
	// first_fem_particle_idx = Particles.Size();
	double Density =0.;
	// double h;
	bool Fixed = false;	//Always are fixed ...
	// contact_surf_id = id;
	// trimesh = &mesh;
	
	for ( int e = 0; e < mesh.element.size(); e++ ){
		Vector pos = mesh.element[e]->centroid;
		double h = hfac * mesh.element[e]->radius;
		Particles.push_back(new Particle(id,pos,Vector(0,0,0),0.0,Density,h,Fixed));
		// Particles[first_fem_particle_idx + e] -> normal  = mesh.element[e] -> normal;
		// Particles[first_fem_particle_idx + e] -> element = e; 
	}
	// cout << Particles.Size() - first_fem_particle_idx << "particles added with ID " << contact_surf_id <<endl;
	// cout << first_fem_particle_idx << " is the first solid particle index."<<endl;
}

int Domain::AssignZone(Vector &start, Vector &end, int &id){
  int partcount = 0;
  for (int a=0; a<Particles.size(); a++){
    bool included=true;
    for (int i=0;i<3;i++){
      if (Particles[a]->x(i) < start(i) || Particles[a]->x(i) > end(i))
        included = false;
    }
    if (included){
      cout << "particle "<<a<< ",ID "<<id<<endl;
      Particles[a]->ID=id;
      Particles[a]->not_write_surf_ID = true;		      
      partcount++;
    }
  }
  return partcount;
}


int ComputeCylinderParticles( double Rxy, double Lz, double r) {
  
  int part = 0;
  

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy=1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
   // if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = r;

		while (zp <= (Lz-r)) {
			j = 0;
			yp = - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp =  - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
          part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
			zp += 2.*r;
		}
	//}//Dim 3
  return part;

}

// ROWS are BC (internal)
void Domain::AddCylUniformLength(int tag, double Rxy, double Lz, 
																				double r, double Density, double h, double ang, int rows) {
	//Util::Stopwatch stopwatch;
	std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

	int PrePS = Particles.size();
	double xp,yp;
	int i,j;
	double qin = 0.03;
	srand(100);
	
	double Lx, Ly;
	
  std::pair <int,Vector> opp_sym; //Position & ID of particles

	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc;
	
	int id_part=0;
	int ghost_rows = 2;
	
	double z0;
	//if (symlength) 	z0 = r;
	//else						
    z0 = /*-Lz/2. + */ r; //CHECK: -Lz/2. - r or -Lz/2.?
	
  int radcount = Rxy / (2. * r ); //center particle is at (r,r,r)
  cout << "Radial Particle count " <<radcount<<endl;
  
	int part_per_row=0;
  std::vector <int> symm_x;
  std::vector <int> symm_y;
  int x_ghost_per_row = 0;
  //Cal
  int tgcount;
  
  std::vector <int> bc_1, bc_2; //IF apply to only boundary

  cout << "Tg Particle count " <<tgcount<<endl;
  
  int part_count = 0;
  
  if (Dimension==3) {
    	//Cubic packing
		double zp;
		int k=0;
		zp = z0;
		//Calculate row count for non ghost particles
		while (zp <= (z0 + Lz -r)){
			k++; zp += 2.0 * r;			
		}
		//cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		k = 0;zp = z0;
    cout << "Length particle count "<<last_nonghostrow+1<<endl;
    

    //First increment is in radius
		while (zp <= ( z0 + Lz - r)) {
      int rcount = 0; //Used only for ghost count
      for (double ri = 0. ; ri < Rxy; ri += 2.*r){
        //cout << "ri "<<ri<<endl;
        
        double dalpha;
        if (ri == 0.) {tgcount =1; dalpha = 0.;}
        else {

          tgcount = (ceil)((2.*M_PI* ri )/(2. * r));  
          dalpha = ang / (tgcount);         

        }
        for (int alphai = 0; alphai < tgcount; alphai++ ){
          int id = tag;
            if (ang < 2.0*M_PI){
              if (alphai<rows){
                id = 2;
                bc_1.push_back(id);
              } else if (alphai>tgcount - 1 - rows){
                id = 3;
                bc_1.push_back(id);
              }
            }
            
          xp =  /*r +*/ ri * cos (alphai*dalpha);
          yp =  /*r +*/ ri * sin (alphai*dalpha);
          cout << "XY "<<xp << ", " << yp <<endl;
          Particles.push_back(new Particle(id,Vector(xp,yp,zp),Vector(0,0,0),0.0,Density,h,false));            
        
          part_count++;
        }
        rcount++;
      } //alpha
			k++;
			zp += 2.0 * r;
		}

		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.size();
		
		cout << "Total Particle count: " << Particles.size() <<endl;
		cout << "Particle mass: " << Mass <<endl;
    
    part_count = Particles.size();

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;									
}



}; // namespace SPH
