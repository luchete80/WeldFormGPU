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

//#include "cuda/Domain_d.cuh" 

#define TAU		0.005
#define VMAX	1.0

//#include "Vector.h"


void UserAcc(SPH::Domain & domi)
{
	// double vtraction;

	// if (domi.getTime() < TAU ) 
		// vtraction = VMAX/TAU * domi.getTime();
	// else
		// vtraction = VMAX;
	
	// #pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	// #ifdef __GNUC__
	// for (size_t i=0; i<domi.Particles.size(); i++)
	// #else
	// for (int i=0; i<domi.Particles.size(); i++)
	// #endif
	
	// {
		// if (domi.Particles[i]->ID == 3)
		// {
			// domi.Particles[i]->a		= Vector(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vector(0.0,0.0,vtraction);
			// domi.Particles[i]->va		= Vector(0.0,0.0,vtraction);
			// domi.Particles[i]->vb		= Vector(0.0,0.0,vtraction);
// //			domi.Particles[i]->VXSPH	= Vector(0.0,0.0,0.0);
		// }
		// if (domi.Particles[i]->ID == 2)
		// {
			// domi.Particles[i]->a		= Vector(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vector(0.0,0.0,0.0);
			// domi.Particles[i]->vb		= Vector(0.0,0.0,0.0);
			// domi.Particles[i]->VXSPH	= Vector(0.0,0.0,0.0);
		// }
	// }
}


using std::cout;
using std::endl;

//__host__		SPH::Domain dom;
	
int main(int argc, char **argv) //try
{
	cout << "Initializing"<<endl;
	SPH::Domain* dom_d;//Cannot be defined as _device
	// //OR cudamalloc((void**)&correctBool, sizeof(int));
	cudaMallocManaged(&dom_d, sizeof(SPH::Domain));
	new(dom_d) SPH::Domain();

  //SPH::Domain	dom;

  dom_d->Dimension	= 3;
  dom_d->Kernel_Set(Qubic_Spline);

//  dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

	double dx,h,rho,K,G,Cs,Fy;
	double H,L,n;

	H	= 1.;
	n	= 15.0;

	rho	= 1000.0;
	dx	= H / n;
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

     double timestep;

	// cout<<"t  = "<<timestep<<endl;
	// cout<<"Cs = "<<Cs<<endl;
	// cout<<"K  = "<<K<<endl;
	// cout<<"G  = "<<G<<endl;
	// cout<<"Fy = "<<Fy<<endl;
	
	// dom_d->GeneralAfter = & UserAcc;
	// dom_d->DomMax(0) = H;
	// dom_d->DomMin(0) = -H;
	cout << "Creating Domain"<<endl;
	dom_d->AddBoxLength(1 ,Vector ( -H/2.0 -H/20., -H/2.0 -H/20., -H/2.0 -H/20. ), H + H/20., H +H/20.,  H + H/20. , dx/2.0 ,rho, h, 1 , 0 , false, false );
		// // std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	// // double x;

    	// // for (size_t a=0; a<dom.Particles.Size(); a++)
    	// // {
    		// // x = dom.Particles[a]->x(0);
			// // dom.Particles[a]->k_T			=	3000.;
			// // dom.Particles[a]->cp_T			=	1.;
			// // dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// // dom.Particles[a]->T_inf 		= 500.;
			// // dom.Particles[a]->T				= 20.0;			
    		// // if ( x < -H/2.0 ) {
    			// // dom.Particles[a]->ID 			= 2;
    			// // dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// // // cout << "Particle " << a << "is convection BC" <<endl;
			// // }
    	// // }

        // // timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		// // cout << "Time Step: "<<timestep<<endl;
		// // //timestep=1.e-6;
		// // //0.3 rho cp h^2/k
	
		
	dom_d->WriteCSV("maz");
	
	WriteCSV_kernel<<<1,1>>>(dom_d);
// // //    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		// dom_d->Solve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.1,"test06",999);
		
        // return 0;
}
//MECHSYS_CATCH
