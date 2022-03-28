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

#include "cuda/Domain_d.cuh" 
#include "cuda/Mechanical.cu" 

#define TAU		0.005
#define VMAX	1.0

#include <sstream>
#include <fstream> 
#include <iostream>

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

void report_gpu_mem()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    std::cout << "Free = " << free << " Total = " << total <<std::endl;
}


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
//https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

using std::cout;
using std::endl;

//__host__		SPH::Domain dom;

void WriteCSV(char const * FileKey, double3 *x, double3 *varv, int count){
	std::ostringstream oss;
	std::string fn(FileKey);
	
	oss << "X, Y, Z, dx, dy, dz"<<endl;;
	
	//#pragma omp parallel for schedule(static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	// #else
	for (int i=0; i<count; i++)//Like in Domain::Move
	//#endif
	{
			oss << x[i].x<<", "<<x[i].y<<", "<<x[i].z<<", "<<varv[i].x<<
", "<<varv[i].y<<
", "<<varv[i].z<<			endl;
		
		//Particles[i]->CalculateEquivalentStress();		//If XML output is active this is calculated twice
		//oss << Particles[i]->Sigma_eq<< ", "<< Particles[i]->pl_strain <<endl;
	}

	std::ofstream of(fn.c_str(), std::ios::out);
	of << oss.str();
	of.close();
}

int main(int argc, char **argv) //try
{
	
	
	cout << "Initializing"<<endl;
	SPH::Domain dom;//Cannot be defined as _device
	// //OR cudamalloc((void**)&correctBool, sizeof(int));
	// cudaMallocManaged(&dom, sizeof(SPH::Domain));
	// new(dom) SPH::Domain();
	
	SPH::Domain_d *dom_d;
	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(SPH::Domain)) );
	report_gpu_mem();

  dom.Dimension	= 3;
  dom.Nproc	= 4;
  //dom.Kernel_Set(Qubic_Spline);

//  dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

	double dx,h,rho,K,G;
	double R,L,n;

	R	= 0.15;
	L	= 0.56;
	n	= 30.0;		//in length, radius is same distance

	rho	= 2700.0;
	K	= 6.7549e10;
	G	= 2.5902e10;
	
	dx = 0.030; //THIS IS FOR TESTING Original 6,5mm, 8mm 10mm, 12,5 and 15mm
	h	= dx*1.2; //Very important

	double Cs	= sqrt(K/rho);

  double timestep = (0.2*h/(Cs));

	// cout<<"t  = "<<timestep<<endl;
	// cout<<"Cs = "<<Cs<<endl;
	// cout<<"K  = "<<K<<endl;
	// cout<<"G  = "<<G<<endl;
	// cout<<"Fy = "<<Fy<<endl;
	
	// dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = L;
	dom.DomMin(0) = -L;
  dom.GeneralAfter = & UserAcc;
	cout << "Creating Domain"<<endl;
	dom.AddCylinderLength(1, Vector(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 
	cout << "Particle count:" <<dom.Particles.size()<<endl;
	
	dom_d->SetDimension(dom.Particles.size());	 //AFTER CREATING DOMAIN
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x =  new double3 [dom.Particles.size()];
	for (int i=0;i<dom.Particles.size();i++){
		//cout <<"i; "<<i<<endl;
		//x[i] = make_double3(dom.Particles[i]->x);
		x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
	}
	int size = dom.Particles.size() * sizeof(double3);
	cout << "Copying to device..."<<endl;
	cudaMemcpy(dom_d->x, x, size, cudaMemcpyHostToDevice);
	cout << "copied"<<endl;
	
	cout << "Setting values"<<endl;
	dom_d->SetDensity(rho);
	dom_d->Set_h(h);
	cout << "done."<<endl;

	double *m =  new double [dom.Particles.size()];
	for (size_t a=0; a<dom.Particles.size(); a++)
		m[a] = dom.Particles[a]->Mass;
	cudaMemcpy(dom_d->m, m, dom.Particles.size() * sizeof(double), cudaMemcpyHostToDevice);	
		
		// // std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	// // double x;

	//MODIFY
	double *T 			=  new double [dom.Particles.size()];
	int 	*BC_type 	=  new int 		[dom.Particles.size()];
	int bcpart = 0;
	for (size_t a=0; a<dom.Particles.size(); a++){
		double xx = dom.Particles[a]->x(0);
		BC_type[a]=0;
		T[a] = 20.;
		if ( xx < -L/2.0 ) {
			bcpart++;
			BC_type[a]=1;
		}
	}		
	cout << "BC particles"<<bcpart<<endl;
	cudaMemcpy(dom_d->T, T, dom.Particles.size() * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dom_d->BC_T, BC_type, dom.Particles.size() * sizeof(int), cudaMemcpyHostToDevice);
	
	dom_d->Alpha = 0.0;//For all particles		
	dom_d->SetShearModulus(G);	// 
	for (size_t a=0; a<dom.Particles.size(); a++)
	{
		//dom.Particles[a]->G				= G; 
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs			= Cs;

		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
		if ( z < 0 ){
			dom.Particles[a]->ID=2;	
		}
		if ( z > L )
			dom.Particles[a]->ID=3;
	}
	
	dom_d->SetFreePart(dom); //All set to IsFree = true in this example
	dom_d->SetID(dom); 
	dom_d->SetCs(dom);
	
	dom_d->SetSigmay(300.e6);
	

        // // timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		// // cout << "Time Step: "<<timestep<<endl;
		// // //timestep=1.e-6;
		// // //0.3 rho cp h^2/k
	
		
	// dom.WriteCSV("maz");
	
	// WriteCSV_kernel<<<1,1>>>(&dom);

	cout << "Solving "<<endl;
	//CheckData<<<1,1>>>(dom_d);
	//cudaDeviceSynchronize(); //Crashes if not Sync!!!
	
	dom_d->deltat = timestep;
	//dom_d->deltat = 2.0e-8;
	
	cout << "Time Step: "<<dom_d->deltat<<endl;
	WriteCSV("test_inicial.csv", x, dom_d->u_h, dom.Particles.size());
	//dom_d->MechSolve(0.00101 /*tf*//*1.01*/,1.e-4);
	//dom_d->MechSolve(100*timestep + 1.e-10 /*tf*//*1.01*/,timestep);
	dom_d->auto_ts = true;
	dom_d->MechSolve(0.002001,1.0e-4);
	//dom_d->MechSolve(2.*timestep,timestep);

	cudaMemcpy(T, dom_d->T, sizeof(double) * dom.Particles.size(), cudaMemcpyDeviceToHost);	
	
        // return 0;
	//WriteCSV("test.csv", x, dom_d->u_h, dom.Particles.size());
	dom_d->WriteCSV("test.csv");
	
	cudaFree(dom_d);
	//report_gpu_mem();
	cout << "Program ended."<<endl;
}
//MECHSYS_CATCH
