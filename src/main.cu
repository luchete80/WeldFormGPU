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

void report_gpu_mem()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    std::cout << "Free = " << free << " Total = " << total <<std::endl;
}


using std::cout;
using std::endl;

//__host__		SPH::Domain dom;
	
int main(int argc, char **argv) //try
{
	cout << "Initializing"<<endl;
	SPH::Domain dom;//Cannot be defined as _device
	// //OR cudamalloc((void**)&correctBool, sizeof(int));
	// cudaMallocManaged(&dom, sizeof(SPH::Domain));
	// new(dom) SPH::Domain();
	
	SPH::Domain_d *dom_d;
	report_gpu_mem();
	cudaMallocManaged(&dom_d, sizeof(SPH::Domain));
	report_gpu_mem();
	dom_d->SetDimension(dom.Particles.size());

  dom.Dimension	= 3;
  dom.Nproc	= 4;
  //dom.Kernel_Set(Qubic_Spline);

//  dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

	double dx,h,rho;
	double H,L,n;

	H	= 1.;
	n	= 15.0;

	rho	= 1000.0;
	dx	= H / n;
	h	= dx*1.2; //Very important
	//Cs	= sqrt(K/rho);

  double timestep;

	// cout<<"t  = "<<timestep<<endl;
	// cout<<"Cs = "<<Cs<<endl;
	// cout<<"K  = "<<K<<endl;
	// cout<<"G  = "<<G<<endl;
	// cout<<"Fy = "<<Fy<<endl;
	
	// dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = H;
	dom.DomMin(0) = -H;
  dom.GeneralAfter = & UserAcc;
	cout << "Creating Domain"<<endl;
	dom.AddBoxLength(1 ,Vector ( -H/2.0 -H/20., -H/2.0 -H/20., -H/2.0 -H/20. ), H + H/20., H +H/20.,  H + H/20. , dx/2.0 ,rho, h, 1 , 0 , false, false );
	cout << "Particle count:" <<dom.Particles.size()<<endl;
	
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x =  new double3 [dom.Particles.size()];
	for (int i=0;i<dom.Particles.size();i++){
		//cout <<"i; "<<i<<endl;
		//x[i] = make_double3(dom.Particles[i]->x);
		x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
	}
	// int size = dom.Particles.size() * sizeof(double3);
	cout << "Copying to device..."<<endl;
	//cudaMemcpy(dom_d->x, x, size, cudaMemcpyHostToDevice);
	cout << "copied"<<endl;
	
	// //Temporary, NB Search in GPU
	cout << "Cell Initiate..."<<endl; dom.CellInitiate();
	cout << "Generating List..."<<endl;	dom.ListGenerate();
	
	cout << "Nb Searching..."<<endl;	dom.MainNeighbourSearch(); 
	cout << "Done"<<endl;
	
	//Creating 2 arrays of nb (TODO: Which is faster 2D or flattened array?)
	//First, counting size of all nbs
	int MAXNB_PPART = 70;
	int** nb2d = new int*[dom.Particles.size()];
	for(int i = 0; i < dom.Particles.size(); ++i)
		nb2d[i] = new int[MAXNB_PPART];
	cout << "Creating 2d array..."<<endl;
	int nbcount = 0;
	int *nb = new int[dom.Particles.size()]; //nb count of each particle
	for (int i = 0; i < dom.Particles.size(); ++i) nb[i] = 0;
	//THIS WILL BE DONE IN SEARCH
	for ( int k = 0; k < dom.Nproc ; k++) {
		cout<< "pair size"<<dom.SMPairs[k].size()<<endl;
		for (int a=0; a<dom.SMPairs[k].size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a,k: "<<a<<" "<<k<<endl;
			//nb2d[dom.SMPairs[k][a].first] [nb[dom.SMPairs[k][a].first]]  = dom.SMPairs[k][a].second;
			//nb2d[dom.SMPairs[k][a].second][nb[dom.SMPairs[k][a].second]] = dom.SMPairs[k][a].first;
			nb[dom.SMPairs[k][a].first]++;
			nb[dom.SMPairs[k][a].second]++;
			nbcount+=2; // Total nb count (for flattened array)
		}
	}
	cout << "NbCount[0] "<< nb[0]<<endl;
	cout << "Done."<<endl;
	
	cout << "Allocating flattened array..."<<endl;
	//FLATENED ARRAY
	int *nb_part =  new int [nbcount]; //This could be sized only once with max nb count
	int *nb_offs =  new int [dom.Particles.size()+1];
	for (int n=0; n<dom.Particles.size();n++) 	nb_offs[n] = 0;
	for (int n=0; n<nbcount;n++) 								nb_part[n] = 0;
  cout << "Nb array size"<< nbcount<<endl;
	
	int i=0;
	cout << "Creating flattened array..."<<endl;	
	// for (int n=0; n<dom.Particles.size();n++) {
		// for (int k=0; k< nb[n] ;k++) {
			// nb_part[i] = nb2d[n][k];
			// i++;		
		// }
		// nb_offs[n+1]=i;
	// }
  cout<< "done"<<endl;
	
	cout << "Allocating in device.."<<endl;
	//Device side
	//cudaMalloc((void **) dom_d->neib_part, 	(nbcount) * sizeof (int));
	//cudaMemcpy(dom_d->neib_part, nb_part, nbcount * sizeof(int), cudaMemcpyHostToDevice);
	//nb offset or count already initiated
	//cudaMemcpy(dom_d->neib_offs, nb_offs, (dom.Particles.size()+1) * sizeof(int), cudaMemcpyHostToDevice);
	
	cout << "done"<<endl;
	cout << "Setting values"
	dom_d->SetConductivity(3000.);
	dom_d->SetHeatCap(1.);
	cout << "done."<<endl;
	
		// // std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	// // double x;

    	for (size_t a=0; a<dom.Particles.size(); a++)
    	{
    		// x = dom.Particles[a]->x(0);
			// dom.Particles[a]->k_T			=	3000.;
			// dom.Particles[a]->cp_T			=	1.;
			// dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// dom.Particles[a]->T_inf 		= 500.;
			// dom.Particles[a]->T				= 20.0;		
			dom.Particles[a]->IsFree	= true;		
    		// if ( x < -H/2.0 ) {
    			// dom.Particles[a]->ID 			= 2;
    			// dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// // cout << "Particle " << a << "is convection BC" <<endl;
			// }
    	}

        // // timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		// // cout << "Time Step: "<<timestep<<endl;
		// // //timestep=1.e-6;
		// // //0.3 rho cp h^2/k
	
		
	// dom.WriteCSV("maz");
	
	// WriteCSV_kernel<<<1,1>>>(&dom);


	//dom_d->ThermalSolve(/*tf*/1.01);
		
        // return 0;
				
	//cudaFree(dom_d);
	report_gpu_mem();
}
//MECHSYS_CATCH
