#include "Domain_d.cuh"
#include "Functions.cuh"
#include <iostream>

#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock
#include "tensor3.cu" //INLINE
#include "InteractionAlt.cu"

#include "Geometry.cu"

#include "cuNSearch.h"
//For Writing file

//This is temporary since can be used a delta_pl_strain for each particle
#define MIN_PS_FOR_NBSEARCH		1.e-6//TODO: MOVE TO CLASS MEMBER
#include "Mesh.cuh"
#include "Mesh.cuh"
#include "Energy.cu"

#define TAU		0.005
#define VMAX	10.0

using namespace std;
namespace SPH{
  
void Domain_d::MechFraserSolve(const double &tf, const double &dt_out){

	int N = particle_count;
	threadsPerBlock = 256; //Or BlockSize
	blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
  Time =0.;
	
	isfirst_step =true;

	step = 0;						//Par of domain_h
	clock_t clock_beg;
	double time_spent;
	clock_beg = clock();
	
	//TimestepCheck(0.7,h,Cs);

	double t_out;
	t_out = dt_out;
	
	double stress_time,forces_time,accel_time,pressure_time,move_time;
	
	clock_t clock_beg_int;

	stress_time = forces_time = accel_time = pressure_time = move_time = 0.;

	cudaMemcpy(x_h, x, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);		


	//Make the nb search at first
	
	vector < Real3> pos;
  // positions.reserve(dom.Particles.size());
	for (unsigned int i = 0; i < particle_count; i++) {
    std::array<Real, 3> x ={{ x_h[i].x,
                              x_h[i].y,
                              x_h[i].z
                            }};
		pos.push_back(x);
	}
	
	cout << "Initializing nb search data.."<<endl;
	double radius = 2.0*h_glob;
	//This bypass the original constructor 
  //TODO: make this
	
	nb_search.deviceData = std::make_unique<cuNSearchDeviceData>(radius);
	nb_search.set_radius(radius);
	
	cuNSearch::NeighborhoodSearch nsearch(radius);
	cout << "Done."<<endl;

	auto pointSetIndex = nsearch.add_point_set(pos.front().data(), pos.size(), true, true);
	auto &pointSet = nsearch.point_set(0);
	
	int *nb_part_h =  new int [particle_count * 100]; //This could be sized only once with max nb count
	int *nb_offs_h =  new int [particle_count + 1];
	
	auto points = pointSet.GetPoints();
	
	int ts_i=0;
	int ts_nb_inc = 5;
	
	bool is_yielding = false;
	double max_pl_strain = 0.;
  cout << "First Rigid Contact Particle: "<<first_fem_particle_idx<<endl;
  
	//First time find nbs
	for (int i=0; i <particle_count;i++){
	((Real3*)points)[i][0] = x_h[i].x;
	((Real3*)points)[i][1] = x_h[i].y;
	((Real3*)points)[i][2] = x_h[i].z;
	}		
	// TODO: FIX THIS! 
	//zsort is much faster than traditional, but particle order and nb changes
	//nsearch.z_sort();
	//nsearch.point_set(pointSetIndex).sort_field((Real3*)nsearch.point_set(pointSetIndex).GetPoints());
	nsearch.find_neighbors();	
		// testNeighboursKernel<<< blocksPerGrid,threadsPerBlock >>>(	0,
		// CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
		// CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
		// CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)
		// );
    
  int count = 1; //step
  
  //totmass = 1.;
  
  this->id_free_surf = 1;
  
  cout << "particle_count "<<particle_count<<endl;
  SetVelKernel<<<blocksPerGrid,threadsPerBlock >>>(this,make_double3(0.,0.,0.));
  cudaDeviceSynchronize(); 
  
  cout << "Main Loop "<<endl;

  int pcount = particle_count;
  if (contact){
    pcount = solid_part_count;
    //THIS IS DONE OUTSIDE MAIN LOOP IN WELDFORM CPU VERSION
    if (this->trimesh != NULL){
      printf("Calculating plane coefficients...\n");
      CalcSpheresKernel<<<blocksPerGrid,threadsPerBlock >>>(this->trimesh);
      cudaDeviceSynchronize();
    }
    else printf("MESH NOT DEFINED\n");
  }
  
  //ONLY FOR TESTING 
  test_h = new double [particle_count];
  while (Time<tf) {
	
		if ( ts_i == 0 && is_yielding ){
			//cout << "Searching nbs"<<endl; 
			/////////////////////////////////////////
			// UPDATE POINTS POSITIONS
			//TODO: THIS HAS TO BE DONE WITH KERNEL
			for (int i=0; i <particle_count;i++){
			((Real3*)points)[i][0] = x_h[i].x;
			((Real3*)points)[i][1] = x_h[i].y;
			((Real3*)points)[i][2] = x_h[i].z;
			}		
			// TODO: FIX THIS! 
			//zsort is much faster than traditional, but particle order and nb changes
			//nsearch.z_sort();
			//nsearch.point_set(pointSetIndex).sort_field((Real3*)nsearch.point_set(pointSetIndex).GetPoints());
			nsearch.find_neighbors();

		}//ts_i == 0
	
		//cout << "
		
		//cout<<"--------------------------- BEGIN STEP "<<step<<" --------------------------"<<endl; 
		//This was in Original LastCompAcceleration
		clock_beg_int = clock();

    ////////TODO: MOVE TO A SINGLE HOST DOMAIN FUNCTION
    CalcDensIncKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)		
		);
    cudaDeviceSynchronize(); //REQUIRED!!!!
       
    ///////Like CPU version, this is calculated with density. 
    UpdateDensityKernel<<<blocksPerGrid,threadsPerBlock >>>(this,deltat);		
    cudaDeviceSynchronize(); //REQUIRED!!!!

    CalcRateTensorsKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)		
		);
    cudaDeviceSynchronize(); //REQUIRED!!!!

		//If kernel is the external, calculate pressure
		//Calculate pressure!
		PressureKernelExt<<<blocksPerGrid,threadsPerBlock >>>(p,PresEq,Cs,P0,rho,rho_0,particle_count);
		cudaDeviceSynchronize();

		clock_beg_int = clock();
		StressStrainKickDriftKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
		cudaDeviceSynchronize();
		stress_time += (double)(clock() - clock_beg_int) / CLOCKS_PER_SEC;
    
		CalcAccelKernel	<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)		
		);
    cudaDeviceSynchronize(); //REQUIRED!!!!
    
    if (contact){
      CalculateSurfaceKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors),		    
      /*id,*/
      totmass);
      cudaDeviceSynchronize(); //REQUIRED!!!!
      //
      CalcContactNbKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)    
      );
      cudaDeviceSynchronize(); //REQUIRED!!!!    

      
      CalcContactForcesKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors) 
      );
      cudaDeviceSynchronize();
    }
    

					
		//TODO: CHANGE this to an interleaved reduction or something like that (see #84)
		if (!is_yielding){
			cudaMemcpy(pl_strain_h, pl_strain, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);
			for (int i=0;i<particle_count;i++){
				if ( pl_strain_h[i] > max_pl_strain )
					max_pl_strain = pl_strain_h[i];
			}
			
			if ( max_pl_strain > MIN_PS_FOR_NBSEARCH ){
				is_yielding = true;
				cout << "Now is yielding"<<endl;
			}
		}
    
    GeneralAfter(*this); //SET BCS
    
		UpdatePosFraserKernel<<<blocksPerGrid,threadsPerBlock >>>(this,deltat);
    cudaDeviceSynchronize(); //REQUIRED!!!!

    UpdateVelKernel<<<blocksPerGrid,threadsPerBlock >>>(this,deltat);
    cudaDeviceSynchronize();
    
    GeneralAfter(*this); //REINFORCE BCs AGAIN
    
    //TODO: Remaining plastic work heat here
    
    ThermalCalcs(deltat);
    
    CalcIntEnergyKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
		cudaDeviceSynchronize();
    
    if (contact){
      if (this->trimesh != NULL){
      MeshUpdateKernel<<<blocksPerGrid,threadsPerBlock >>>(this->trimesh, deltat);
      cudaDeviceSynchronize();
      } else {
        cout << "No contact mesh defined."<<endl;
      }
    }
    
    if (contact) {
      UpdateContactParticlesKernel<<< blocksPerGrid,threadsPerBlock >>>(this);
      cudaDeviceSynchronize();
    }
      
		if (isfirst_step) isfirst_step = false;
		Time +=deltat;		
	
		if (Time >= t_out) {		
			cudaMemcpy(ID_h, ID, sizeof(int) * particle_count, cudaMemcpyDeviceToHost);	
			cudaMemcpy(x_h, x, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
			cudaMemcpy(u_h, u, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
			cudaMemcpy(v_h, v, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
			cudaMemcpy(a_h, a, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
      cudaMemcpy(nb_h, CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts), 
                          sizeof(unsigned int) * particle_count, cudaMemcpyDeviceToHost);	
			
			cudaMemcpy(p_h, p, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);	
			
			cudaMemcpy(rho_h, rho, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);
			cudaMemcpy(sigma_eq_h, sigma_eq, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);	
			cudaMemcpy(pl_strain_h, pl_strain, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);
      
      cudaMemcpy(contneib_count_h,contneib_count, sizeof(int) * particle_count, cudaMemcpyDeviceToHost);
      if (contact)
        cudaMemcpy(contforce_h, contforce, sizeof(double3) * pcount, cudaMemcpyDeviceToHost);	
			cudaMemcpy(normal_h, normal, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
      
      if (thermal_solver)
        cudaMemcpy(T_h, T, sizeof(double) * particle_count, cudaMemcpyDeviceToHost);	      

      //TODO: MAKE A LIST WITH CONTACT PARTICLES
      double cont_force_sum = 0.;
      for (int i=0;i<particle_count;i++)
        if (ID_h[i]==id_free_surf)
          cont_force_sum+=length(contforce_h[i]);
			
			char str[10];
			sprintf(str, "out_%d.csv", count);
      count++;
			WriteCSV(str);
			
			t_out += dt_out;
			time_spent = (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
			cout << "-------------------------\nTime "<<Time<<", GPU time "<<time_spent<<endl;
			cout << "Current time step: "<< deltat << endl;
			cout << "Forces calc: "			<<forces_time<<endl;
			cout << "Stresses calc: "		<<stress_time<<endl;
			
			double3 max= make_double3(0.,0.,0.);
      double max_ps = 0.;
			for (int i=0;i<particle_count;i++){
				//cout << "Particle " << i << "Vel: "<< v_h[i].x<<", "<<v_h[i].y<< ", "<< v_h[i].z<<endl;
				//cout << "Particle " << i << "Acc: "<< a_h[i].x<<", "<<a_h[i].y<< ", "<< a_h[i].z<<endl;
				if (u_h[i].x>max.x) max.x = u_h[i].x;
				if (u_h[i].y>max.y) max.y = u_h[i].y;
				if (u_h[i].z>max.z) max.z = u_h[i].z;
        if ( pl_strain_h[i] > max_ps) max_ps = pl_strain_h[i];
			}
			cout << "Max disp: "<< max.x<<", "<<max.y<<", "<<max.z<<endl;
      cout << "Max pl_strain: "<<max_ps<<endl;
      if (contact)
        cout << "Contact Force Sum: "<< cont_force_sum << endl;
      //cout << "Int Energy "<< int_energy_sum<<endl;
		}
		time_spent = (double)(clock() - clock_beg) / CLOCKS_PER_SEC;	
		step ++;
		//cout<<"--------------------------- END STEP, Time"<<Time <<", --------------------------"<<endl; 

		ts_i ++;
		if ( ts_i > (ts_nb_inc - 1) ) 
			ts_i = 0;
		
	}//while <tf


	time_spent = (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
	
	printf("Total steps: %d, time spent %f\n",step, time_spent);
	
	delete nb_part_h;
	delete nb_offs_h;

}

};//SPH