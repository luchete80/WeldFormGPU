#include "Domain_d.cuh"
#include "Functions.cuh"
#include <iostream>

#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock
#include "tensor3.cu" //INLINE
#include "Interaction.cu"

#include "InteractionAlt.cuh"

#include "Geometry.cu"
#include "Contact.cu"

#include "cuNSearch.h"
#include "Material.cuh"
//For Writing file

//This is temporary since can be used a delta_pl_strain for each particle
#define MIN_PS_FOR_NBSEARCH		1.e-6//TODO: MOVE TO CLASS MEMBER
#include "Mesh.cuh"

using namespace std;

namespace SPH{

__global__ void WholeVelocityKernel(Domain_d *dom_d){
	dom_d->WholeVelocity();
}
//Called by __global__
void __device__ Domain_d::WholeVelocity() {
		int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
	}
}

//////////////////////////////////////////////////
///// THIS IS ONLY WHEN THERE ARE FIXED PARTICLES!
//////////////////////////////////////////////////
__device__ void Domain_d::PrimaryComputeAcceleration (/*int i*/) {
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {
		double3 xij;
		double h_,K;
	int Dimension = 3;
	
	// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
	//Same Material pairs
	int neibcount;
	#ifdef FIXED_NBSIZE
	neibcount = neib_offs[i];
	#else
	neibcount =	neib_offs[i+1] - neib_offs[i];
	#endif
	// printf("neibcount %d\n",neibcount);
	// printf("Nb indexed,i:%d\n",i);
	for (int k=0;k < neibcount;k++) { //Or size
		// P1	= FSMPairs[k][a].first;
		// P2	= FSMPairs[k][a].second;
		int j = NEIB(i,k); //TODO; MAKE A FIXED PAIR
		xij	= x[i]-x[j];
		h_	= (h[i]+ h[j])/2.0;
		double nxij = length(xij);
		//Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);
		//(size_t const & Dim, size_t const & KT, double const & q, double const & h);
		K	= Kernel(Dimension, 0, nxij/h_, h_);
		if ( !IsFree[i] ) {
				SumKernel[i] += K;
				p[i]	+= p[j] * K /*+ dot(Gravity,xij)*rho[j]*K*/;
				sigma[i] 	 = sigma[i] + K * sigma[j];
				if (NoSlip[i])		NSv[i] 	+= v[j] * K;
		} else {
				SumKernel[j] += K;
				p[j]	+= p[i] * K /*+ dot(Gravity,xij)*rho[i]*K*/;
				sigma[j]	 = sigma[j] + K * sigma[i];
				if (NoSlip[j])		NSv[j] 	+= v[i] * K;

		}	
	}//FIXED neibcount k

	////////////////////////////////////////////////////////////////
	// // Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
	////////////////////////////////////////////////////////////////
	
	// #pragma omp parallel for schedule (static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<FixedParticles.Size(); i++)
	// #else
	// for (int i=0; i<FixedParticles.Size(); i++)
	// #endif
		// if (Particles[FixedParticles[i]]-> ID != contact_surf_id)  //ADDED TO Prevent adding surface (rigid contact) particles
		// if (Particles[FixedParticles[i]]->SumKernel!= 0.0) {
			// size_t a = FixedParticles[i];
			// Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
			// Particles[a]->Sigma	= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
			// if (Particles[a]->NoSlip)	Particles[a]->NSv	= Particles[a]->NSv/Particles[a]->SumKernel;

			// // Tensile Instability for fixed soil and solid particles
			// if (Particles[a]->TI > 0.0)
			// {
				// // XY plane must be used, It is very slow in 3D
				// if (Dimension == 2) {
					// double teta, Sigmaxx, Sigmayy, C, S;
					// if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0)
						// teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1)));
					// else
						// teta = M_PI/4.0;

					// C = cos(teta);
					// S = sin(teta);
					// Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
					// Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
					// if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
					// if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
					// Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					// Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					// Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				// }
				// else {
					// Mat3_t Vec,Val,VecT,temp;
					// Rotation(Particles[a]->Sigma,Vec,VecT,Val);
					// double pc_ti_inv_d2=Particles[a]->TI/(Particles[a]->Density*Particles[a]->Density);//Precompute some values
					// // if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
					// // if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
					// // if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
					// if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					// Mult(Vec,Val,temp);
					// Mult(temp,VecT,Particles[a]->TIR);
				// }
			// }
		// }

	}//i<part_count
}

__device__ inline void Domain_d::UpdateVel(double dt){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {
    	v[i] += dt*a[i];   
  }
}

__global__ inline void UpdateVelKernel(Domain_d *dom, double dt) {
  dom->UpdateVel(dt);
}

__device__ inline void Domain_d::SetVel(double3 vi){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	if ( i < particle_count ) {
    	v[i] = vi;   
  }
}

__global__ void SetVelKernel(Domain_d *dom,  double3 v){
  dom->SetVel(v);
}

__device__ inline void Domain_d::UpdatePos(double dt){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {
    	x[i] += dt*v[i];   
      u[i] += dt*v[i];
  }
}

__global__ void UpdatePosKernel(Domain_d *dom, double dt) {
  dom->UpdatePos(dt);
}

//THIS UPDATES POSITION FROM VELOCITY AND ACCEL, NOT ONLY VELOCITY
__device__ inline void Domain_d::UpdatePosFraser(double dt){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {
    double3 u_ = dt * v[i] + 0.5 * a[i] * dt * dt;
    	x[i] += u_;   
      u[i] += u_;
  }
}

__global__ void UpdatePosFraserKernel(Domain_d *dom, double dt) {
  dom->UpdatePosFraser(dt);
}


__global__ void UpdateDensityKernel(Domain_d *dom, double dt) {
  dom->UpdateDensity(dt);
}

__device__ inline void Domain_d::UpdateDensity(double dt){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	
	if ( i < solid_part_count ) {
    	rho[i] += dt*drho[i];   
  } 
}

__host__ inline void Domain_d::CalcDensity(const double &dt, int blocksPerGrid,int threadsPerBlock){
    // CalcDensIncKernel<<<blocksPerGrid,threadsPerBlock >>>(this,
      // CudaHelper::GetPointer(nsearch.deviceData->d_NeighborCounts),
      // CudaHelper::GetPointer(nsearch.deviceData->d_NeighborWriteOffsets),
      // CudaHelper::GetPointer(nsearch.deviceData->d_Neighbors)		
		// );
    // cudaDeviceSynchronize(); //REQUIRED!!!!
    
}

//IS NOT NECESSARY TO PASS ENTIRE DOMAIN!
void __global__ MoveKernelExt(double3 *v, double3 *va, double3 *vb,
													double *rho, double *rhoa, double *rhob,double *drho,
													double3 *x, double3 *a,
													double3 *u, /*Mat3_t I, */double dt,
													bool FirstStep, int particle_count)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {
	if (FirstStep) {
		//printf("First Step\n");
		rhoa[i] = rho[i] - dt/2.0*drho[i];
		va[i] = v[i] - dt/2.0*a[i];
		//printf("First Step!\n");
	}
	rhob[i] = rhoa[i];
	rhoa[i] += dt*drho[i];
	rho[i] = (rhoa[i]+rhob[i])/2.0;
	// if (i==1250){
		// printf("Move. particle 1250 rho %f rhoa %f rhob %f\n",rho[i],rhoa[i],rhob[i]);
	// }
	vb[i] = 	va[i];
	va[i] += dt*a[i];
	v[i] = (va[i] + vb[i])/2.0;
	x[i] += dt*va[i];
	
	u[i] += dt*va[i];

    // Mat2Leapfrog(dt);
	// if (FirstStep) FirstStep = false;
	}
}

// __device__ PartData_d::PrimaryComputeAcceleration(){
	// // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
	
	
	// // Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
// }

// void __global__ MechSolveKernel (double dt, PartData_d *partdata) {
	// int i = threadIdx.x+blockDim.x*blockIdx.x;
	// dTdt[i] = 0.;
	
	// int neibcount;
	// #ifdef FIXED_NBSIZE
	// neibcount = neib_offs[i];
	// #else
	// neibcount =	neib_offs[i+1] - neib_offs[i];
	// #endif
	// printf("Solving\n");
	// for (int k=0;k < neibcount;k++) { //Or size
		// //if fixed size i = part * NB + k
		// //int j = neib[i][k];
		// int j = NEIB(i,k);
		// printf("i,j\n",i,j);
		// double3 xij; 
		// xij = x[i] - x[j];
		// double h_ = (h[i] + h[j])/2.0;
		// double nxij = length(xij);
		
		// double GK	= GradKernel(3, 0, nxij/h_, h_);
		// //		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
		// //   mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
		// dTdt[i] += m[j]/rho[j]*( 4.0*k_T[i]*k_T[j]/(k_T[i]+k_T[j]) * (T[i] - T[j])) * dot( xij , GK*xij )/(nxij*nxij);
	// }
	// dTdt[i] *=1/(rho[i]*cp[i]);
// }

__device__ __forceinline__ void Domain_d::LastComputeAcceleration(){
	//partdata->();
}

__global__ void PressureKernelExt(double *p, double *PresEq, double *Cs, double *P0,double *Density, double *RefDensity, int particle_count){
	
	int i = threadIdx.x + blockDim.x*blockIdx.x;	
	// if (i == 1250){
		// printf("PresEq[i], Cs[i], P0[i],Density[i], RefDensity[i]: %f %f %f %f %f \n",PresEq[i], Cs[i], P0[i],Density[i], RefDensity[i]);
	// }
	if ( i < particle_count ) {	
		p[i] = EOS(PresEq[i], Cs[i], P0[i],Density[i], RefDensity[i]); //CALL BEFORE!
	}
}

__global__ void StressStrainExtKernel(double *sigma,	//OUTPUT
																								double *strain,double *straina,double *strainb, //OUTPUT
																								//INPUT
																								double *p, double *rotrate, 
																								double *shearstress,double *shearstressa, double *shearstressb,
																								
																								double dt, int particle_count) {
	// int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	// if ( i < particle_count ) {	
		// //Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// // Jaumann rate terms
		// tensor3 RotationRateT,SRT,RS;
		// tensor3 RotationRate;
		
		// double temprr[6],tempss[6];
		// for (int k=0;k<6;k++){ //First the diagonal
			// temprr[k]=rotrate[6*i+k];
			// tempss[k]=shearstress[6*i+k];
		// }
		
		// RotationRate.FromFlatSym(tempss);
		// RotationRate.FromFlatAntiSym(temprr);
		// RotationRateT = temprr.Trans();
		// SRT = ;
		// Trans(RotationRate,RotationRateT);
		// Mult(ShearStress,RotationRateT,SRT);
		// Mult(RotationRate,ShearStress,RS);
		// double dep =0.;
		// double prev_sy;
		// double Et;
		
		// // Elastic prediction step (ShearStress_e n+1)
		// if (FirstStep)
			// ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
		// ShearStressb	= ShearStressa;
		// ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;	

		// //Fail, TODO
		
		// ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
		// Sigma = -p[i] * OrthoSys::I + ShearStress;	//Fraser, eq 3.32
		
		// if (FirstStep)
			// Straina	= -dt/2.0*StrainRate + Strain;
		// Strainb	= Straina;
		// Straina	= dt*StrainRate + Straina;
		// Strain	= 1.0/2.0*(Straina+Strainb);
	// }
}

__device__ void Domain_d::StressStrain(int i) {
	
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	double dep = 0.;
	
	if ( i < particle_count ) {	
		//Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// Jaumann rate terms
		tensor3 RotationRateT,SRT,RS;
		tensor3 RotationRate;
		tensor3 StrainRate;
		tensor3 ShearStress,ShearStressa,ShearStressb;
		tensor3 Sigma;
		tensor3 Strain,Straina,Strainb;
		
		double temprr[6],tempss[6],tempsr[6];
		double tempssa[6],tempssb[6];
		for (int k=0;k<6;k++){ //First the diagonal
			temprr[k] = rotrate[6*i+k];
			tempss[k] = shearstress[6*i+k];
			tempsr[k] = strrate[6*i+k];
			tempssa[k]= shearstressa[6*i+k];
			tempssb[k]= shearstressb[6*i+k];			
		}
		ShearStress   = FromFlatSym (tempss);
		ShearStressa  = FromFlatSym(tempssa);
		ShearStressb  = FromFlatSym(tempssb);
		
		StrainRate    = FromFlatSym(tempsr);
		RotationRate  = FromFlatAntiSym(temprr);

		
		RotationRateT = Trans(RotationRate);
		
		SRT = ShearStress * RotationRateT;
		RS = RotationRate * ShearStress;

		// if (i==1250){
			// printf("Stress Kernel, StrainRate\n");print(StrainRate);
			// printf("Stress Kernel, Identity() before calc\n");print((StrainRate.xx+StrainRate.yy+StrainRate.zz)* Identity());
			// printf("G, %f\n",G[i]);}
			
		// Elastic prediction step (ShearStress_e n+1)
		if (isfirst_step){
			ShearStressa	= -deltat/2.0*(2.0*G[i]*(StrainRate - 1.0/3.0*(StrainRate.xx+StrainRate.yy+StrainRate.zz)* Identity()) + SRT+RS) + ShearStress;
		}
		ShearStressb	= ShearStressa;
		ShearStressa	= deltat*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate.xx+StrainRate.yy+StrainRate.zz)*Identity())+SRT+RS) + ShearStressa;	

		// if (i==1250){
		// printf("Stress Kernel ShearStressA\n");print(ShearStressa);}
			
		// //Fail, TODO
				
		double J2	= 0.5*(ShearStressa.xx*ShearStressa.xx + 2.0*ShearStressa.xy*ShearStressa.yx +
					2.0*ShearStressa.xz*ShearStressa.zx + ShearStressa.yy*ShearStressa.yy +
					2.0*ShearStressa.yz*ShearStressa.zy + ShearStressa.zz*ShearStressa.zz);

    //Scale back, Fraser Eqn 3-53
		double sig_trial = sqrt(3.0*J2); 
    if ( sigma_y[i] < sig_trial ) ShearStressa = sigma_y[i]/sig_trial * ShearStressa; //Yielding      
    //std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;
		
		sigma_eq[i] = sig_trial;	
		
		if ( sig_trial > sigma_y[i]) {
			dep=( sig_trial - sigma_y[i])/ (3.*G[i] /*+ Ep*/);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			pl_strain[i] += dep;	
      //printf("Particle %d, dep %.1e, sigtrial %.1e\n",i,dep,sig_trial);
			sigma_eq[i] = sigma_y[i];
		}

    
		ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
		Sigma = -p[i] * Identity() + ShearStress;	//Fraser, eq 3.32
		// if (i == 1250){
			// printf("Time %.4e Particle 1250, pressure %f , ShearStresszz %f Sigma \n",Time, p[i], ShearStress.zz);
			// print(Sigma);
		// }
		
		if (isfirst_step)
			Straina	= -deltat/2.0*StrainRate + Strain;
		Strainb	= Straina;
		Straina	= deltat*StrainRate + Straina;
		Strain	= 1.0/2.0*(Straina+Strainb);

		///// OUTPUT TO Flatten arrays
		ToFlatSymPtr(Sigma, sigma,6*i);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM
		
		ToFlatSymPtr(Strain, 	strain,6*i);
		ToFlatSymPtr(Straina, straina,6*i);
		ToFlatSymPtr(Strainb, strainb,6*i);
		
		ToFlatSymPtr(ShearStress, shearstress,6*i);
		ToFlatSymPtr(ShearStressa, shearstressa,6*i);
		ToFlatSymPtr(ShearStressb, shearstressb,6*i);
		
		// if (i==1250){
			// printf("Stress Strain kernel, particle 1250 Sigma\n");
			// print(Sigma);
		// }
	}//particle count
}

//////// DO NOT USE A AND B, THIS IS USED BY KICKDRIFT /////////////
__device__ void Domain_d::StressStrainOne(int i) {
	
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	double dep = 0.;
	
	if ( i < solid_part_count ) {	
		//Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// Jaumann rate terms
		tensor3 RotationRateT,SRT,RS;
		tensor3 RotationRate;
		tensor3 StrainRate;
		tensor3 ShearStress;
		tensor3 Sigma;
		tensor3 Strain,Straina,Strainb;
		
		double temprr[6],tempss[6],tempsr[6];
		double tempssa[6],tempssb[6];
		for (int k=0;k<6;k++){ //First the diagonal
			temprr[k] = rotrate[6*i+k];
			tempss[k] = shearstress[6*i+k];
			tempsr[k] = strrate[6*i+k];
		}
		ShearStress   = FromFlatSym (tempss);	
		StrainRate    = FromFlatSym(tempsr);
		RotationRate  = FromFlatAntiSym(temprr);
   // printf(" Strain rate xx %f \n",StrainRate.zz);
		RotationRateT = Trans(RotationRate);
		
		SRT = ShearStress * RotationRateT;
		RS = RotationRate * ShearStress;

		ShearStress	= deltat*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate.xx+StrainRate.yy+StrainRate.zz)*Identity())+SRT+RS) + ShearStress;	

    eff_strain_rate[i] = sqrt ( 	0.5*( (StrainRate.xx-StrainRate.yy)*(StrainRate.xx-StrainRate.yy) +
                                (StrainRate.yy-StrainRate.zz)*(StrainRate.yy-StrainRate.zz) +
                                (StrainRate.zz-StrainRate.xx)*(StrainRate.zz-StrainRate.xx)) + 
                          3.0 * (StrainRate.xy*StrainRate.xy + StrainRate.yz*StrainRate.yz + StrainRate.zx*StrainRate.zx)
                        );
                        
		double J2	= 0.5*(ShearStress.xx*ShearStress.xx + 2.0*ShearStress.xy*ShearStress.yx +
					2.0*ShearStress.xz*ShearStress.zx + ShearStress.yy*ShearStress.yy +
					2.0*ShearStress.yz*ShearStress.zy + ShearStress.zz*ShearStress.zz);

    //Scale back, Fraser Eqn 3-53
		double sig_trial = sqrt(3.0*J2); 
    if ( sigma_y[i] < sig_trial ) ShearStress = sigma_y[i]/sig_trial * ShearStress; //Yielding      
    //std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;

    if       (mat[i]->Material_model == HOLLOMON )       sigma_y[i] = CalcHollomonYieldStress(pl_strain[i],mat[i]); 
	else if  (mat[i]->Material_model == JOHNSON_COOK )   sigma_y[i] = CalcJohnsonCookYieldStress(pl_strain[i],eff_strain_rate[i],T[i],mat[i]);
		
		sigma_eq[i] = sig_trial;	
		
		if ( sig_trial > sigma_y[i]) {
      if              (mat[i]->Material_model == HOLLOMON ) {
		double Et;
				Et =CalcHollomonTangentModulus(pl_strain[i], mat[i]); //Fraser 3.54
				// Et_m = Et;        
      } else if       (mat[i]->Material_model == JOHNSON_COOK ) {
        //Et = mat->CalcTangentModulus(pl_strain, eff_strain_rate, T); //Fraser 3.54
      } else if       (mat[i]->Material_model == BILINEAR ) {
        //Ep = mat->Elastic().E()*Et/(mat->Elastic().E()-Et);
      }
      //if (Ep<0) Ep = 1.*mat->Elastic().E();
      
			dep=( sig_trial - sigma_y[i])/ (3.*G[i] /*+ Ep*/);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			pl_strain[i] += dep;	
      //printf("Particle %d, dep %.1e, sigtrial %.1e\n",i,dep,sig_trial);
			sigma_eq[i] = sigma_y[i];
		}//plastic

    // if (mat[i]->Material_model == BILINEAR )
      // sigma_y[i] += dep*Ep;
    
		Sigma = -p[i] * Identity() + ShearStress;	//Fraser, eq 3.32

		Strain	= deltat * StrainRate + Strain;
    
    // if (mat[i]->Material_model==JOHNSON_COOK){
      // printf("JOHNSON_COOK!\n"); //test
    // }

		///// OUTPUT TO Flatten arrays
		ToFlatSymPtr(Sigma, sigma,6*i);  //TODO: CHECK IF RETURN VALUE IS SLOWER THAN PASS AS PARAM		
		ToFlatSymPtr(Strain, 	strain,6*i);		
		ToFlatSymPtr(ShearStress, shearstress,6*i);

	}//particle count
}

__global__ void StressStrainKernel(Domain_d *dom){
  printf("DEPRECATED. SSA AND SSB NOT ALLOCATED ANYMORE\n");
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom->StressStrain(i);
}

__global__ void StressStrainKickDriftKernel(Domain_d *dom){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom->StressStrainOne(i);
}

#define TAU		0.005
#define VMAX	10.0

// TODO #53 Make generic function pointer
// THISIS ONLY AN EXAMPLE	
__device__ void Domain_d::ApplyBCVel(int bcid, 
																		double3 bcv){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	
	if ( i < particle_count ) {	
		//printf("particle %d bc \n",i);
		if (ID[i]==bcid){
			a[i]		= make_double3(0.0);
			v[i]		= bcv;
			va[i]		= bcv;
			
		}
	}
}

__global__ void ApplyBCVelKernel (Domain_d *dom, int bcid, double3 bcv) {
	
	dom->ApplyBCVel (bcid,bcv);
}

__global__ void ApplyBCVelExtKernel(	double *v, //Output
																double *va,
																int *ID, 	//Input
																int bcid, 
																double bcv,
																double Time,
																int particle_count) {
	
	int i = threadIdx.x + blockDim.x*blockIdx.x;	
	
	//VMAX/TAU * domi.getTime();
	
	if ( i < particle_count ) {	
		//if (ID[i]==bcid)
			
	}
}

__global__ void TimestepCheckKernel(const double &CFL,
																double *h,
																double *Cs){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	
	
	//VMAX/TAU * domi.getTime();
	
	// if ( i < particle_count ){
		
		// t1 = CFL*h[i]/Cs[i];//Or is Cij??
	// }															

}

#include "cuda_helper.h"
#include "cuNSearch.h"

#include<array>
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;
using namespace cuNSearch;
using Real3 = std::array<Real, 3>;

__global__ void testNeighboursKernel(	const uint particle,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	if (i == particle){
		printf("Particle %d nbs\n", particle);
		const uint writeOffset = neighborWriteOffsets[i];

		for (int j=0; j< particlenbcount[i];j++){
			printf("%d ", neighbors[writeOffset + j]);
		}

	}
}



};//SPH