#include "Domain_d.cuh"
#include "Functions.cuh"
#include <iostream>
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
	}
	rhob[i] = rho[i];
	rhoa[i] += dt*drho[i];
	rho[i] = (rhoa[i]+rhob[i])/2.0;
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
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
	if ( i < particle_count ) {	
		//Pressure = EOS(PresEq, Cs, P0,Density, RefDensity); //CALL BEFORE!

		// Jaumann rate terms
		tensor3 RotationRateT,SRT,RS;
		tensor3 RotationRate;
		
		double temprr[6],tempss[6];
		for (int k=0;k<6;k++){ //First the diagonal
			temprr[k]=rotrate[6*i+k];
			tempss[k]=shearstress[6*i+k];
		}
		
		RotationRate.FromFlatSym(tempss);
		RotationRate.FromFlatAntiSym(temprr);
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
	}
}

__device__ void Domain_d::StressStrain() {
	
	int i = threadIdx.x + blockDim.x*blockIdx.x;
		
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
		for (int k=0;k<6;k++){ //First the diagonal
			temprr[k]=rotrate[6*i+k];
			tempss[k]=shearstress[6*i+k];
			tempsr[k]=strrate[6*i+k];
		}
		ShearStress.FromFlatSym(tempss);
		StrainRate.FromFlatSym(tempsr);
		RotationRate.FromFlatAntiSym(temprr);
		StrainRate.FromFlatSym(tempsr);
		
		RotationRateT = RotationRate.Trans();
		
		SRT = ShearStress * RotationRateT;
		RS = RotationRate * ShearStress;
		
		// Elastic prediction step (ShearStress_e n+1)
		if (isfirst_step)
			ShearStressa	= -deltat/2.0*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))* Identity())+SRT+RS) + ShearStress;
		ShearStressb	= ShearStressa;
		ShearStressa	= deltat*(2.0*G[i]*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*Identity())+SRT+RS) + ShearStressa;	

		// //Fail, TODO
		
		ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
		Sigma = -p[i] * Identity() + ShearStress;	//Fraser, eq 3.32
		
		if (isfirst_step)
			Straina	= -deltat/2.0*StrainRate + Strain;
		Strainb	= Straina;
		Straina	= deltat*StrainRate + Straina;
		Strain	= 1.0/2.0*(Straina+Strainb);
		
		///// OUTPUT TO Flatten arrays
		Sigma.ToFlatSymPtr(sigma,i);
		Strain.ToFlatSymPtr(strain,i);
		Straina.ToFlatSymPtr(straina,i);
		Strainb.ToFlatSymPtr(strainb,i);
		ShearStress.ToFlatSymPtr(shearstress,i);
		ShearStressa.ToFlatSymPtr(shearstressa,i);
		ShearStressb.ToFlatSymPtr(shearstressb,i);
		
	}//particle count
}

__global__ void StressStrainKernel(Domain_d *dom){
	dom->StressStrain();
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

void Domain_d::MechSolve(const double &tf){

	int N = particle_count;
	int threadsPerBlock = 256; //Or BlockSize
	int blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
  Time =0.;
	
	isfirst_step =true;

	while (Time<tf) {
		
		//This was in Original LastCompAcceleration
		CalcForcesKernel	<<<blocksPerGrid,threadsPerBlock >>>(this);
		cudaDeviceSynchronize(); //REQUIRED!!!!
		
		//IMPOSE BC!
		ApplyBCVelKernel	<<<blocksPerGrid,threadsPerBlock >>>(this, 2, make_double3(0.,0.,0.));
		cudaDeviceSynchronize();
		double vbc = VMAX/TAU*Time;
		
		ApplyBCVelKernel	<<<blocksPerGrid,threadsPerBlock >>>(this, 3, make_double3(0.,0.,-vbc));
		cudaDeviceSynchronize();

		//Move particle and then calculate streses and strains ()
		MoveKernelExt<<<blocksPerGrid,threadsPerBlock >>> (v, va,vb,
														rho, rhoa, rhob, drho,
														x, a,
														u, /*Mat3_t I, */deltat,
														isfirst_step, particle_count);	
		cudaDeviceSynchronize(); //REQUIRED!!!!
		

		//If kernel is the external, calculate pressure
		//Calculate pressure!
		PressureKernelExt<<<blocksPerGrid,threadsPerBlock >>>(p,PresEq,Cs,P0,rho,rho_0,particle_count);
		cudaDeviceSynchronize();
		// StressStrainExtKernel(sigma,	//OUTPUT
																									// double *strain,*straina,*strainb, //OUTPUT
																									// //INPUT
																									// double *p, double *rotrate, 
																									// double* shearstress,double* shearstressa, double* shearstressb,
												
																									// double dt, int particle_count);
		StressStrainKernel<<<blocksPerGrid,threadsPerBlock >>>(this);
		cudaDeviceSynchronize();
		if (isfirst_step) isfirst_step = false;
		Time +=deltat;
		
		cudaMemcpy(u_h, u, sizeof(double3) * particle_count, cudaMemcpyDeviceToHost);	
		double3 max= make_double3(0.,0.,0.);
		for (int i=0;i<particle_count;i++){
			if (u_h[i].x>max.x) max.x = u_h[i].x;
			if (u_h[i].y>max.y) max.y = u_h[i].y;
			if (u_h[i].z>max.z) max.z = u_h[i].z;
		}
		cout << "Max disp "<< max.x<<", "<<max.y<<", "<<max.z<<endl;
		
		cout << "Time "<<Time<<endl;
		
		//TODO: Pass toPartData
		//CalcForcesMember	<<<blocksPerGrid,threadsPerBlock >>>(partdata);
		//MechSolveKernel<<< >>>();
	
	}//while <tf
}

};//SPH