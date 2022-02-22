#include "Domain_d.cuh"
#include "Functions.cuh"
#include "tensor.cuh"

namespace SPH {
__global__ void CalcForce2233(PartData_d *partdata){
	
	partdata->CalcForce2233();
}
//Be a part data member???
//CALLED BY GLOBAL
__device__ inline void PartData_d::CalcForce2233()
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	int neibcount;
	#ifdef FIXED_NBSIZE
	neibcount = neib_offs[i];
	#else
	neibcount =	neib_offs[i+1] - neib_offs[i];
	#endif
	printf("Solving\n");
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];
		int j = NEIB(i,k);
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = x[i] - partdata->x[j];
		double nxij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		
		//Artifficial visc
		// double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		// double Beta	= (P1->Beta + P2->Beta)/2.0;
		
		if (!IsFree[i]) {
			di = DensitySolid(PresEq[i], P2->Cs, P2->P0,p[i], P2->RefDensity);
			mi = FPMassC[i] * m[j];
		} else {
			di = rho[i];
			mi = m[i];
		}
		if (!IsFree[j]) {
			dj = DensitySolid (PresEq[i], P1->Cs, P1->P0,p[j], P1->RefDensity);
			mj = P2->FPMassC * P1->Mass;
		} else {
			dj = rho[j];
			mj = m[j];
		}

		double3 vij	= v[i] - v[j];
		double h_ = (h[i] + h[j])/2.0;
		double nxij = length(xij);
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, 0, nxij/h_, h_);
		double K	= Kernel(3, 0, rij/h, h);
		
		////// Artificial Viscosity
		tensor3 PIij;
		//set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!IsFree[i]) Ci = SoundSpeed(PresEq[j], Cs[j], di, RefDensity[j]); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!IsFree[j]) Cj = SoundSpeed(PresEq[j], Cs[i], dj, RefDensity[i]); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			Cij = 0.5*(Ci+Cj);
			
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * Identity();		///<(2.74) Li, Liu Book
		}
		
		
		tensor3 Sigmaj,Sigmai;
		// set_to_zero(Sigmaj);
		// set_to_zero(Sigmai);
		
		//TODO: CONVERT FLATTENED ARRAY TO TENSOR
		Sigmai = Sigma[i];
		Sigmaj = Sigma[j];

		//THIS IS COMMENTED IN THE ORIGINAL CODE
//		if (IsFree[i]) Sigmai = P1->Sigma; else  Sigmai = P2->Sigma;
//		if (IsFree[j]) Sigmaj = P2->Sigma; else  Sigmaj = P1->Sigma;
		
		// Tensile Instability //////////////////////
		tensor3 TIij;
		//set_to_zero(TIij);
		if (TI[i] > 0.0 || TI[j] > 0.0) 
			TIij = pow((K/Kernel(Dimension, KernelType, (TIInitDist[i] + TIInitDist[j])/(2.0*h_), h_)),(P1->TIn+P2->TIn)/2.0)*(TIR[i]+TIR[j]);
			//TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR); //COMMENTED IN ORIGINAL CODE
		
		// NoSlip BC velocity correction 		////////////////////////////////
		float3 vab = make_float3(0.0);
		if (IsFree[i]*IsFree[j]) {
			vab = vij;
		} else {
			if (P1->NoSlip || P2->NoSlip) {
				// No-Slip velocity correction
				if (IsFree[i])	vab = v[i] - (2.0f*v[j]-P2->NSv); else vab = (2.0f*v[i]-P1->NSv) - v[j];
			}
			// Please check
			if (!(P1->NoSlip || P2->NoSlip)) {
				if (IsFree[i]) vab = v[i] - v[j]b; else vab = v[i] - v[j];
//				if (IsFree[i]) vab.x = v[i](0) + v[j]b(0); else vab.x = -v[i]b(0) - v[j](0);
			}
		} //Are not both fixed
		
		tensor3 StrainRate,RotationRate;
		// set_to_zero(StrainRate);
		// set_to_zero(RotationRate);

		////////////////////////////////////
		// // Calculation strain rate tensor
		////////////////////////////////////
		StrainRate(0,0) = 2.0*vab.x*xij.x;
		StrainRate(0,1) = vab.x*xij.y+vab.y*xij.x;
		StrainRate(0,2) = vab.x*xij.z+vab.z*xij.x;
		StrainRate(1,0) = StrainRate(0,1);
		StrainRate(1,1) = 2.0*vab.y*xij.y;
		StrainRate(1,2) = vab.y*xij.z+vab.z*xij.y;
		StrainRate(2,0) = StrainRate(0,2);
		StrainRate(2,1) = StrainRate(1,2);
		StrainRate(2,2) = 2.0*vab.z*xij.z;
		StrainRate	= -0.5 * GK * StrainRate;

		// // Calculation rotation rate tensor
		RotationRate(0,1) = vab.x*xij.y-vab.y*xij.x;
		RotationRate(0,2) = vab.x*xij.z-vab.z*xij.x;
		RotationRate(1,2) = vab.y*xij.z-vab.z*xij.y;
		RotationRate(1,0) = -RotationRate(0,1);
		RotationRate(2,0) = -RotationRate(0,2);
		RotationRate(2,1) = -RotationRate(1,2);
		RotationRate	  = -0.5 * GK * RotationRate;

		// XSPH Monaghan
		if (XSPH != 0.0  && (IsFree[i]*IsFree[j])) {
			//omp_set_lock(&P1->my_lock);
			P1->VXSPH += XSPH*mj/(0.5f*(di+dj))*K*(-vij);
			//omp_unset_lock(&P1->my_lock);

			//omp_set_lock(&P2->my_lock);
			P2->VXSPH += XSPH*mi/
			(0.5*(di+dj))*
			K*vij;
			//omp_unset_lock(&P2->my_lock);
		}		
		
	}//neibcount
	
}

}; //SPH