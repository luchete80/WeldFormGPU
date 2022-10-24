#include "Domain_d.cuh"

namespace SPH{

__global__ inline void CalcAccelKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom_d->CalcAccel(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
	0,0.0);
}


__device__ inline void Domain_d::CalcAccel(
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors,
	int KernelType,
	float XSPH
	)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < solid_part_count ) {
	int Dimension = 3; //TODO, put in another 
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];
	
	//printf("Solving\n");
	tensor3 StrainRate, RotationRate;
	tensor3 StrainRateSum,RotationRateSum;
	
	a[i]		=	make_double3(0.,0.,0.);
	
	clear(RotationRateSum);
	clear(StrainRateSum);
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];
		int j = neighbors[writeOffset + k];
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		
		if (!IsFree[i]) {
			di = DensitySolid(PresEq[i], Cs[j], P0[j],p[i], rho_0[j]);
			mi = FPMassC[i] * m[j];
		} else {
			di = rho[i];
			mi = m[i];
		}
		if (!IsFree[j]) {
			dj = DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
			mj = FPMassC[j] * m[i];
		} else {
			dj = rho[j];
			mj = m[j];
		}

		double3 vij	= v[i] - v[j];
		double h_ = (h[i] + h[j])/2.0;
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, KernelType, rij/h_, h_);
		double K	= Kernel(3, 0, rij/h_, h_);
		
		////// Artificial Viscosity
		tensor3 PIij;
    clear(PIij);
		//set_to_zero(PIij);

		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h_*dot(vij,xij)/(rij*rij+0.01*h_*h_);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!IsFree[i]) Ci = SoundSpeed(PresEq[j], Cs[j], di, rho_0[j]); else Ci = SoundSpeed(PresEq[i], Cs[i], di, rho_0[i]);
			if (!IsFree[j]) Cj = SoundSpeed(PresEq[j], Cs[i], dj, rho_0[i]); else Cj = SoundSpeed(PresEq[j], Cs[j], dj, rho_0[j]);
			Cij = 0.5*(Ci+Cj);
			
			//printf("C %f %f\n",Ci,Cj);
			if (dot(vij,xij)<0) 
				PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * Identity();		///<(2.74) Li, Liu Book
		}
    
		tensor3 Sigma,Sigmaj,Sigmai;
	
		//TODO: CONVERT FLATTENED ARRAY TO TENSOR
		//TODO: Avoid temp array conversion and test
		double tempi[6],tempj[6];
		for (int k=0;k<6;k++){ //First the diagonal
			tempi[k]=sigma[6*i+k];
			tempj[k]=sigma[6*j+k];
		}
		
		Sigmai = FromFlatSym(tempi);
		Sigmaj = FromFlatSym(tempj);
		
		// Tensile Instability //////////////////////
		tensor3 TIij;
		tensor3 TIRi, TIRj;

		// NoSlip BC velocity correction 		////////////////////////////////
		double3 vab = make_double3(0.0);
		if (IsFree[i]*IsFree[j]) {
			vab = vij;
		} else {
			if (NoSlip[i] || NoSlip[j] ) {
				// No-Slip velocity correction
				if (IsFree[i])	vab = v[i] - (2.0f*v[j]- NSv[j]); 
				else vab = (2.0f*v[i]- NSv[i]) - v[j];
			}
			// Please check
			if (!(NoSlip[i] || NoSlip[j])) {
				if (IsFree[i]) vab = v[i] - v[j]; else vab = v[i] - v[j];
//				if (IsFree[i]) vab.x = v[i](0) + v[j]b(0); else vab.x = -v[i]b(0) - v[j](0);
			}
		} //Are not both fixed

		// XSPH Monaghan
		if (XSPH != 0.0  && (IsFree[i]*IsFree[j])) {
			//omp_set_lock(&P1->my_lock);
			VXSPH[i] += XSPH*mj/(0.5f*(di+dj))*K*(-vij);
		}		

		// Locking the particle 1 for updating the properties
		a[i] 		+= mj * ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij /* TIij */) * (GK*xij);
		}//neibcount

	}//i < partcount
}


__global__ inline void CalcDensIncKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom_d->CalcDensInc(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
	0);
}


__device__ inline void Domain_d::CalcDensInc(
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors,
	int KernelType
	)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < solid_part_count ) {
	int Dimension = 3; //TODO, put in another 
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];
	
	drho[i]	= 0.0;
	
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];
		int j = neighbors[writeOffset + k];
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		
		if (!IsFree[i]) {
			di = DensitySolid(PresEq[i], Cs[j], P0[j],p[i], rho_0[j]);
			mi = FPMassC[i] * m[j];
		} else {
			di = rho[i];
			mi = m[i];
		}
		if (!IsFree[j]) {
			dj = DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
			mj = FPMassC[j] * m[i];
		} else {
			dj = rho[j];
			mj = m[j];
		}

		double3 vij	= v[i] - v[j];
		double h_ = (h[i] + h[j])/2.0;
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, KernelType, rij/h_, h_);
		double K	= Kernel(3, 0, rij/h_, h_);
    
		// NoSlip BC velocity correction 		////////////////////////////////
		double3 vab = make_double3(0.0);
		if (IsFree[i]*IsFree[j]) {
			vab = vij;
		} else {
			if (NoSlip[i] || NoSlip[j] ) {
				// No-Slip velocity correction
				if (IsFree[i])	vab = v[i] - (2.0f*v[j]- NSv[j]); 
				else vab = (2.0f*v[i]- NSv[i]) - v[j];
			}
			// Please check
			if (!(NoSlip[i] || NoSlip[j])) {
				if (IsFree[i]) vab = v[i] - v[j]; else vab = v[i] - v[j];
//				if (IsFree[i]) vab.x = v[i](0) + v[j]b(0); else vab.x = -v[i]b(0) - v[j](0);
			}
		} //Are not both fixed

		double temp1 = 0.0;
    
		//if (Dimension == 2) temp(2) = 0.0;
		temp1 = dot( vij , GK*xij );
		drho[i]	+= mj * (di/dj) * temp1;
		}//neibcount

	}//i < partcount
}

///////////////////////////////////
//////// CalcRateTensors: FOR KICKDRIFT SOLVER WHERE TENSOR ARE CALCULATED AFTER 

__global__ inline void CalcRateTensorsKernel(Domain_d *dom_d,
	const uint *particlenbcount,
	const uint *neighborWriteOffsets,
	const uint *neighbors){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom_d->CalcRateTensors(
	particlenbcount,
	neighborWriteOffsets,
	neighbors);
}

__device__ /*__forceinline__*/inline void Domain_d::CalcRateTensors(const uint *particlenbcount,
                                                        const uint *neighborWriteOffsets,
                                                        const uint *neighbors){
                                                          
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < solid_part_count ) {
	int Dimension = 3; //TODO, put in another 
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];
  
	tensor3 StrainRate,RotationRate;
	tensor3 StrainRateSum,RotationRateSum;
	
	clear(StrainRateSum);
	clear(RotationRateSum);
	
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
    int j = neighbors[writeOffset + k];

		double3 xij = x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;

		if (!IsFree[i]) {
			di = DensitySolid(PresEq[i], Cs[j], P0[j],p[i], rho_0[j]);
			mi = FPMassC[i] * m[j];
		} else {
			di = rho[i];
			mi = m[i];
		}
		if (!IsFree[j]) {
			dj = DensitySolid (PresEq[i], Cs[i], P0[i],p[j], rho_0[i]);
			mj = FPMassC[j] * m[i];
		} else {
			dj = rho[j];
			mj = m[j];
		}

		double3 vij	= v[i] - v[j];
		double h_ = (h[i] + h[j])/2.0;
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, 0/*KernelType*/, rij/h_, h_);
		double K	= Kernel(3, 0, rij/h_, h_);
		
		tensor3 Sigma,Sigmaj,Sigmai;
		// set_to_zero(Sigmaj);
		// set_to_zero(Sigmai);
		
		//TODO: CONVERT FLATTENED ARRAY TO TENSOR
		//TODO: Avoid temp array conversion and test
		double tempi[6],tempj[6];
		for (int k=0;k<6;k++){ //First the diagonal
			tempi[k]=sigma[6*i+k];
			tempj[k]=sigma[6*j+k];
		}
		
		Sigmai = FromFlatSym(tempi);
		Sigmaj = FromFlatSym(tempj);

		double3 vab = make_double3(0.0);
		//if (IsFree[i]*IsFree[j]) {
			vab = vij;

		////////////////////////////////////
		// // Calculation strain rate tensor
		////////////////////////////////////
		StrainRate.xx = 2.0*vab.x*xij.x;
		StrainRate.xy = vab.x*xij.y+vab.y*xij.x;
		StrainRate.xz = vab.x*xij.z+vab.z*xij.x;
		StrainRate.yx = StrainRate.xy;
		StrainRate.yy = 2.0*vab.y*xij.y;
		StrainRate.yz = vab.y*xij.z+vab.z*xij.y;
		StrainRate.zx = StrainRate.xz;
		StrainRate.zy = StrainRate.yz;
		StrainRate.zz = 2.0*vab.z*xij.z;
		//StrainRate	= (-0.5) * GK * StrainRate;
		StrainRate	= StrainRate * ((-0.5) * GK);
		
		clear(RotationRate);
		// // Calculation rotation rate tensor
		RotationRate.xy = vab.x*xij.y-vab.y*xij.x;
		RotationRate.xz = vab.x*xij.z-vab.z*xij.x;
		RotationRate.yz = vab.y*xij.z-vab.z*xij.y;
		RotationRate.yx = -RotationRate.xy;
		RotationRate.xz = -RotationRate.xz;
		RotationRate.zy = -RotationRate.yz;
		//RotationRate	  	= -0.5 * GK * RotationRate; //THIS OPERATOR FAILS
		RotationRate	  	= RotationRate * (-0.5 * GK);

		double3 temp = make_double3(0.0);
		double temp1 = 0.0;
		
			double mj_dj= mj/dj;

			StrainRateSum 	= StrainRateSum + mj_dj * StrainRate;
			RotationRateSum = RotationRateSum + mj_dj * RotationRate;
      
		}//neibcount
		///// OUTPUT TO Flatten arrays
		ToFlatSymPtr(RotationRateSum, rotrate,6*i);
		ToFlatSymPtr(StrainRateSum, strrate,6*i);	//Is the same for antisymm, stores upper diagonal

	}//i < partcount                                                          
                                                          
}



}; //SPH