#include "Domain_d.cuh"
#include "Functions.cuh"
#include "tensor.cuh"

namespace SPH {
// THIS CAUSED UNRESOLVED EXTERNAL, WITH INLINE TENSOR FUNCTIONS 
__global__ void CalcForcesMember(PartData_d *partdata){
	
	partdata->CalcForce2233(0,0.0);
}

__global__ void CalcForcesKernel(Domain_d *dom_d){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;
	dom_d->CalcForce2233(0,0.0);
}

#define NEIBS(i, k) partdata->neib_part[partdata->neib_offs[i]+k]
__global__ void CalcForce2233(PartData_d *partdata){
	
	//partdata->CalcForce2233(0,0.0);
	CalcForcesExt(partdata);
}

// __device__ tensor3 Identity(){
	// tensor3 ret;
	// ret(0,0) = ret(1,1) = ret(2,2) = 1.;
	// //ret[1][1]=ret[2][2]=1.;
	
	// return ret;
// }

// __device__ tensor3::tensor3(){

// }


//THIS IS NOT A FUNCTION MEMBER
// THIS ALLOWS FOR COMPILING VIA NON SEPARABLE COMPILATION!!!!
__device__ void CalcForcesExt(PartData_d *partdata){
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	int Dimension = 3; //TODO, put in another 
	int neibcount;
	#ifdef FIXED_NBSIZE
	neibcount = partdata->neib_offs[i];
	#else
	neibcount =	partdata->neib_offs[i+1] - partdata->neib_offs[i];
	#endif
	printf("Solving\n");
	for (int k=0;k < neibcount; k++) { //Or size
		int j = NEIBS(i,k);
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = partdata->x[i] - partdata->x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
		if (!partdata->IsFree[i]) {
			di = DensitySolid(partdata->PresEq[i], partdata->Cs[j], partdata->P0[j],partdata->p[i], partdata->rho_0[j]);
			mi = partdata->FPMassC[i] * partdata->m[j];
		} else {
			di = partdata->rho[i];
			mi = partdata->m[i];
		}
		if (!partdata->IsFree[j]) {
			dj = DensitySolid (partdata->PresEq[i], partdata->Cs[i], partdata->P0[i],partdata->p[j], partdata->rho_0[i]);
			mj = partdata->FPMassC[j] * partdata->m[i];
		} else {
			dj = partdata->rho[j];
			mj = partdata->m[j];
		}	
		
		double3 vij	= partdata->v[i] - partdata->v[j];
		double h_ = (partdata->h[i] + partdata->h[j])/2.0;
			
		//double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double GK	= GradKernel(3, 0, rij/h_, h_);
		double K	= Kernel(3, 0, rij/h_, h_);

		////// Artificial Viscosity
		tensor3 PIij;
		//set_to_zero(PIij);
		if (partdata->Alpha!=0.0 || partdata->Beta!=0.0)
		{
			double MUij = h_*dot(vij,xij)/(rij*rij+0.01*h_*h_);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!partdata->IsFree[i]) Ci = SoundSpeed(partdata->PresEq[j], partdata->Cs[j], di, partdata->rho_0[j]); 
			else 											Ci = SoundSpeed(partdata->PresEq[i], partdata->Cs[i], di, partdata->rho_0[i]);
			if (!partdata->IsFree[j]) Cj = SoundSpeed(partdata->PresEq[j], partdata->Cs[i], dj, partdata->rho_0[i]); 
			else 											Cj = SoundSpeed(partdata->PresEq[j], partdata->Cs[j], dj, partdata->rho_0[j]);
			Cij = 0.5*(Ci+Cj);
			
			if (dot(vij,xij)<0) PIij = (partdata->Alpha*Cij*MUij + partdata->Beta*MUij*MUij)/(0.5*(di+dj)) * Identity();		///<(2.74) Li, Liu Book
		}
		
	}//neibcount	
}

//TODO; COMPARE WITH ORIGINAL IN PARTDATA
//This exludes thermal, 
//__global__ void CalcForce2233(PartMassVolInfo
																//PartMechData *pmd
																//){
	
//Be a part data member???
//CALLED BY GLOBAL
//TODO; DIVIDE PARTDATA INTO DIFFERENT FIELDS
__device__ /*inline*/ void PartData_d::CalcForce2233(
	/* const double &Dimension*/
	int KernelType,
	float XSPH)
{
}



__device__ /*inline*/ void Domain_d::CalcForce2233(
	/* const double &Dimension*/
	int KernelType,
	float XSPH)
{
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	
	if ( i < particle_count ) {
	int Dimension = 3; //TODO, put in another 
	int neibcount;
	#ifdef FIXED_NBSIZE
	neibcount = neib_offs[i];
	#else
	neibcount =	neib_offs[i+1] - neib_offs[i];
	#endif
	//printf("Solving\n");
	tensor3 StrainRate,RotationRate;
	
	for (int k=0;k < neibcount; k++) { //Or size
		//if fixed size i = part * NB + k
		//int j = neib[i][k];
		int j = NEIB(i,k);
		//double h	= partdata->h[i]+P2->h)/2;
		double3 xij = x[i] - x[j];
		double rij = length(xij);
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		
		//Artifficial visc
		// double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		// double Beta	= (P1->Beta + P2->Beta)/2.0;
		
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
		
		//printf("i %d, Ti %f\n",i, T[i]);
		
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
		
		Sigmai.FromFlatSym(tempi);
		Sigmaj.FromFlatSym(tempj);
		//Sigmai = Sigma[i];
		//Sigmaj = Sigma[j];

		//THIS IS COMMENTED IN THE ORIGINAL CODE
//		if (IsFree[i]) Sigmai = P1->Sigma; else  Sigmai = P2->Sigma;
//		if (IsFree[j]) Sigmaj = P2->Sigma; else  Sigmaj = P1->Sigma;
		
		// Tensile Instability //////////////////////
		tensor3 TIij;
		tensor3 TIRi, TIRj;
		// for (int k=0;k<6;k++) {
			// TIRi.FromFlatSymPtr(TIR[6*i+k]);
			// TIRj.FromFlatSymPtr(TIR[6*j+k]);
		// }
		//TODO: CONVERT TIR FROM FLATTENED ARRAY TO TENSOR
		//set_to_zero(TIij);
		if (TI[i] > 0.0 || TI[j] > 0.0) 
			TIij = pow((K/Kernel(Dimension, KernelType, (TIInitDist[i] + TIInitDist[j])/(2.0*h_), h_)),(TIn[i] + TIn[j])/2.0)*(TIRi+TIRj);
			//TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR); //COMMENTED IN ORIGINAL CODE
		
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
		
		// set_to_zero(StrainRate);
		// set_to_zero(RotationRate);

		////////////////////////////////////
		// // Calculation strain rate tensor
		////////////////////////////////////
		if (i==1250 || j==1250)
			printf("Time, i,j,vab: %.4e %d %d %f %f %f\n",Time, i,j,vab.x,vab.y,vab.z);
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
		
		if (i==1250 || j==1250)
			printf("Strain Rate i %f %f %f\n",StrainRate(0,0),StrainRate(1,1),StrainRate(2,2));

		// // Calculation rotation rate tensor
		RotationRate(0,1) = vab.x*xij.y-vab.y*xij.x;
		RotationRate(0,2) = vab.x*xij.z-vab.z*xij.x;
		RotationRate(1,2) = vab.y*xij.z-vab.z*xij.y;
		RotationRate(1,0) = -RotationRate(0,1);
		RotationRate(2,0) = -RotationRate(0,2);
		RotationRate(2,1) = -RotationRate(1,2);
		RotationRate	  = -0.5 * GK * RotationRate;
		
		//printf("Particle %d strain rate: %f %f %f\n",i,StrainRate(0,0),StrainRate(1,1),StrainRate(2,2));

		// XSPH Monaghan
		if (XSPH != 0.0  && (IsFree[i]*IsFree[j])) {
			//omp_set_lock(&P1->my_lock);
			VXSPH[i] += XSPH*mj/(0.5f*(di+dj))*K*(-vij);
			//omp_unset_lock(&P1->my_lock);
	
		
			//NOT WRITE IN THE OTHER PART!
			//omp_set_lock(&P2->my_lock);
			// VXSPH[j] += XSPH*mi/
			// (0.5*(di+dj))*
			// K*vij;
			//omp_unset_lock(&P2->my_lock);
		}		
		
		
		double3 temp = make_double3(0.0);
		double temp1 = 0.0;
		
		//if (GradientType == 0)
		if (i == 1250)
			printf("Particle 1250 Time %.4e, Sigmaizz %f , Sigmajzz %f\n",Time, Sigmai(2,2),Sigmaj(2,2));
		temp = ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij /*+ TIij */) * (GK*xij);
			//Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj /*+ PIij + TIij */) , temp); //TODO: TIR AND ARTIFF VISC
		// else
			// Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp);

		//if (Dimension == 2) temp(2) = 0.0;
		temp1 = dot( vij , GK*xij );

		// Locking the particle 1 for updating the properties
		a[i] 		+= mj * temp;
		drho[i]	+= mj * (di/dj) * temp1;

		if (IsFree[i]) {
			float mj_dj= mj/dj;
			//P1->ZWab	+= mj_dj* K;
			//printf("mj /dj %f\n",mj_dj);
			StrainRate = StrainRate + mj_dj * StrainRate;
			RotationRate = RotationRate + mj_dj * RotationRate;
			
			//P1->RotationRate = P1->RotationRate + mj_dj*RotationRate;
		}
		// else
			// P1->ZWab	= 1.0;

		// if (P1->Shepard)
			// if (P1->ShepardCounter == P1->ShepardStep)
				// P1->SumDen += mj*    K;

	
		// THIS IS THE ORIGINAL
		// Locking the particle 2 for updating the properties
		//omp_set_lock(&P2->my_lock);
			// P2->a		-= mi * temp;
			// P2->dDensity	+= mi * (dj/di) * temp1;
			// if (IsFree[j]) {
				// float mi_di = mi/di;
				// //P2->ZWab	+= mi_di* K;
				// //StrainRate[j]	  += StrainRate[j] + mi_di*StrainRate;
				// //RotationRate[j] += P2->RotationRate + mi_di*RotationRate;

			// }
			// else
				// P2->ZWab	= 1.0;

			// if (P2->Shepard)
				// if (P2->ShepardCounter == P2->ShepardStep)
					// P2->SumDen += mi*    K;

	
		
		}//neibcount

		///// OUTPUT TO Flatten arrays
		RotationRate.ToFlatSymPtr(rotrate,6*i);
		StrainRate.ToFlatSymPtr(strrate,6*i);	//Is the same for antisymm, stores upper diagonal
		if (i==1250)
			printf("TOTAL (SUM) Strain Rate part %d %f %f %f\n",i, StrainRate(0,0),StrainRate(1,1),StrainRate(2,2));
	}//i < partcount
}

}; //SPH