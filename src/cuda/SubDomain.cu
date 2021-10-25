#include "SubDomain.cuh"
#include "vector_math.h"

namespace SPH {

inline void SubDomain::DelParticles (int const & Tags)
{
    // Array<int> idxs; // indices to be deleted

	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<particlecount; i++)	//Like in Domain::Move
	// #else
	// for (int i=0; i<particlecount; i++)//Like in Domain::Move
	// #endif
    // {
        // if (Particles[i]->ID==Tags)
		// {
			// //omp_set_lock(&dom_lock);
        	// idxs.Push(i);
			// //omp_unset_lock(&dom_lock);
		// }
    // }
    // if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particles to delete");
    // Particles.DelItems (idxs);

    // std::cout << "\n" << "Particle(s) with Tag No. " << Tags << " has been deleted" << std::endl;
}



inline void __device__ SubDomain::StartAcceleration (float3 const & a) {

	// #pragma omp parallel for schedule(static) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<particlecount; i++)	//Like in Domain::Move
	// #else
	// for (int i=0; i<particlecount; i++)//Like in Domain::Move
	// #endif
	for (int i=0; i<particlecount;i++) {
	    if (Particles[i]->IsFree){
			// Tensile Instability for all soil and solid particles
			if (Particles[i]->TI > 0.0)
        		{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2)
				{
					double teta, Sigmaxx, Sigmayy, C, S;

					if ((Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[i]->Sigma(0,1)/(Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[i]->Sigma(0,0) + 2.0*C*S*Particles[i]->Sigma(0,1) + S*S*Particles[i]->Sigma(1,1);
					Sigmayy = S*S*Particles[i]->Sigma(0,0) - 2.0*C*S*Particles[i]->Sigma(0,1) + C*C*Particles[i]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[i]->TI * Sigmaxx/(Particles[i]->Density*Particles[i]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[i]->TI * Sigmayy/(Particles[i]->Density*Particles[i]->Density); else Sigmayy = 0.0;
					Particles[i]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[i]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[i]->TIR(0,1) = Particles[i]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else
				{
					tensor3 Vec,Val,VecT,temp;
					double pc_ti_inv_d2=Particles[i]->TI/(Particles[i]->Density*Particles[i]->Density);//Precompute some values
					Rotation(Particles[i]->Sigma,Vec,VecT,Val);
					//Before
					// if (Val(0,0)>0) Val(0,0) = -Particles[i]->TI * Val(0,0)/(Particles[i]->Density*Particles[i]->Density); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -Particles[i]->TI * Val(1,1)/(Particles[i]->Density*Particles[i]->Density); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -Particles[i]->TI * Val(2,2)/(Particles[i]->Density*Particles[i]->Density); else Val(2,2) = 0.0;
					if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					//Mult(Vec,Val,temp);
					temp=Vec*Val;
					//Mult(temp,VecT,Particles[i]->TIR);
					Particles[i]->TIR = temp,VecT;
				}
			}
	    	}
	    	else
	    	{
	       		// Reset the pressure and the induced velocity for solid boundaries
	    		Particles[i]->NSv =make_float3(0.0,0.0,0.0);
	    		Particles[i]->Pressure = 0.0;
				Particles[i]->Sigma = 0.;
				Particles[i]->Sigma=0.0;
				Particles[i]->ShearStress=0.0;
	    	}



		//Reset to zero for all particles
		Particles[i]->a		= a;
		Particles[i]->SatCheck	= false;
		Particles[i]->dDensity	= 0.0;
		Particles[i]->VXSPH	= make_float3(0.0);
		Particles[i]->ZWab	= 0.0;
		Particles[i]->SumDen	= 0.0;
		Particles[i]->SumKernel	= 0.0;
		if (Dimension == 2) Particles[i]->v.z =0.0;
		Particles[i]->StrainRate=0.0f;
		Particles[i]->RotationRate=0.0;
		
	}
}

inline __device__ void SubDomain::PrimaryComputeAcceleration () {
	size_t P1,P2;
	float3 xij;
	float h,K;		//TODO: cHANGE Change to double
	// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
	for (size_t a=0; a<FSMPairscount;a++) {
		P1	= FSMPairs[a][0];
		P2	= FSMPairs[a][1];
		xij	= Particles[P1]->x-Particles[P2]->x;
		h	= (Particles[P1]->h+Particles[P2]->h)/2.0;

		Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);

		K	= Kernel(Dimension, KernelType, length(xij)/h, h);
		if ( !Particles[P1]->IsFree ) {
			//omp_set_lock(&Particles[P1]->my_lock);
									Particles[P1]->SumKernel+= K;
				Particles[P1]->Pressure	+= Particles[P2]->Pressure * K /*+ dot(Gravity,xij)*Particles[P2]->Density*K*/;
				Particles[P1]->Sigma 	 = Particles[P1]->Sigma + K * Particles[P2]->Sigma;
				if (Particles[P1]->NoSlip)		Particles[P1]->NSv 	+= K * Particles[P2]->v;
			//omp_unset_lock(&Particles[P1]->my_lock);
		} else {
			//omp_set_lock(&Particles[P2]->my_lock);
									Particles[P2]->SumKernel+= K;
				Particles[P2]->Pressure	+= Particles[P1]->Pressure * K /*+ dot(Gravity,xij)*Particles[P1]->Density*K*/;
				Particles[P2]->Sigma	 = Particles[P2]->Sigma + K * Particles[P1]->Sigma;
				if (Particles[P2]->NoSlip)		Particles[P2]->NSv 	+= K * Particles[P1]->v;
			//omp_unset_lock(&Particles[P2]->my_lock);
		}
	}
	

	// Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
	for (int i=0; i<FixedParticlescount; i++){

		if (Particles[FixedParticles[i]]->SumKernel!= 0.0) {
			size_t a = FixedParticles[i];
			Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
			Particles[a]->Sigma	= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
			if (Particles[a]->NoSlip)	Particles[a]->NSv	= Particles[a]->NSv/Particles[a]->SumKernel;

			// Tensile Instability for fixed soil and solid particles
			if (Particles[a]->TI > 0.0)
			{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2) {
					double teta, Sigmaxx, Sigmayy, C, S;
					if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
					Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
					Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else {
					tensor3 Vec,Val,VecT,temp;
					Rotation(Particles[a]->Sigma,Vec,VecT,Val);
					double pc_ti_inv_d2=Particles[a]->TI/(Particles[a]->Density*Particles[a]->Density);//Precompute some values
					// if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
					if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					// Mult(Vec,Val,temp);
					// Mult(temp,VecT,Particles[a]->TIR);
					temp=Vec*Val;
					Particles[a]->TIR = temp * VecT;
				}
			}
		}
	}

}
/*
inline __device__ void SubDomain::Periodic_X_Correction(float3 & x, double const & h, Particle * P1, Particle * P2)
{
	if (Domsize(0)>0.0) {if (x(0)>2*Cellfac*h || x(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? x(0) -= Domsize(0) : x(0) += Domsize(0);}}
	if (Domsize(1)>0.0) {if (x(1)>2*Cellfac*h || x(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? x(1) -= Domsize(1) : x(1) += Domsize(1);}}
	if (Domsize(2)>0.0) {if (x(2)>2*Cellfac*h || x(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? x(2) -= Domsize(2) : x(2) += Domsize(2);}}
}
*/


inline void __device__ SubDomain::LastComputeAcceleration () {
	for (int i=0; i<SMPairscount;i++)
		CalcForce2233(Particles[SMPairs[i][0]],Particles[SMPairs[i][1]]);

	for (int i=0; i<FSMPairscount;i++)
		CalcForce2233(Particles[FSMPairs[i][1]],Particles[FSMPairs[i][1]]);


	//LUCIANO: THIS SHOULD BE PERFORMED OUTSIDE
	// for (int i=0 ; i<Nproc ; i++)
	// {
		// SMPairs[i].Clear();
		// FSMPairs[i].Clear();
		// NSMPairs[i].Clear();
	// }

	//Min time step check based on the acceleration
	double test	= 0.0;
	deltatmin	= deltatint;

	for (int i=0; i<particlecount; i++) {
		if (Particles[i]->IsFree) {
			test = sqrt(Particles[i]->h/length(Particles[i]->a));
			if (deltatmin > (sqrt_h_a*test)) {
				//omp_set_lock(&dom_lock);
					deltatmin = sqrt_h_a*test;
				//omp_unset_lock(&dom_lock);
			}
		}
	}
}



__global__ void StartAccelerationKernel(SubDomain &sd){


}


/*inline */ __host__ void StartAcceleration (SubDomain &sd) {

	StartAccelerationKernel<<<3,4 >>>(sd);

}

/*inline*/ __host__ void PrimaryComputeAcceleration(SubDomain &sd){

}
__global__ void PrimaryComputeAccelerationKernel(SubDomain &sd){

}


/*inline*/ __host__ void LastComputeAcceleration(SubDomain &sd){

} 
__global__ void LastComputeAccelerationKernel(SubDomain &sd){

}



//New, for Bonet gradient correction
// inline void Domain::CalcGradCorrMatrix () {
	// double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	// std::tensor3 < tensor3> temp(particlecount);
	// tensor3 m,mt;

	// //#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	// for ( size_t k = 0; k < Nproc ; k++) {
		// Particle *P1,*P2;
		// tensor3 xij;
		// double h,GK;
		// //TODO: DO THE LOCK PARALLEL THING
		// for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			// P1	= Particles[SMPairs[k][a].first];
			// P2	= Particles[SMPairs[k][a].second];
			// xij	= P1->x - P2->x;
			// h	= (P1->h+P2->h)/2.0;
			// GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			// di = P1->Density; mi = P1->Mass;
			// dj = P2->Density; mj = P2->Mass;
		
			// Dyad (tensor3(GK*xij),xij,m);
			// mt = mj/dj * m;
			// ////omp_set_lock(&P1->my_lock);
			// temp[SMPairs[k][a].first] = temp[SMPairs[k][a].first]  + mt;
			// temp[SMPairs[k][a].second]= temp[SMPairs[k][a].second] - mt;
		// }
	// }//Nproc

	// #pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	// for (int i=0; i<particlecount; i++){
		// //cout << "temp "<<temp[i]<<endl;
		// /** Inverse.*/
		// //inline void Inv (tensor3 const & M, tensor3 & Mi, double Tol=1.0e-10)}	
		// Inv(temp[i],m);		
		// Particles[i] ->gradCorrM = m;
	// }
	
// }

// inline void SubDomain::ClearNbData(){
	
	// for (int i=0 ; i<Nproc ; i++) { //In the original version this was calculated after
		// SMPairs[i].Clear();
		// FSMPairs[i].Clear();
		// NSMPairs[i].Clear();
	// }
	// CellReset();
	// ListGenerate();
	// m_isNbDataCleared = true;
// }

// FOR THE MOMENT IS SET IN CPU
// inline void __device__ SubDomain::WholeVelocity() {
    // //Apply a constant velocity to all particles in the initial time step
    // if (BC.allv.norm()>0.0 || BC.allDensity>0.0) {
    	// float3 vel = make_float3(0.0,0.0,0.0);
    	// double den = 0.0;

	// //#pragma omp parallel for schedule (static) private(vel,den) num_threads(Nproc)
    	// for (int i=0 ; i<particlecount ; i++) {
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
// }

};//SPH
