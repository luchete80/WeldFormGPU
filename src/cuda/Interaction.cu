#include "Domain_d.cuh"
#include "Functions.cuh"

//Be a part data member???
//CALLED BY GLOBAL
__device__ inline void ParticleData_d::CalcForce2233(ParticleData_d *partdata)
{
	int i = threadIdx.x+blockDim.x*blockIdx.x;
	
	int neibcount;
	#ifdef FIXED_NBSIZE
	neibcount = neib_offs[i];
	#else
	neibcount =	neib_offs[i+1] - neib_offs[i];
	#endif
	printf("Solving\n");
	for (int k=0;k < neibcount;k++) { //Or size
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
		// if (!IsFree[j]) {
			// dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
			// mj = P2->FPMassC * P1->Mass;
		// } else {
			// dj = P2->Density;
			// mj = P2->Mass;
		// }
	}
}