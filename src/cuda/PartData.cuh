#ifndef SPH_PARTDATA_CUH
#define SPH_PARTDATA_CUH

namespace SPH{
//#ifdef FIXED_NBSIZE //fixed nb per part (row), filled with zeroes
//#define MAXNB_PPART
//#define NEIB(i, k) neib_part [ MAXNB_PPART * i + k]  
//#else
#define NEIB(i, k) neib_part[neib_offs[i]+k]  //Just the right amount of indices, non filled with zeroes
//#endif

class PartData_d{
	
	public:
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neibs;	//1D array, faster
	int *neib_offs;	//1D array, faster
	int *neib_part;	//1D array, faster
	int *neibcount;	//Useful??
	
	bool *isFree, *noSlip;
	
	double3 *x,*v,*a;
	//General
	double *rho, *m;	//Mass and density

	double 	rho_0;	///< Reference Density of Particle

	//THERMAL
	double *T, *Ta, *Tb, *dTdt;
	double *k_T, *cp_T,*h_conv, *T_inf;
	
	//Mechanical
	double *p;
	double *sigma; //To convert after to tensor;
	
	//Be in another class
	double  *FPMassC;        ///< Mass coefficient for fixed particles to avoid leaving particles
		
	//ARTIFFICIAL VISC	
	double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
	double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
	
	//TENSILE INSTABILITY
	double	*TI;		///< Tensile instability factor
	double	*TIn;		///< Tensile instability power
	double 	*TIInitDist;	///< Initial distance of particles for calculation of tensile instability
	
	double3 VXSPH;
	
	//BOUNDARY
	bool *IsFree;
	
	
	__device__ inline void CalcForce2233();
	
};

__global__ void CalcForce2233(PartData_d *partdata);


};
#endif