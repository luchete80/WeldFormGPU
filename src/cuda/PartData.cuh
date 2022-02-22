#ifndef SPH_PARTDATA_CUH
#define SPH_PARTDATA_CUH

class PartData_d{
	
	public:
	//cuNSearch::NeighborhoodSearch neib;
	//Structured in AOS
	int **neib;	//array of lists
	int *neibs;	//1D array, faster
	int *neib_offs;	//1D array, faster
	int *neibcount;	//Useful??
	
	bool *isFree, *noSlip;
	
	//General
	double *rho, *m;	//Mass and density
	//THERMAL
	double *T, *Ta, *Tb, *dTdt;
	double *k_T, *cp_T,*h_conv, *T_inf;
	
	//Mechanical
	double *p;
	double *sigma; //To convert after to tensor;
	
	//Be in another class
	double  *FPMassC;        ///< Mass coefficient for fixed particles to avoid leaving particles
	
};

#endif