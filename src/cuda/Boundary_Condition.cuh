#ifndef _BOUNDARY_CONDITION_CUH
#define _BOUNDARY_CONDITION_CUH

namespace SPH {
struct boundaryCondition {
	int 	zoneId;
	int 	type;	// ENUM TYPE Velocity, Force, Temperature
	bool 	free;	//is necessary??
	int 	valueType;		//0: Constant, 1 amplitude table
	double3 value;       //If constant
  double3 value_ang;       //Angular value
	int 	ampId;			//if valuetype == 1
	double 	ampFactor;		//if valuetype == 1
};

}; //SPH


#endif