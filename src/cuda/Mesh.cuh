#ifndef _MESH_CUH_
#define _MESH_CUH_

namespace SPH{

//Element is not anymore here, is everything flattened in mesh_d class
class TriMesh_d{
	
	public:

	//Element 						elem_data;
	double3 						*node,node_v; //Positions and veloc, 
	int									*elnode;			//3 per element
	double 							*pplane;
	
	double							v;						//Constant Uniform v
	TriMesh_d();
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 &p2, const int &dens);
	inline void ApplyConstVel(const double3 &v);
	inline void CalcCentroidVelFromNodes();
	inline void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline void CalcNormals();
	inline void CalcSpheres();
	void CalcCentroids();
};
};

#endif