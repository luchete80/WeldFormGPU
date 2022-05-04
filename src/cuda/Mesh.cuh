#ifndef _MESH_CUH_
#define _MESH_CUH_

namespace SPH{

class Element_d{
	public:
	Element(){}
	Element(const int &n1, const int &n2, const int &n3);
	
	//SPHERE
	double3 *centroid;	
	double3 *normal;
	double3 *v;					//At centroid
	double 	radius;
	int 		*node;		//3 per element
	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
	int 		nfar;						//farthest away node from baricenter
	//Sphere* centroid;
	//Mesh*		mesh;
};

class TriMesh_d{
	
	public:

	Element 						elem_data;
	double3 						*node,node_v; //Positions and veloc
	double3 						*node_v;				//Node velocities
	
	double							v;						//Constant Uniform v
	TriMesh();
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const Vec3_t p2, const int &dens);
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