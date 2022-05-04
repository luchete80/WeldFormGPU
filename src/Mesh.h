#ifndef _MESH_H_
#define _MESH_H_

#include "Vector.h"
#include <vector>

namespace SPH{
	
// class Sphere{
	// public:
	// Sphere();
	// Sphere (Element *elem);
	
	// Vector  pos;
	// //Vector* x[3];
	// Element *element;
	// int node[3];	//CONST 
	// double radius;
	
// };


// It is like triangular mesh
class Element{
	public:
	Element(){}
	Element(const int &n1, const int &n2, const int &n3);
	
	//SPHERE
	Vector 	centroid;	
	Vector 	normal;
	Vector 	v;					//At centroid
	double 	radius;
	int 		node[3];
	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
	int 		nfar;						//farthest away node from baricenter
	//Sphere* centroid;
	//Mesh*		mesh;
};


class TriMesh{
	
	public:

	std::vector <Element* > 	element;
	std::vector <Vector* > 		node;
	std::vector <Vector* > 		node_v;				//Node velocities
	
	Vector							v;						//Constant Uniform v
	TriMesh();
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const Vector p1, const Vector p2, const int &dens);
	inline void ApplyConstVel(const Vector &v);
	inline void CalcCentroidVelFromNodes();
	inline void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline void CalcNormals();
	inline void CalcSpheres();
	void CalcCentroids();
};

};
#include "Mesh.cpp"

#endif