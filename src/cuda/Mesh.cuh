#ifndef _MESH_CUH_
#define _MESH_CUH_

#ifndef __CUDACC__
#include <cmath>
#endif

#include "cuda_runtime.h"

#include "NastranReader.h"


namespace SPH{

struct elem_data{
    double3             *centroid,*normal;
};
//Element is not anymore here, is everything flattened in mesh_d class
class TriMesh_d{
	
	public:

	//Element 						elem_data;
	int 								id;
	double3 						*node,*node_v; //Positions and veloc, 
	int									*elnode;			//3 per element
	double 							*pplane; //Only for testing
  //Element data, TODO: PASS TO ELEMDATA
  double3             *centroid,*normal;
  int                 *nfar;
  double3             m_v,m_w;
  int nodecount, elemcount;
  
	
	//double							v;						//Constant Uniform v
	TriMesh_d(){}
	void ReadFromNastran(NastranReader &nr, bool flipnormals = false);
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens);
  void __device__ __host__ SetVel(const double3 &v) {m_v = v;} //Like in WeldForm CPU version
  void __device__ __host__ SetID(int i) {id = i;} //Like in WeldForm CPU version
  void __device__ __host__ GetID(int *i) {*i=id;} //Like in WeldForm CPU version
	inline void ApplyConstVel(const double3 &v);
	inline void CalcCentroidVelFromNodes();
	inline __host__ __device__ void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline __device__ void Move(double dt);
  inline __host__ __device__ void CalcNormals();
	inline __host__ __device__ void CalcSpheres();
	inline __device__ void CalcCentroids();
  inline __device__ void CheckNormals();
  inline __host__ __device__ void Move(const double3 &v);
  
  
};

__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt);

__global__ inline void CalcSpheresKernel(TriMesh_d *mesh_d);

};//SPH

#include "Mesh.cu"

#endif
