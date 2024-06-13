#include "Mesh.cuh"

namespace SPH {

void TriMesh_d::ReadFromNastran(NastranReader &nr, bool flipnormals){
  
  
  nodecount = nr.node_count;
  elemcount = nr.elem_count;
  
  cout << "Creating nodes.."<<endl;

  cudaMalloc((void **)&node   , 	nodecount * sizeof (double3));
  cudaMalloc((void **)&node_v , 	nodecount * sizeof (double3));
  
  double3 *node_h, *node_vh;
  node_h  =  new double3 [nodecount];
  node_vh =  new double3 [nodecount];
  
//  dimension = nr.dim;
//  //Insert nodes
 for (int n=0;n<nr.node_count;n++){

    // if (!flipnormals)
      // node.Push(new Vec3_t(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]));
    // else 
      // node.Push(new Vec3_t(nr.node[3*n+1],nr.node[3*n],nr.node[3*n+2]));

  node_h[n]		=make_double3(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]);
  node_vh[n]	=make_double3(0.,0.,0.);
 }
  cudaMemcpy(node, node_h,    nodecount * sizeof (double3), cudaMemcpyHostToDevice);
  cudaMemcpy(node_v, node_vh, nodecount * sizeof (double3), cudaMemcpyHostToDevice);  

  /////////////////// ELEMENTS ////////////////////
  cudaMalloc((void **)&centroid , 	elemcount * sizeof (double3));
  cudaMalloc((void **)&normal 	, 	elemcount * sizeof (double3));
  cudaMalloc((void **)&elnode 	, 	3 * elemcount * sizeof (int));
	
  int *elnode_h = new int[3*elemcount];
  double3 *centroid_h = new double3[elemcount];
  double3 *normal_h   = new double3[elemcount];
//  cout << "Generated "<<node.Size()<< " trimesh nodes. "<<endl;
//  //cout << "Normals"<<endl;
//  cout << "Writing elements..."<<endl;
 for (int e=0;e<nr.elem_count;e++){
    for (int n=0;n<3;n++) elnode_h[3*e+n] = nr.elcon[3*e+n];
    // element.Push(new Element(nr.elcon[3*e],nr.elcon[3*e+1],nr.elcon[3*e+2]));
    double3 v = (node_h[nr.elcon[3*e]] + node_h[nr.elcon[3*e+1]] + node_h[nr.elcon[3*e+2]] ) / 3.0;    
    
    centroid_h[e] = v;
    double3 v1, v2;
    //In COUNTERCLOCKWISE
    v1 = node_h[nr.elcon[3*e+1]] - node_h[nr.elcon[3*e]];
    v2 = node_h[nr.elcon[3*e+2]] - node_h[nr.elcon[3*e]];
    normal_h[e] = cross(v1,v2);
    if (flipnormals)
      normal_h[e] = -normal_h[e];
    
    normal_h[e] = normal_h[e]/length(normal_h[e]);
    printf ("Normal %d %f %f %f \n",e, normal_h[e].x,normal_h[e].y,normal_h[e].z);
    // Vec3_t v;
		// if (dimension ==3) v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]] + *node[nr.elcon[3*e+2]] ) / 3. ;
    // else               v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]])  / 2. ;
    // element[e] -> centroid = v;
    // //TODO: CHANGE FOR CALCNORMALS
    // if (dimension==3){
      // Vec3_t v1, v2;
      // //In COUNTERCLOCKWISE
      // v1 = *node[nr.elcon[3*e+1]] - *node[nr.elcon[3*e]];
      // v2 = *node[nr.elcon[3*e+2]] - *node[nr.elcon[3*e]];
      // element[e] ->normal = cross (v1,v2);

      // element[e] ->normal /= Norm(element[e] ->normal);
      // //cout << "v1 "<< v1<< ", v2 " <<v2<< ", normal "<<element[e]->normal <<endl;
    // } else { //See calc normals
        // Vec3_t u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        // v[0] = -u[1];
        // v[1] =  u[0];
        // v[2] =  0.0;
        // element[e] -> normal = v/norm(v);
    // }  
 }
//  cout << "Generated "<<element.Size()<< " trimesh elements. "<<endl;  
//  
  m_v = m_w = make_double3(0.,0.,0.);
  
  cudaMemcpy(elnode, elnode_h, 3 * elemcount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(centroid, centroid_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);
  cudaMemcpy(normal, normal_h, elemcount * sizeof(double3), cudaMemcpyHostToDevice);

  cudaMalloc((void **)&pplane , 	elemcount * sizeof (double));
  cudaMalloc((void **)&nfar   , 	elemcount * sizeof (int));

  delete node_h;
  delete node_vh;
  delete elnode_h;
  delete centroid_h;
  delete normal_h; 
}

};