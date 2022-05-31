// TODO: extend to all dirs
//NOTE: DENSITY IS OF ELEMENTS
//This also will be passed to device
namespace SPH{
//NOW THIS IS ZORIENTED, CHANGE TO EVERY PLANE
inline void TriMesh_d::AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens){
	
	double x1,x2,x3;
	double l1,l2;
	double3 p = p2-p1;
	int dir[3];
	if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	else									{dir[0] = 0; dir[1] = 1;}
	
	dir [2] = axis; //dir2 is which remains constant
	
	x3 = p1.z;
	x2 = p1.y; 
  
  //TODO: CORRECT
  //x3 = p1.y;
  //x2 = p1.y;
  
	//double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
  double dl = p.x/dens;
  int nodecount = (dens+1)*(dens+1);
  // node = new double3 [nodecount];
  // node_v = new double3 [nodecount];
  
  //Is it necessary to paralellize mesh nodes??
  cudaMalloc((void **)&node   , 	nodecount * sizeof (double3));
  cudaMalloc((void **)&node_v , 	nodecount * sizeof (double3));
  
  double3 *node_h, *node_vh;
  node_h  =  new double3 [nodecount];
  node_vh =  new double3 [nodecount];
  
	//cout <<"dens: "<<dens<<endl;
	//Plane is in 0 and 1 dirs
  cout << "Creating nodes.."<<endl;
	int vi=0;
	int test =dens+1;
	for (int j=0; j<test; j++) {
		//x1 = p1(dir[0]);
    x1 = p1.x;
		for (int i=0; i<test; i++){
			double3 v;
			v.x=x1;v.y=x2;v.z=x3;
			//cout << "i,j" << i << ", " << j<<endl; 
			//node.Push(new double3(x1,x2,x3));
			node_h[vi]		=make_double3(v.x,v.y,v.z);
			node_vh[vi]	=make_double3(0.,0.,0.);
      vi++;
			// node.Push(new double3(v(0),v(1),v(2)));
			// node_v.Push(new double3(0.,0.,0.));
			//cout << "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			x1+=dl;
		}
		x2+=dl;
	}
  cudaMemcpy(node_h, node, nodecount, cudaMemcpyHostToDevice);
  cudaMemcpy(node_vh, node_v, nodecount, cudaMemcpyHostToDevice);

  cout << "Element count: "<<elcount << endl;  
  cout << "done. Creating elements... ";
	int n[4];
	int el =0;
	int i;
	
	int elcount = dens * dens * 2;
	cudaMalloc((void **)&centroid , 	elcount * sizeof (double3));
	cudaMalloc((void **)&normal 	, 	elcount * sizeof (double3));
	cudaMalloc((void **)&elnode 	, 	3 * elcount * sizeof (int));	
  int *elnode_h = new int[3*elcount];
  double3 *centroid_h = new double3[elcount];
  double3 *normal_h   = new double3[elcount];
	
	for (size_t j = 0 ;j  < dens; j++ ) {
				// cout <<"j, dens" <<j<<", "<<dens<<endl;
				// cout <<"j<dens"<< (j  < dens)<<endl;
		for ( i = 0; i < dens; i++ ){
				// cout <<"i, dens" <<i<<", "<<dens<<endl;
				// cout <<"i <dens"<< (i  < dens)<<endl;
				n[0] = (dens + 1)* j + i; 		n[1] = n[0] + 1; 
				n[2] = (dens + 1)* (j+1) + i; n[3] = n[2] + 1;
			//cout <<" jj" << jj<<endl;
			int elcon[2][3];	// TODO: check x, y and z normals and node direction 
												// For all plane orientations
			//If connectivity  is anticlockwise normal is outwards
			if (positaxisorent) {
				elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
				elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
			} else {
				elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
				elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
			}
			for ( int e= 0; e<2;e++) { // 2 triangles
				int elnodeid = 3*el;
				//element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
				elnode_h[elnodeid + 0] = elcon[e][0]; 
				elnode_h[elnodeid + 1] = elcon[e][1]; 
				elnode_h[elnodeid + 2] = elcon[e][2];
				//cout << "Element "<< el <<": ";
				// for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
				// cout <<endl;
				
				double3 v = ( node_h[elcon[e][0]] + node_h[elcon[e][1]] + node_h[elcon[e][2]] ) / 3. ;
				//element[el] -> centroid = v; 
				centroid_h[el] = v;
				//cout << "Centroid" << element[el] -> centroid << endl;
				el++;
			}
		}// i for
		
	}

	///////////////////////////////////////////
	//// MESH GENERATION END
	cout << endl<<"Done. Creating normals"<<endl;
	for (int e = 0; e < elcount; e++){ 
		double f=-1.;
		if (positaxisorent) f= 1.;
		//element[e] -> normal (axis) = f;
		if (axis == 0)			normal_h[e].x = f;
		else if (axis == 1)	normal_h[e].y = f;
		else 								normal_h[e].z = f;
	}

  cudaMemcpy(elnode_h, elnode, elcount, cudaMemcpyHostToDevice);
  cudaMemcpy(centroid_h, centroid, elcount, cudaMemcpyHostToDevice);
  cudaMemcpy(normal_h, normal, elcount, cudaMemcpyHostToDevice);

  delete node_h;
  delete elnode_h;
  delete centroid_h;
  delete normal_h;  
}

inline __device__ void TriMesh_d::UpdatePlaneCoeff(){
	//Update pplan
	//	pplane[i] = dot(*node [element[e] -> node[element[e] ->nfar]],element[e] -> normal);

}

};