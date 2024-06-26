#include <iostream>
#include "Mesh.h"
using namespace std;
namespace SPH {

	
TriMesh::TriMesh(){
	
	
}

void TriMesh::ReadFromNastran(NastranReader &nr, bool flipnormals){
  //dimension = nr.dim;
  //Insert nodes
  for (int n=0;n<nr.node_count;n++){
    if (!flipnormals)
      node.push_back(new Vector(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]));
    else 
      node.push_back(new Vector(nr.node[3*n+1],nr.node[3*n],nr.node[3*n+2]));
    
		node_v.push_back(new Vector(0.,0.,0.));
  }
  cout << "Generated "<<node.size()<< " trimesh nodes. "<<endl;
  //cout << "Normals"<<endl;
  cout << "Writing elements..."<<endl;
  for (int e=0;e<nr.elem_count;e++){
    element.push_back(new Element(nr.elcon[3*e],nr.elcon[3*e+1],nr.elcon[3*e+2]));		  
    Vector v;
		//if (dimension ==3) 
      v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]] + *node[nr.elcon[3*e+2]] ) / 3. ;
    //else               v = ( *node[nr.elcon[3*e]] + *node[nr.elcon[3*e+1]])  / 2. ;
    element[e] -> centroid = v;
    //TODO: CHANGE FOR CALCNORMALS
//    if (dimension==3){
      Vector v1, v2;
      //In COUNTERCLOCKWISE
      v1 = *node[nr.elcon[3*e+1]] - *node[nr.elcon[3*e]];
      v2 = *node[nr.elcon[3*e+2]] - *node[nr.elcon[3*e]];
      element[e] ->normal = cross (v1,v2);

      element[e] ->normal = element[e] ->normal/element[e] ->normal.norm();
      //cout << "v1 "<< v1<< ", v2 " <<v2<< ", normal "<<element[e]->normal <<endl;
    // } else { //See calc normals
        // Vec3_t u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
        // v[0] = -u[1];
        // v[1] =  u[0];
        // v[2] =  0.0;
        // element[e] -> normal = v/norm(v);
    // }
  }
  cout << "Generated "<<element.size()<< " trimesh elements. "<<endl;  
  
  // m_v = 0.;
  // m_w = 0.;  
}

Element::Element(const int &n1, const int &n2, const int &n3){
	
	//centroid = Vector();
	node[0] = n1; node [1] = n2; node[2] = n3;
	
}

void TriMesh::CalcCentroids(){
	
	for (int e=0;e<element.size();e++)
		element[e]-> centroid = ( *node[element[e]->node[0]] + *node[element[e]->node[1]] + *node[element[e]->node[2]] ) / 3.; 
	
}
// TODO: extend to all dirs
void TriMesh::AxisPlaneMesh(const int &axis, bool positaxisorent, const Vector p1, const Vector p2,  const int &dens){
	int elemcount = dens * dens;
	
	double x1,x2,x3;
	double l1,l2;
	Vector p = p2-p1;
	int dir[3];
	if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	else									{dir[0] = 0; dir[1] = 1;}
	
	dir [2] = axis; //dir2 is which remains constant
	
	x3 = p1(dir[2]);

	x2=p1(dir[1]); 
	double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
	//cout <<"dens: "<<dens<<endl;
	//Plane is in 0 and 1 dirs
	int test =dens+1;
	for (int j=0; j<test; j++) {
		x1 = p1(dir[0]);
		for (int i=0; i<test; i++){
			Vector v;
			v(dir[0])=x1;v(dir[1])=x2;v(dir[2])=x3;
			//cout << "i,j" << i << ", " << j<<endl; 
			//node.push_back(new Vector(x1,x2,x3));
			node.push_back(new Vector(v(0),v(1),v(2)));
			node_v.push_back(new Vector(0.,0.,0.));
			//cout << "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			x1+=dl;
		}
		x2+=dl;
	}

	int n[4];
	int el =0;
	int i;
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
			//cout << "elnodes"<<endl;
			for ( int e= 0; e<2;e++) { // 2 triangles
				element.push_back(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
				//cout << "Element "<< el <<": ";
				// for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
				// cout <<endl;
				
				Vector v = ( *node[elcon[e][0]] + *node[elcon[e][1]] + *node[elcon[e][2]] ) / 3. ;
				element[el] -> centroid = v; 
				//cout << "Centroid" << element[el] -> centroid << endl;
				el++;
			}
		}// i for
		
	}
	///////////////////////////////////////////
	//// MESH GENERATION END
	cout << "Creating normals"<<endl;
	for (int e = 0; e < element.size(); e++){ 
		double f=-1.;
		if (positaxisorent) f= 1.;
		element[e] -> normal (axis) = f;
	}

  cout << "Created Mesh with "<< node.size()<< " nodes. "<<endl;
  // if (node.size() == 0)
    // throw new Fatal("ATTENTION! Check mesh generation");
}

//This is done once, Since mesh is rigid
//Calculate radius and plane coefficient
 void TriMesh::CalcSpheres(){
	double max;
	for (int e = 0; e < element.size(); e++){ 
		max = 0.;
		Vector rv;
		for (int n = 0 ;n < 3; n++){
			rv = *node [element[e]->node[n]] - element[e] -> centroid;
			if (rv.norm() > max) max = rv.norm();
			element[e]-> nfar = n;
		}
		element[e]-> radius = max;	//Fraser Eq 3-136
	}
	UpdatePlaneCoeff();
	
}

 void TriMesh::UpdatePlaneCoeff(){
	//Update pplane
	for (int e = 0; e < element.size(); e++) 
		element[e]-> pplane = dot(*node [element[e] -> node[element[e] ->nfar]],element[e] -> normal);

}
 void TriMesh::CalcNormals(){
	Vector u, v, w;
	for (int e = 0; e < element.size(); e++) {
			u = *node [element[e]->node[1]] - *node [element[e]->node[0]];
			v = *node [element[e]->node[2]] - *node [element[e]->node[0]];
			w = cross(u,v);
			element[e] -> normal = w/w.norm();
			//Fraser Eqn 3.34
			//Uj x Vj / |UjxVj|
	}
}

 void TriMesh::ApplyConstVel(const Vector &v){
		for (int n=0;n<node.size();n++)
			*node_v[n] = v;
}

 void TriMesh::CalcCentroidVelFromNodes(){
	
	
}

 void TriMesh::UpdatePos(const double &dt){
	
	//Seems to be More accurate to do this by node vel
	//This is used by normals
	for (int n=0;n<node.size();n++){
		*node[n] = *node[n] + (*node_v[n])*dt;
	}
	
}

void TriMesh::Scale(const double &f){
	//Seems to be More accurate to do this by node vel
	//This is used by normals

	for (int n=0;n<node.size();n++){
		*node[n] *= f;
	} 
  cout << "calc centroids"<<endl;
  CalcCentroids();
  cout << "calc normals"<<endl;
  CalcNormals();        //From node positions
  cout << "generate plane coeffs"<<endl;
  //UpdatePlaneCoeff();   //pplane
}

void TriMesh::Move(const Vector &v){
	//Seems to be More accurate to do this by node vel
	//This is used by normals

	for (int n=0;n<node.size();n++){
		*node[n] += v;
	} 
  cout << "calc centroids"<<endl;
  CalcCentroids();
  cout << "calc normals"<<endl;
  CalcNormals();        //From node positions
  cout << "generate plane coeffs"<<endl;
  //UpdatePlaneCoeff();   //pplane
}
};
