/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Domain.h"
#include "Input.h"
#include <fstream>
#include <iostream>


#include "cuda/Domain_d.cuh" 
#include "cuda/Mechanical.cu" 
#include "cuda/SolverFraser.cu"

#include "cuda/Mesh.cuh"
//#include "cuda/Mesh.cu"
#include "cuda/Boundary_Condition.cuh"
#include "NastranReader.h" 

// #include "InteractionAlt.cpp"
// #include "Mesh.h"

// #include "SolverFraser.cpp"
// #include "Geometry.cpp"
// #include "SolverKickDrift.cpp"
void report_gpu_mem()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);
    std::cout << "Free = " << free << " Total = " << total <<std::endl;
}


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
//https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define TAU		0.005
#define VMAX	10.0

#define PRINTVEC(v)	cout << v[0]<<", "<<v[1]<<", "<<v[2]<<endl;

using namespace std;
using namespace SPH;

// template <typename T > 
// void ReadParameter(T in, nlohmann::json in) {
  // cout << "Reading Configuration parameters..."<<endl; 
  // cout << "Time step size: ";
  // readValue(in["timeStepSize"], /*scene.timeStepSize*/ts);
  // cout << ts << endl;
// }


void UserAcc(SPH::Domain_d & domi) {

	// for (int i=0; i<domi.particle_count; i++) { 
		// // for (int bc=0;bc<domi.bConds.size();bc++){
			// // if (domi.Particles[i]->ID == domi.bConds[bc].zoneId ) {
				// // if (domi.bConds[bc].type == 0 ){ //VELOCITY
					// // if (domi.bConds[bc].valueType == 0) {
            // // domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
            // // domi.Particles[i]->v		= domi.bConds[bc].value;          
          // // } else if (domi.bConds[bc].valueType == 1) {///amplitude
            // // for (int j=0;j<domi.amps.size();j++){
              // // if(domi.amps[j].id == domi.bConds[bc].ampId){
                // // double val = domi.bConds[bc].ampFactor * domi.amps[j].getValAtTime(domi.getTime());
                // // Vec3_t vec = val * domi.bConds[bc].value;
                // // domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
                // // domi.Particles[i]->v		= vec;
              // // }//if if match
            // // }//for amps
          // // } //VALUE TYPE == AMPLITUDE 
				// // }//TYPE == VELOCITY
        
			// // }//ZoneID 			
		// // }//BC


    for (int bc=0;bc<domi.bConds.size();bc++) {
      ApplyBCVelKernel	<<<domi.blocksPerGrid,domi.threadsPerBlock >>>(&domi, domi.bConds[bc].zoneId, domi.bConds[bc].value);
      cudaDeviceSynchronize();
    }

   //cout << "Applying BC"<<endl;
		// ApplyBCVelKernel	<<<domi.blocksPerGrid,domi.threadsPerBlock >>>(&domi, 1, make_double3(0.,0.,0.));
		// cudaDeviceSynchronize();
    // double vbc;
    // if (domi.Time < TAU) vbc = VMAX/TAU*domi.Time;
    // else            vbc = VMAX;
    // //cout << "vbc "<<vbc<<endl;
		// ApplyBCVelKernel	<<<domi.blocksPerGrid,domi.threadsPerBlock >>>(&domi, 2, make_double3(0.,0.,-vbc));
		// cudaDeviceSynchronize();
    
  if (domi.contact){
    for (int bc=0;bc<domi.bConds.size();bc++){
      for (int m=0;m<domi.trimesh_count;m++){
        int id;
        // THIS IS NOT WORKING
        // getTrimeshIDKernel<<<1,1>>>(&domi,m,&id);
        // cudaDeviceSynchronize();
        //cout << "mesh id "<<id<<endl;
        // if (domi.trimesh[m]->id == domi.bConds[bc].zoneId)
        // //if ( (getTrimeshIDKernel<<<1,1>>>(&domi,m)) == domi.bConds[bc].zoneId)
          // if (domi.bConds[bc].valueType == 0) { ///constant
            // //OLD, when trimesh was not a vector
            // //domi.trimesh[m]->SetVel(domi.bConds[bc].value);
            
            //SetMeshVelKernel<<<1,1>>>(&domi,m, domi.bConds[bc].value);
            // //domi.trimesh->SetRotAxisVel(domi.bConds[bc].value_ang);
          // }//BCOND 

      //}//mesh
    }//bcs
    }//trimesh
  }//contact
}

int main(int argc, char **argv)
{

	if (argc > 1) {
		string inputFileName=argv[1];	
		std::ifstream ifs(argv[1]);
		json j;
		ifs >> j;
		
		nlohmann::json config 		= j["Configuration"];
		nlohmann::json material 	= j["Materials"];
		nlohmann::json domblock 	= j["DomainBlocks"];
		nlohmann::json domzones 	= j["DomainZones"];
		nlohmann::json amplitudes 	= j["Amplitudes"];
		nlohmann::json rigbodies 		= j["RigidBodies"];
    nlohmann::json contact_ 		= j["Contact"];
		nlohmann::json bcs 			= j["BoundaryConditions"];
		nlohmann::json ics 			= j["InitialConditions"];

		
		SPH::Domain	dom; //TODO: DELETE THIS AND PASS TO DOMAIN
    SPH::Domain_d *dom_d;
    std::vector<TriMesh *> mesh; ////// TODO: ALLOW FOR MULTIPLE MESH CONTACT
    //std::vector<SPH::TriMesh_d *> mesh_d;
    
    SPH::TriMesh_d *mesh_d;
    
    SPH::TriMesh_d *unifmesh_d;
    
    report_gpu_mem();
    gpuErrchk(cudaMallocManaged(&dom_d, sizeof(SPH::Domain_d)) );
    report_gpu_mem();
    dom_d->dom_bid_type = None;
		dom.Dimension	= 3;
    
    dom_d->trimesh_count = 0;
		
		// string kernel;
    double ts;
    
    cout << "--------------------------------------------"<<endl;
    cout << "----------------- WELDFORM GPU-----------------"<<endl;
    cout << "----------------- v. 0.0.1 -----------------"<<endl;
    cout << "--------------------------------------------"<<endl<<endl<<endl;
    
    // cout << "Reading Configuration parameters..."<<endl; 
    
    // int np = 4;
    // readValue(config["Nproc"], np);   
    // dom.Nproc	= np;
        
    // string sumType = "Nishimura";
    // readValue(config["sumType"], /*scene.timeStepSize*/sumType);
    // if (sumType == "Locking") dom.nonlock_sum = false;
    // else if (sumType == "Nishimura") dom.nonlock_sum = true;
    // else cout << "sumType value not valid. Options are \"Locking\" and \"Nishimura\". "<<endl; 
    
    
		cout << "Time step size: ";
    readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
    // cout << ts << endl;
    // dom.Kernel_Set(Qubic_Spline);
    
    
    // string solver = "Mech";
    // readValue(config["solver"],solver);
    
    // if (solver=="Mech-Thermal")
      // dom.thermal_solver = true;
		
		// readValue(config["integrationMethod"], dom.Scheme); //0:Verlet, 1:LeapFrog, 2: Modified Verlet

     	// //dom.XSPH	= 0.5; //Very important

    double dx,r,h;
    bool h_update = false;
    double hfactor = 0.0;
    r = dx = 0.0; //If not set, h is calculated from mass and density
    
    cout << "Reading Config section ... "<<endl;
		readValue(config["particleRadius"], r);
		readValue(config["hFactor"], hfactor);



		
	
		// //////////////
		// // MATERIAL //
		// //////////////
		double rho,E,nu,K,G,Cs,Fy;
    double Et, Ep;  //Hardening (only for bilinear and multilear)
    std::vector<double> c;
    // c.resize(10);
    string mattype = "Bilinear";
    cout << "Reading Material.."<<endl;
    cout << "Type.."<< endl; readValue(material[0]["type"], 		mattype);
    readValue(material[0]["density0"], 		rho);
    readValue(material[0]["youngsModulus"], 	E);
    readValue(material[0]["poissonsRatio"], 	nu);
    readValue(material[0]["yieldStress0"], 	Fy);
    readArray(material[0]["const"], 		c);

    
    
    // // THERMAL PROPERTIES

    // double k_T, cp_T;
    // readValue(material[0]["thermalCond"], 	  k_T);
    // readValue(material[0]["thermalHeatCap"], 	cp_T);    
    
    // cout << "Done. "<<endl;
       
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));

		dx 	= 2.*r;
    bool h_fixed = false;
    h  = dx*hfactor; //Very important
    if (h>0.0) {
      cout << "Smoothing Length forced to: " << h <<endl;
      h_fixed = true;
    } else {
      cout << "Smoothing length calculated from Mass and density"<<endl;
    }

    h	= dx*hfactor; //Very important
    Cs	= sqrt(K/rho);
    
    cout << "Cs: "<<sqrt(K/rho);
        // double timestep,cflFactor;
		// int cflMethod;
    double sim_time, output_time;
    // string cont_alg = "Fraser";
    bool auto_ts[] = {true, false, false}; //ONLY VEL CRITERIA


		// readValue(config["cflMethod"], cflMethod);
		// if (cflMethod == 0)
			// readValue(config["timeStepSize"], timestep);
		// else {
			// readValue(config["cflFactor"], cflFactor);
			// timestep = (cflFactor*h/(Cs));
		// }
    readValue(config["outTime"], output_time);
    readValue(config["simTime"], sim_time);
    readBoolVector(config["autoTS"], auto_ts);
    // double alpha = 1.;
    // double beta = 0.;
    // bool h_upd = false;
    // double tensins = 0.3;
    // bool kernel_grad_corr = false;
    readValue(config["artifViscAlpha"],dom_d->Alpha); //TODO: ARTIFF VISC PER PARTICLE
    readValue(config["artifViscBeta"],dom_d->Beta);
    // readValue(config["contAlgorithm"],cont_alg);
    // readValue(config["kernelGradCorr"],kernel_grad_corr);
    // readValue(config["smoothlenUpdate"],h_upd);
    dom_d->auto_ts = auto_ts[0];
    // dom.auto_ts_acc = auto_ts[1];
    // dom.auto_ts_cont = auto_ts[2];
    // readValue(config["tensileInstability"],tensins);
    // if (h_upd) //default is false...
      // dom.h_update = true;
    // /////////////-/////////////////////////////////////////////////////////////////////////////////
		// // DOMAIN //
		// ////////////
		double3 start,L;
    int id;
		string domtype = "Box";
    int matID;
    string gridCS = "Cartesian";
    double slice_ang = 2.0000001 * M_PI;
    bool sym[] = {false,false,false};
		readValue(domblock[0]["id"], 	id);
		readVector(domblock[0]["start"], 	start);
		cout << "Reading Domain dim" << endl;  readVector(domblock[0]["dim"], 	L);
		cout << "Reading Domain type" << endl; readValue(domblock[0]["type"], 	domtype); //0: Box
    cout << "Reading Domain mat id" << endl;  readValue(domblock[0]["matID"], 	matID); //0: Box
    cout << "Grid Coordinate System" << endl;  readValue(domblock[0]["gridCoordSys"], 	gridCS); //0: Box
    cout << "Slice Angle " << endl;  readValue(domblock[0]["sliceAngle"], 	slice_ang); //0: Box
    readBoolVector(domblock[0]["sym"], 	sym); //0: Box
     for (int i=0;i<3;i++) {//TODO: Increment by Start Vector
			// dom.DomMax(0) = L[i];
			// dom.DomMin(0) = -L[i];
		}		
    

		
		// // inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// // double r, double Density, double h, bool Fixed) {
												
		// //dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10. + dx, r, rho, h, false); 
		
    // if (abs(L[2]) < h ) {
      // dom.Dimension = 2;
      // cout << "Z Value is less than h. Dimension is set to 2. "<<endl;
      // cout << "Dimension also could be set in config section." <<endl;
    // }
    
		cout << "Dimensions: "<<endl;
		//PRINTVEC(L)
		if (domtype == "Box"){
      cout << "Adding Box ..."<<endl;      
			dom.AddBoxLength(0 ,start, L.x , L.y,  L.z , r ,rho, h, 1 , 0 , false, false );		
      dom_d->particle_count = dom.Particles.size(); ///// IN THE FUTURE DDOMAIN_D WILL MAKE 
		}   else if (domtype == "Cylinder"){
      cout << "Adding Cylinder";      
			// if (sym[0] && sym[1]){
        // cout << " with symmetry..."<<endl;
        // dom.AddXYSymCylinderLength(0, L[0]/2., L[2], r, rho, h, false, sym[2]); 
      // }
      // else {
        // cout << "..."<<endl;
        if ( gridCS == "Cartesian"){
          cout << "Reserved "<<ComputeCylinderParticles (L.x/2., L.z, r)<<" particles."<<endl;
        cout << "Length " << L.x<<", "<< L.y<<", "<< L.z<<", "<<endl;
          dom.AddCylinderLength(0, start, L.x/2., L.z, r, rho, h, false);  /////// GENERATED AT HOST TO THEN COPY
//void Domain::AddCylUniformLength(int tag, double Rxy, double Lz, 
//																				double r, double Density, double h) 
          // dom.AddCylUniformLength(0, L.x/2.0, L.z, 
																				// r, rho, h, M_PI/16.0, 1, L.x/4.0); 

          // dom.AddCylUniformLength(0, L.x/2.0, L.z, 
																				// r, rho, h); 


        }else if (gridCS == "Cylindrical"){
          // dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h);
          if (slice_ang==2.0 * M_PI){
          dom.AddCylUniformLength(0, L.x/2.0, L.z, 
																				r, rho, h); 
          } else {
            dom.AddCylUniformLength(0, L.x/2.0, L.z, 
																				r, rho, h, M_PI/16.0, 1, L.x/4.0); 
            dom_d->dom_bid_type = AxiSymm_3D;
          }          
       }//Cylindrical
        dom_d->particle_count = dom.Particles.size(); ///// IN THE FUTURE DDOMAIN_D WILL MAKE 
        
        
    } else if (domtype == "File"){ //DECIDE ACCORDING TO EXTENSION
        string filename = "";
        readValue(domblock[0]["fileName"], 	filename); 
        cout << "Reading Particles Input file " << filename <<endl;  
        dom_d->ReadFromLSdyna(filename.c_str());
        
        double tot_mass = 0.;
        for (int p=0;p<dom_d->particle_count;p++){
          double x,y,z;
          x =dom_d->x_h[p].x;
          y =dom_d->x_h[p].y;
          z = dom_d->x_h[p].z;
          dom.Particles.push_back(new SPH::Particle(0,Vector(x,y,z),Vector(0,0,0),0.0,rho,h,false));
          dom.Particles[p]->Mass = dom_d->m_h[p];
          //if (!h_fixed) dom.Particles[p]->h = pow(dom.Particles[p]->Mass/rho,0.33333) * 1.2; ////// THIS CRASHES
          if (dom_d->realloc_ID)dom.Particles[p]->ID = dom_d->ID_h[p];
          tot_mass+=dom_d->m_h[p];
        }
        delete dom_d->x_h,dom_d->m_h;
        printf( "Total Mass Readed from LS-Dyna: %fn", tot_mass);
    }
	

		cout << "Particle count: "<<dom.Particles.size()<<endl;

    cout << "Domain Zones "<<domzones.size()<<endl;		
		for (auto& zone : domzones) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid;
      Vector vstart, vend;
			readValue(zone["id"], 		zoneid);
			readVector(zone["start"], 	vstart);
			readVector(zone["end"], 	vend);
      cout << "Zone id "<<zoneid<<endl;
			cout << "Dimensions: "<<endl;
			cout << "start"<< vstart(0)<<"; "<< vstart(1)<<"; "<< vstart(2)<<"; "<<endl;
			cout << "end"<< vend(0)<<"; "<< vend(1)<<"; "<< vend(2)<<"; "<<endl;
      
			int partcount =dom.AssignZone(vstart,vend,zoneid); ////IN DEVICE DOMAINf
      std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
      if (partcount == 0) cout << "-----WARNING------: Zone " <<zoneid<<" has ZERO particles"<<endl;
		}
    
    // //////////////////////////////////////////////////////////
    // ////////////////// RIGID BODIES //////////////////////////
    string rigbody_type;
    bool contact = false;
    if (readValue(rigbodies[0]["type"],rigbody_type))
      contact = true;
    double3 dim;
    
	readVector(rigbodies[0]["start"], 	start);       
	readVector(rigbodies[0]["dim"], 	dim); 
    bool flipnormals = false;
    readValue(rigbodies[0]["flipNormals"],flipnormals);
    
    // double heatcap = 1.;
    // readValue(rigbodies[0]["thermalHeatCap"],heatcap);
    // //TODO: WRitE TO PArTiclES
    // if (rigbody_type == "File"){
      // // string filename = "";
      // // readValue(rigbodies[0]["fileName"], 	filename); 
      // // cout << "Reading Mesh input file..." << endl;
      // // SPH::NastranReader reader("Tool.nas", flipnormals);
    // }
    // else {
      // if (dim (0)!=0. && dim(1) != 0. && dim(2) !=0. && rigbody_type == "Plane")
        // throw new Fatal("ERROR: Contact Plane Surface should have one null dimension");
    // }


    //BEFORE CONTACT

    dom_d->solid_part_count = dom_d->particle_count;  //AFTER SET DIMENSION
    
    int bc_count = 0;
    std::vector<boundaryCondition> bcondvec;
		for (auto& bc : bcs) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype,var,ampid;

			double ampfactor;
			bool free=true;
			SPH::boundaryCondition bcon;
      bcon.type = 0;        //DEFAULT: VELOCITY
      bcon.valueType = 0;   //DEFAULT: CONSTANT
      bcon.value_ang = make_double3 (0.0);
			readValue(bc["zoneId"], 	bcon.zoneId);
      //type 0 means velocity vc
			readValue(bc["valueType"], 	bcon.valueType);
			if (bcon.valueType == 0){//Constant
        readVector(bc["value"], 	      bcon.value);      //Or value linear
        readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value
      } else 
        if ( bcon.valueType == 1){ //Amplitude
				readValue(bc["amplitudeId"], 		bcon.ampId);
				readValue(bc["amplitudeFactor"], 	bcon.ampFactor);
			}
				
			readValue(bc["free"], 	bcon.free);
			dom_d->bConds.push_back(bcon);
      bcondvec.push_back(bcon);
      bc_count++;
			
      std::cout<< "BCs "<<  ", Zone ID: "<<bcon.zoneId<<", Value :" <<bcon.value.x<<", "<<bcon.value.y<<", "<<bcon.value.z<<std::endl;
		}//Boundary Conditions
    //dom_d->bConds

    std::vector<int> mesh_id(dom_d->particle_count);
    for (int p=0;p<dom.Particles.size();p++) mesh_id[p]=0; //THIS HAS TO BE CHANGED, ONLKY TO RIGID PARTICLES
    
    cout << "Set contact to ";
    int solid_part_count = dom_d->solid_part_count;
    int last_pcount = dom_d->solid_part_count;;
    
    
    /////////////////////////////////////////// MESH READING /////////////////////////////////////
    ////////THIS IS FOR FLATTEN MESH
    int tot_node_count = 0;
    int tot_elem_count = 0;
    
    if (contact){
      cout << "true."<<endl;
      cout << "Rigid body Count: " << rigbodies.size() << endl;
  		dom_d->contact = true; //ATTENTION: SetDimension sets contact to OFF so...
      cout << "Reading contact mesh..."<<endl;
      
      ///// ORGINAL
      // SPH::TriMesh_d *mesh_d;
      //gpuErrchk(cudaMallocManaged(&mesh_d, sizeof(SPH::TriMesh_d)) );
      //cudaMalloc((void**)&dom_d->trimesh,         rigbodies.size()* sizeof(SPH::TriMesh_d));
      //mesh_d.resize(rigbodies.size());
      
      
      //For m
      gpuErrchk(cudaMallocManaged(&mesh_d, rigbodies.size()*sizeof(SPH::TriMesh_d)) );
      cudaMalloc((void**)&dom_d->trimesh,         rigbodies.size()* sizeof(SPH::TriMesh_d*));
      cudaMalloc((void**)&dom_d->contact_surf_id, rigbodies.size()* sizeof(int));     
      cudaMalloc((void**)&dom_d->first_fem_particle_idx, rigbodies.size()* sizeof(int));   

      //dom_d->trimesh = new SPH::TriMesh_d* [rigbodies.size()];
      
      ////// first_fem_particle_idx BEFORE CREATING PARTICLES
      //dom_d->first_fem_particle_idx = new int[rigbodies.size()]; // TODO: THIS SHOULD BE DONE AUTOMATICALLY
      dom_d->first_fem_particle_idx_0 = dom_d->particle_count;
      //cout << "Max Solid Part Count "<< dom_d->first_fem_particle_idx_0;
       
      for (int m=0;m<rigbodies.size();m++){
       
      if (readValue(rigbodies[m]["type"],rigbody_type))
      double3 dim;
    
      readVector(rigbodies[m]["start"], 	start);       
      readVector(rigbodies[m]["dim"], 	dim); 
      bool flipnormals = false;
      readValue(rigbodies[m]["flipNormals"],flipnormals);

      //gpuErrchk(cudaMallocManaged(&mesh_d[m], sizeof(SPH::TriMesh_d)) );   

        
        //TODO: CONVERT TO ARRAY std::vector<SPH::TriMesh_d> *mesh_d;
        //TODO: CHANGE TO EVERY DIRECTION
        int dens = 10;
        readValue(rigbodies[m]["partSide"],dens);
        if (rigbody_type == "Plane"){
          // TODO: CHECK IF MESH IS NOT DEFINED
          //mesh_d.push_back(new TriMesh);
          cout << "Mesh Dimensions: "<< dim.x <<", "<<dim.y<< ", "<<dim.z<<endl;
          mesh.push_back(new TriMesh);
          mesh[m]->AxisPlaneMesh(2, false, start, Vector(start.x + dim.x,start.y + dim.y , start.z),dens);
          cout << "Creating Device mesh"<< endl;
          mesh_d[m].AxisPlaneMesh(2, false, start, make_double3(start.x + dim.x,start.y + dim.y , start.z),dens);
          cout << "Done creating Device and host meshes" <<endl;
        } else if (rigbody_type == "File"){

          Vector md = 0.0;
          string filename = "";
          readValue(rigbodies[m]["fileName"], 	filename); 
          readVector(rigbodies[m]["moveDir"],md);       
          //readValue(rigbodies[0]["scaleFactor"],scalefactor);           
          cout << "Reading Mesh input file " << filename <<endl;
          NastranReader reader((char*) filename.c_str());
          mesh_d[m].ReadFromNastran(reader,false);
          //mesh_d->Move(make_double3(md[0],md[1],md[2]));
          //mesh_d = New
          //mesh.push_back (new SPH::TriMesh(reader,flipnormals ));
          mesh.push_back (new SPH::TriMesh);
          mesh[m]->ReadFromNastran(reader, false);
          
          mesh[m]->Move(md); //TODO: DELETE
        }
        tot_node_count += mesh_d[m].nodecount;
        tot_elem_count += mesh_d[m].elemcount;
                
        double hfac = 1.1;	//Used only for Neighbour search radius cutoff

        
        readValue(rigbodies[m]["zoneId"],id);
        dom.AddTrimeshParticles(*mesh[m], hfac, id); 
        //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
        
        //mesh_d[m].id = id; //ERROR; IS POSITION OR ID??? Check Also AssignTrimeshID Kernel
        //cout << "mesh id "<<mesh_d[m].id<<endl;
        
        
        AddTrimeshParticlesKernel<<<1,1>>>(dom_d, &mesh_d[m], hfac, id);
        gpuErrchk( cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        //BEFORE ALLOCATING 
        cout << "Allocating ..."<<endl;
        cout << "Assigning "<<endl;

        AssignTrimeshAddressKernel<<<1,1 >>>(dom_d,m,&mesh_d[m]);
        gpuErrchk( cudaPeekAtLastError() );
        cudaDeviceSynchronize();
       

        cout << "Ok."<<endl;
        //SetMeshIDKernel<<<1,1>>>(dom_d,m, id);
        //cudaDeviceSynchronize();

        int id_int;
        // getTrimeshIDKernel<<<1,1>>>(dom_d,m,&id_int);
        // cout << "MESH ID: "<<id_int<<endl;
        //cudaDeviceSynchronize();
        //CRASHES
        //cudaMemcpy(&id_int, &dom_d->trimesh[m]->id, sizeof (int), cudaMemcpyDeviceToHost);
        //ADDING MESH ID
        cout << "added particles count"<<dom.Particles.size() - last_pcount<<"with mesh_id ="<< m <<endl;
        for (int p=last_pcount;p<dom.Particles.size();p++) mesh_id.push_back(m);
        last_pcount = dom.Particles.size();

        for (int bc=0;bc<dom_d->bConds.size();bc++){
            if (id == dom_d->bConds[bc].zoneId){
              if (dom_d->bConds[bc].valueType == 0) { ///constant

            //NEW FORMAT, ASSIGNING node_v directly
            SetNodesVelKernel<<<1,1>>>(&mesh_d[m],dom_d->bConds[bc].value);
            cudaDeviceSynchronize();   
            // //OLD, when trimesh was not a vector
            // //domi.trimesh[m]->SetVel(domi.bConds[bc].value);
            
            
            /*
              SetMeshVelKernel<<<1,1>>>(dom_d,m, dom_d->bConds[bc].value);
              cudaDeviceSynchronize();  
                cout << "Mesh ID " << id << ", "<< "Velocity set to : "
                                            <<dom_d->bConds[bc].value.x 
                                            <<", " <<dom_d->bConds[bc].value.y 
                                            <<", " << dom_d->bConds[bc].value.z<<endl;
              */
              //AFTER SET VEL WHICH SETS m_V
              
              
              

              }
            }
        }
        
        cout << "Surf id int: "<<id_int<<endl;
        //dom_d->trimesh[0] = mesh_d; //TODO: CHECK WHY ADDRESS IS LOST
        cout << "Assigned "<<endl;
        
      }//MESH m
      
      /////////////////////////////////////////////////////////////
      ////////////////////// FLATTEN FOR SEVERAL SRFACES /////////
      //////////////////////// SEVERAL DEVICE POINTERS WTH ERRORS/////
      
      
      cout << "Node count: %d, Elem count \n"<< tot_node_count<<", "<<tot_elem_count<<endl;
      
      //NOW CREATES THE MESH
      //gpuErrchk(cudaMallocManaged(&unifmesh_d, rigbodies.size()*sizeof(SPH::TriMesh_d)) );     
      cudaMalloc((void**)&unifmesh_d,         sizeof(SPH::TriMesh_d));
      //unifmesh_d->setDim(tot_node_count,tot_elem_count);
      setDimKernel<<<1,1>>>(unifmesh_d,tot_node_count,tot_elem_count);
      
      //cudaMalloc((void**)&unifmesh_d->node, tot_node_count*sizeof(double3));


      /*
      for (int m=0;m<rigbodies.size();m++){
        //unifmesh_d = ;
        addMeshKernel<<<1,1>>>(unifmesh_d,&mesh_d[m]);
        cudaDeviceSynchronize();
      }
      */
      
      //CHANGING ADDRESS
      //AssignTrimeshAddressKernel<<<1,1 >>>(dom_d,0,unifmesh_d);
      //gpuErrchk( cudaPeekAtLastError() );
      //cudaDeviceSynchronize();
      
      
      //// DELETE PREVIOUS MESHES
      
      
      cout << "Assigning contact params"<<endl;
				
      double penaltyfac = 0.5;
      std::vector<double> fric_sta(1), fric_dyn(1), heat_cond(1);
      readValue(contact_[0]["fricCoeffStatic"], 	fric_sta[0]); 
      readValue(contact_[0]["fricCoeffDynamic"], 	fric_dyn[0]); 
      readValue(contact_[0]["heatCondCoeff"], 	  heat_cond[0]);
      
      // bool heat_cond_ = false;
      // if (readValue(contact_[0]["heatConductance"], 	heat_cond_)){
        // dom.cont_heat_cond = true;
        // dom.contact_hc = heat_cond[0];
      // }
      
      dom_d->friction_dyn = fric_dyn[0];
      dom_d->friction_sta = fric_sta[0];
      cout << "Contact Static Friction Coefficient: "<<dom_d->friction_dyn<<endl;
      // dom.PFAC = 0.8;
      // dom.DFAC = 0.0;
      
		} 
    else 
      cout << "false. "<<endl;

		// std::vector <SPH::amplitude> amps;
		
		// for (auto& ampl : amplitudes) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// // MaterialData* data = new MaterialData();
			// int zoneid,valuetype;
			// std::vector<double> time, value;
			// readValue(ampl["id"], 		zoneid);
			// //readValue(zone["valueType"],zoneid);
			// readArray(ampl["time"], 	time);
			// readValue(ampl["value"], 	value);
			// SPH::amplitude amp;
			// for (int i=0;i<time.size();i++){
				// amp.time.push_back(time[i]);
				// amp.value.push_back(value[i]);
			// }
			// amps.push_back(amp);
			// //std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
		// }
    
    // boundaryCondition *bConds_h    =  new boundaryCondition [bc_count];
    // for (int b=0;b<bc_count;b++)
      // bConds_h[b] = bcondvec[b];
    // cudaMalloc((void**)&dom_d->bConds, bc_count * sizeof(boundaryCondition ));
    // cudaMemcpy(dom_d->bConds, bConds_h, bc_count * sizeof(boundaryCondition), cudaMemcpyHostToDevice);
		
		// double IniTemp = 0.;
		// for (auto& ic : ics){
			// double temp;
			// if (solver == "Mech-Thermal"){
				// readValue(ic["Temp"], IniTemp);
				// cout << "Initial Temp: "<<IniTemp<<endl;
			// }
		// }
    
    // //Add fixed particles, these have priority
    
    // //TODO: CHECK IF DIFFERENT ZONES ARE INTERF
    // //Generate Domain
    // dom.gradKernelCorr = kernel_grad_corr;
    // dom.ts_nb_inc = 5;
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// TODO: THIS SHOULD BE PASSED TO A DOM_D INTIIALIZE FUNCTION ///////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  dom_d->GeneralAfter = & UserAcc;
	dom_d->SetDimension(dom.Particles.size());	 //AFTER CREATING DOMAIN

		 
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x    =  new double3 [dom.Particles.size()];
	int     *mid  =  new int [dom.Particles.size()];  
  
	for (int i=0;i<dom.Particles.size();i++){
		x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
    mid[i] = mesh_id[i];
  }
	int size = dom.Particles.size() * sizeof(double3);
  int size_int = dom.Particles.size() * sizeof(int);
	cout << "Copying to device "<<dom.Particles.size() << " particle properties ..."<<endl;
	cudaMemcpy(dom_d->x, x, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dom_d->mesh_id, mid, size, cudaMemcpyHostToDevice);


	for (int i=0;i<dom.Particles.size();i++){
		x[i] = make_double3(0.,0.,0.);
	}
	cudaMemcpy(dom_d->v, x, size, cudaMemcpyHostToDevice);
  
	cout << "copied"<<endl;

  // ////// MATERIAL  
  // Material_ *mat_h = (Material_ *)malloc(dom_d->solid_part_count * sizeof(Material_ *)); 
  

  
  // Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;
  // if      (mattype == "Bilinear")    {
    // Material_ *material_h  = new Bilinear();
    // Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    // cout << "Material Constants, Et: "<<c[0]<<endl;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear ));
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	
    // cout << "Created Bilinear material"<<endl;
  // } 
  // else if (mattype == "Hollomon")    {
    // Material_ *material_h  = new Hollomon(el,Fy,c[0],c[1]);
    // cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
  // } else if (mattype == "JohnsonCook") {
    // //Order is 
                               // //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
    // Material_ *material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
    
    // //Only 1 material to begin with
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    // cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
  // } else                              printf("ERROR: Invalid material type.");
    
////// MATERIAL  
  //Material_ *material_h = (Material_ *)malloc(dom_d->solid_part_count * sizeof(Material_ *)); 
  

  
  Elastic_ el(E,nu);
  cout << "Mat type  "<<mattype<<endl;
  cout << "E: "<<E<<endl;
  Material_ *material_h  = new Material_(el);
  if      (mattype == "Bilinear")    {

    Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    material_h->SetEp(Ep);
    cout << "Material Constants, Et: "<<c[0]<<endl;
    material_h->Et = c[0];

    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear ));
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	
    cout << "Created Bilinear material"<<endl;


  } 
  else if (mattype == "Hollomon")    {
//    Material_ *material_h  = new Hollomon(el,Fy,c[0],c[1]);
//    cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
  } else if (mattype == "JohnsonCook") {
    //Order is 
      // void Init_JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              // const double &c, const double &eps_0_,
              // const double &m_, const double &T_m_, const double &T_t_)
    //material_h->Init_JohnsonCook(el,c[0],c[1],c[2]);
    material_h->A= c[0];
    material_h->B= c[1];
    material_h->n= c[2];
    material_h->C= c[3];
    material_h->eps_0 = c[4];
    material_h->m= c[5];
    material_h->T_m= c[6];
    material_h->T_t= c[7];    
    
    // material_h->eps_0= c[6];
                               //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
//    Material_ *material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
//    
//    //Only 1 material to begin with
//    cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
//    cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    cout << "JOHNSON COOK Material Constants, A: "<< material_h->A << ", B: "<<c[1]<< ", n: "<<c[2]<<", C: "<<c[3]<< ", eps_0: "<<c[4]<<", m: "<<c[5]<<", T_m: "<<c[6]<<", T_t: "<<c[6]<<endl;
    Fy = CalcJohnsonCookYieldStress(0.0,0.0,0.0,material_h);
  } else                              
  	printf("ERROR: Invalid material type.");
  
  cout << "Elastic Properties "<<endl;
  cout << "Young Modulus: "<<material_h->elastic_m.E()<<endl;
  
    cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_ ));
    cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);
	
	cout << "Setting values"<<endl;
	dom_d->SetDensity(rho);
	dom_d->Set_h(h);
	cout << "done."<<endl;
  
  bool mass_ok = true;
  double totmass = 0.;
	double *m =  new double [dom.Particles.size()];
  //double *h_ = new double [dom.Particles.size()];
	for (size_t a=0; a<dom.Particles.size(); a++){
		m[a] = dom.Particles[a]->Mass;
		totmass +=m[a];
    if (m[a] < 1.0e-10 && a< dom_d->solid_part_count) {
      mass_ok = false;
    }
    //if (!h_fixed) h_[a] = dom.Particles[a]->h;
	}
  //delete m, h_;
  if (!mass_ok) cout << "-----WARNING ---- some particles have ZERO MASS"<<endl;
  cout << "Total Mass: "<<totmass<<endl;
  cudaMemcpy(dom_d->m, m, dom.Particles.size() * sizeof(double), cudaMemcpyHostToDevice);	
  cudaMemcpy(&dom_d->totmass, &totmass, sizeof(double),cudaMemcpyHostToDevice);
  // if (!h_fixed)
    // cudaMemcpy(dom_d->h, h_, dom.Particles.size() * sizeof(double), cudaMemcpyHostToDevice);
  
  dom_d->SetShearModulus(G);	// 
  for (size_t a=0; a<dom.Particles.size(); a++) {
    //dom.Particles[a]->G				= G; 
    dom.Particles[a]->PresEq	= 0;
    dom.Particles[a]->Cs			= Cs;

    dom.Particles[a]->TI		= 0.3;
    dom.Particles[a]->TIInitDist	= dx;
  }// particles
  
  dom_d->SetFreePart(dom); //All set to IsFree = true in this example
  dom_d->SetID(dom); 
  dom_d->SetCs(dom);
  
  dom_d->SetSigmay(Fy);
  cout << "Initial Yield Stress set as :" << Fy<<endl;

  ///////////////////////////////// IF CONTACT 
  ////////////////////////////////////////////
  int nwcount = 0;
  bool *not_write = new bool[solid_part_count];
  for (int i=0;i< solid_part_count;i++){
    not_write[i] = false;
    if (dom.Particles[i]->ID != 0){
      not_write[i] = true;
      //cout << "ID "<<dom.Particles[i]->ID <<endl;
      nwcount++;
    }
  }
  cout << "Set "<<nwcount<<" particles fixed ID"<<endl;
  
  cudaMemcpy(dom_d->not_write_surf_ID, not_write, solid_part_count * sizeof(bool), cudaMemcpyHostToDevice);
  delete not_write;
  
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////// INITIALIZE //////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    ///// REMAINS THERMAL
    
    
    // cout << "Reduction Type is: ";
    // if (dom.nonlock_sum)
      // cout << "Nishimura "<<endl;
    // else
      // cout << "Locking "<<endl;
		// //dom.SolveDiffUpdateLeapfrog(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    // if (solver=="Mech" || solver=="Mech-Thermal")
      // dom.SolveDiffUpdateFraser(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    // else if (solver=="Mech" || solver=="Mech-Thermal-KickDrift")
      // dom.SolveDiffUpdateKickDrift(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    // else if (solver=="Thermal")
      // dom.ThermalSolve(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    // else 
      // throw new Fatal("Invalid solver.");
		// } else {
      // throw new Fatal("Particle Count is Null. Please Check Radius and Domain Dimensions.");
    // }
		// dom.WriteXDMF("maz");

	cout << "Solving "<<endl;
	//CheckData<<<1,1>>>(dom_d);
	//cudaDeviceSynchronize(); //Crashes if not Sync!!!
	
	

	
	cout << "Time Step: "<<dom_d->deltat<<endl;
	//riteCSV("test_inicial.csv", x, dom_d->u_h, dom.Particles.size());
	//dom_d->MechSolve(0.00101 /*tf*//*1.01*/,1.e-4);
	//dom_d->MechSolve(100*timestep + 1.e-10 /*tf*//*1.01*/,timestep);
  

	//dom_d->MechSolve(0.0101,1.0e-4);
  
  //New solver
  dom_d->auto_ts = false;
  //timestep = (1.0*h/(Cs+VMAX));
  dom_d->deltat = 0.4*h/(Cs+VMAX);
  
  
  //dom_d->MechKickDriftSolve(0.0101,1.0e-4);
  //LEAPFROG IS WORKING WITH ALPHA = 1
  //KICKDRIFT IS NOT 
  //dom_d->MechLeapfrogSolve(0.0101,1.0e-4);
  //dom_d->MechFraserSolve(5*timestep,timestep);
  dom_d->MechFraserSolve(sim_time,output_time);

	dom_d->WriteCSV("test.csv");
	
	cudaFree(dom_d);
	//report_gpu_mem();
	cout << "Program ended."<<endl;

	}	//Argc > 0
  else {cout << "No input file found. Please specify input file."<<endl;}
	
    return 0;
}


