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
#include "cuda/Mesh.cu"
#include "cuda/Boundary_Condition.cuh"

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
    
  // if (domi.contact){
    // for (int bc=0;bc<domi.bConds.size();bc++){
      // for (int m=0;m<domi.trimesh.size();m++){
        // if (domi.trimesh[m]->id == domi.bConds[bc].zoneId)
          // domi.trimesh[m]->SetVel(domi.bConds[bc].value);
      // }//mesh
    // }//bcs
  // }//contact

  if (domi.contact){
    for (int bc=0;bc<domi.bConds.size();bc++){
        if (domi.trimesh->id == domi.bConds[bc].zoneId)
          domi.trimesh->SetVel(domi.bConds[bc].value);
    }//bcs
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
    report_gpu_mem();
    gpuErrchk(cudaMallocManaged(&dom_d, sizeof(SPH::Domain)) );
    report_gpu_mem();
  
		dom.Dimension	= 3;
		
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
    
    
    string solver = "Mech";
    readValue(config["solver"],solver);
    
    if (solver=="Mech-Thermal")
      dom_d->thermal_solver = true;
		
		// readValue(config["integrationMethod"], dom.Scheme); //0:Verlet, 1:LeapFrog, 2: Modified Verlet

     	// //dom.XSPH	= 0.5; //Very important

    double dx,r,h;
    bool h_update = false;
		double hfactor;
    
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
    bool kernel_grad_corr = false;
    readValue(config["artifViscAlpha"],dom_d->Alpha); //TODO: ARTIFF VISC PER PARTICLE
    readValue(config["artifViscBeta"],dom_d->Beta);
    // readValue(config["contAlgorithm"],cont_alg);
    readValue(config["kernelGradCorr"],kernel_grad_corr);
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
    bool sym[] = {false,false,false};
		readValue(domblock[0]["id"], 	id);
		readVector(domblock[0]["start"], 	start);
		cout << "Reading Domain dim" << endl;  readVector(domblock[0]["dim"], 	L);
		cout << "Reading Domain type" << endl; readValue(domblock[0]["type"], 	domtype); //0: Box
    cout << "Reading Domain mat id" << endl;  readValue(domblock[0]["matID"], 	matID); //0: Box
    cout << "Grid Coordinate System" << endl;  readValue(domblock[0]["gridCoordSys"], 	gridCS); //0: Box
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
			dom.AddBoxLength(id ,start, L.x , L.y,  L.z , r ,rho, h, 1 , 0 , false, false );		
		}   else if (domtype == "Cylinder"){
      cout << "Adding Cylinder";      
			// if (sym[0] && sym[1]){
        // cout << " with symmetry..."<<endl;
        // dom.AddXYSymCylinderLength(0, L[0]/2., L[2], r, rho, h, false, sym[2]); 
      // }
      // else {
        // cout << "..."<<endl;
        // if ( gridCS == "Cartesian")
          cout << "Reserved "<<ComputeCylinderParticles (L.x/2., L.z, r)<<" particles."<<endl;
        cout << "Length " << L.x<<", "<< L.y<<", "<< L.z<<", "<<endl;
          dom.AddCylinderLength(0, start, L.x/2., L.z, r, rho, h, false);  /////// GENERATED AT HOST TO THEN COPY
        // else if (gridCS == "Cylindrical")
          // dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h);
          
      // }
    }
      // else {
        // cout << "..."<<endl;
        // if ( gridCS == "Cartesian")
          // dom.AddCylinderLength(0, start, L[0]/2., L[2], r, rho, h, false, sym[2]); 
        // else if (gridCS == "Cylindrical")
          // dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h);
          
      // }
    // }

        // cout <<"t  			= "<<timestep<<endl;
        // cout <<"Cs 			= "<<Cs<<endl;
        // cout <<"K  			= "<<E<<endl;
        // cout <<"G  			= "<<nu<<endl;
        // cout <<"Fy 			= "<<Fy<<endl;
		// cout <<"dx 			= "<<dx<<endl;
		// cout <<"h  			= "<<h<<endl;
		// cout <<"-------------------------"<<endl;
		// cout <<	"Dim: "<<dom.Dimension<<endl;				

    // for (int i=0;i<3;i++) {//TODO: Increment by Start Vector
			// dom.DomMax(0) = L[i];
			// dom.DomMin(0) = -L[i];
		// }		

		cout << "Particle count: "<<dom.Particles.size()<<endl;

    cout << "Domain Zones "<<domzones.size()<<endl;		
		for (auto& zone : domzones) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid;
      Vector vstart, vend;
			readValue(zone["id"], 		zoneid);
			readVector(zone["start"], 	vstart);
			readVector(zone["end"], 	vend);
      //cout << "Zone id "<<zoneid<<endl;
			// cout << "Dimensions: "<<endl;
			//cout << "start"<< vstart(0)<<"; "<< vstart(1)<<"; "<< vstart(2)<<"; "<<endl;
			//cout << "end"<< vend(0)<<"; "<< vend(1)<<"; "<< vend(2)<<"; "<<endl;
      
			int partcount =dom.AssignZone(vstart,vend,zoneid); ////IN DEVICE DOMAIN
      std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
		}
    
    // //////////////////////////////////////////////////////////
    // ////////////////// RIGID BODIES //////////////////////////
    string rigbody_type;
    bool contact = false;
    if (readValue(rigbodies[0]["type"],rigbody_type))
      contact = true;
    Vector dim;
    
		Vector vstart;
    readVector(rigbodies[0]["start"], 	vstart);       
		readVector(rigbodies[0]["dim"], 	dim); 
    bool flipnormals = false;
    readValue(rigbodies[0]["flipNormals"],flipnormals);
    
    // double heatcap = 1.;
    // readValue(rigbodies[0]["thermalHeatCap"],heatcap);
    // //TODO: WRitE TO PArTiclES
    if (rigbody_type == "File"){
      // string filename = "";
      // readValue(rigbodies[0]["fileName"], 	filename); 
      // cout << "Reading Mesh input file..." << endl;
      // SPH::NastranReader reader("Tool.nas", flipnormals);
    }
    else {
      if (dim (0)!=0. && dim(1) != 0. && dim(2) !=0. && rigbody_type == "Plane")
        printf("ERROR: Contact Plane Surface should have one null dimension\n");
    }
    std::vector<TriMesh *> mesh; ////// TODO: ALLOW FOR MULTIPLE MESH CONTACT
    SPH::TriMesh_d *mesh_d;
  
    cout << "Set contact to ";
    if (contact){
      cout << "true."<<endl;
      dom_d->contact = true;
      cout << "Reading contact mesh..."<<endl;
      //TODO: CHANGE TO EVERY DIRECTION
      int dens = 10;
      readValue(rigbodies[0]["partSide"],dens);
      if (rigbody_type == "Plane"){

        SPH::TriMesh mesh;
        //mesh[0]->AxisPlaneMesh(2, false, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), start(2)),dens);
        mesh.AxisPlaneMesh(2,false,vstart,
        Vector(vstart(0)+dim(0),vstart(1)+dim(1), vstart(2)),
        dens);

        double hfac = 1.1;
        dom_d->first_fem_particle_idx = dom.Particles.size(); // TODO: THIS SHOULD BE DONE AUTOMATICALLY
        int solid_count = dom.Particles.size(); //BEFORE ADDING CONTACT MESH
        
        dom.AddTrimeshParticles(mesh, hfac, 11); //TO SHARE SAME PARTICLE NUMBER
        dom_d->contact_surf_id = 11; //TO DO: AUTO! From Domain_d->AddTriMesh
        
        //TODO: Mesh has to be deleted
        SPH::TriMesh_d *mesh_d;
        gpuErrchk(cudaMallocManaged(&mesh_d, sizeof(SPH::TriMesh_d)) );
        mesh_d->AxisPlaneMesh(2,false,make_double3(vstart(0),vstart(1),vstart(2)),make_double3(vstart(0)+dim(0),vstart(1)+dim(1),vstart(2)),dens);
        
        cout << "Domain Size "<<dom.Particles.size()<<endl;
        //BEFORE ALLOCATING 
        int particlecount = dom.Particles.size();
        // //cout << "Particles "<<
        dom_d->SetDimension(particlecount);	 //AFTER CREATING DOMAIN
        dom_d->solid_part_count = solid_count;  //AFTER SET DIMENSION
        dom_d->trimesh = mesh_d; //TODO: CHECK WHY ADDRESS IS LOST
        if (dom_d->trimesh ==NULL)
          cout << "ERROR. No mesh defined"<<endl;
        
      } else if (rigbody_type == "File"){
        // string filename = "";
        // readValue(rigbodies[0]["fileName"], 	filename); 
        // cout << "Reading Mesh input file " << filename <<endl;
        // SPH::NastranReader reader(filename.c_str());
          // mesh.push_back (new SPH::TriMesh(reader,flipnormals ));
      }
      // cout << "Creating Spheres.."<<endl;
      // //mesh.v = Vec3_t(0.,0.,);
      // mesh[0]->CalcSpheres(); //DONE ONCE
      // double hfac = 1.1;	//Used only for Neighbour search radius cutoff
      // cout << "Adding mesh particles ...";
      // int id;
      // readValue(rigbodies[0]["zoneId"],id);
      // dom.AddTrimeshParticles(mesh[0], hfac, id); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
        
      
      std::vector<double> fric_sta(1), fric_dyn(1), heat_cond(1);
      readValue(contact_[0]["fricCoeffStatic"], 	fric_sta[0]); 
      readValue(contact_[0]["fricCoeffDynamic"], 	fric_dyn[0]); 
      // readValue(contact_[0]["heatCondCoeff"], 	  heat_cond[0]);
      
      // bool heat_cond_ = false;
      // if (readValue(contact_[0]["heatConductance"], 	heat_cond_)){
        // dom.cont_heat_cond = true;
        // dom.contact_hc = heat_cond[0];
      // }
      
      dom_d->friction_dyn = fric_dyn[0];
      dom_d->friction_sta = fric_sta[0];
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
			
      std::cout<< "BCs "<<  ", Zone ID: "<<bcon.zoneId<<", Value :" <<bcon.value.x<<", "<<bcon.value.y<<", "<<bcon.value.z<<std::endl;
		}//Boundary Conditions
		
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
    dom_d->gradKernelCorr = kernel_grad_corr;
    // dom.ts_nb_inc = 5;
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// TODO: THIS SHOULD BE PASSED TO A DOM_D INTIIALIZE FUNCTION ///////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  dom_d->GeneralAfter = & UserAcc;
	dom_d->SetDimension(dom.Particles.size());	 //AFTER CREATING DOMAIN
  //SPH::Domain	dom;
	//double3 *x =  (double3 *)malloc(dom.Particles.size());
	double3 *x =  new double3 [dom.Particles.size()];
	for (int i=0;i<dom.Particles.size();i++){
		//cout <<"i; "<<i<<endl;
		//x[i] = make_double3(dom.Particles[i]->x);
		x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
	}
	int size = dom.Particles.size() * sizeof(double3);
	cout << "Copying to device "<<dom.Particles.size() << " particle properties ..."<<endl;
	cudaMemcpy(dom_d->x, x, size, cudaMemcpyHostToDevice);


	for (int i=0;i<dom.Particles.size();i++){
		x[i] = make_double3(0.,0.,0.);
	}
	cudaMemcpy(dom_d->v, x, size, cudaMemcpyHostToDevice);
  
	cout << "copied"<<endl;

  ////// MATERIAL  
  Material_ *mat_h = (Material_ *)malloc(dom_d->solid_part_count * sizeof(Material_ *)); 
  

  
  Elastic_ el(E,nu);
  cout << "Mat type  "<<mattype<<endl;
  
  
  //cudaMalloc(&dom_d->materials, sizeof(Material_ ));
  
  cudaMalloc((void**)&dom_d->materials_ptr, sizeof(Material_ *)); //https://forums.developer.nvidia.com/t/virtual-funtions-in-kernels/22117/3
  Material_ *material_h; //(Material_ *)malloc(1 * sizeof(Material_ ));
  
    Material_ mat(el);
    mat.Ep = 10.0;
  //TODO: MATERIALS SHOULD BE A VECTOR   
  if      (mattype == "Bilinear")    {
    if (c.size()>0)
      Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
    else
      cout << "ERROR. MATERIAL CONSTANTS UNDEFINED"<<endl;
    cout << "Material Ep"<<Ep<<endl;
    material_h  = new Material_(el);
    material_h->Ep = Ep;
    material_h->Material_model = BILINEAR;
    // cout << "Material Constants, Et: "<<c[0]<<endl;
    // material_h->Material_model = BILINEAR;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Bilinear )); //
    // cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Bilinear), cudaMemcpyHostToDevice);	

    cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_ )); //
     cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);	
    
  } 
  else if (mattype == "Hollomon")    {
    // material_h  = new Hollomon(el,Fy,c[0],c[1]);
    // cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    // cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Hollomon));
    
    material_h  = new Material_(el);
    material_h->InitHollomon(el,Fy,c[0],c[1]);
    material_h->Material_model = HOLLOMON;
    cudaMalloc((void**)&dom_d->materials, 1 * sizeof(Material_));
    
    //init_hollomon_mat_kernel<<<1,1>>>(dom_d); //CRASH
    //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Hollomon*), cudaMemcpyHostToDevice);	
    cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(Material_), cudaMemcpyHostToDevice);	 //OR sizeof(Hollomon)??? i.e. derived class
    
  
  } else if (mattype == "JohnsonCook") {
    //Order is 
                               //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
   //Material_ *material_h  = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
    
    //Only 1 material to begin with
    //cudaMalloc((void**)&dom_d->materials, 1 * sizeof(JohnsonCook ));
    //cudaMemcpy(dom_d->materials, material_h, 1 * sizeof(JohnsonCook), cudaMemcpyHostToDevice);	
    cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
  } else                              printf("ERROR: Invalid material type.");
    
  // InitMatHollomonKernel<<<1,1>>>(dom_d,*material_h); //BECAUSE MEMCPY IS NOT WORKING
  // cudaDeviceSynchronize(); 
	
	cout << "Setting values"<<endl;
	dom_d->SetDensity(rho);
	dom_d->Set_h(h);
	cout << "done."<<endl;

	double *m =  new double [dom.Particles.size()];
	for (size_t a=0; a<dom.Particles.size(); a++)
		m[a] = dom.Particles[a]->Mass;
	cudaMemcpy(dom_d->m, m, dom.Particles.size() * sizeof(double), cudaMemcpyHostToDevice);	
  
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


