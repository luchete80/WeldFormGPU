//TODO:!!! Pending contact
//contactforce = 0
//In global initialize 
// domain max_contact_force
//#include "Mesh.cuh"

#define MAX_NB_COUNT    20

#include "Mesh.cuh"

namespace SPH{
__global__ inline void CalcContactNbKernel(Domain_d *dom_d,
const uint *particlenbcount,
const uint *neighborWriteOffsets,
const uint *neighbors) {
 	
  dom_d->CalcContactNb(
	particlenbcount,
	neighborWriteOffsets,
	neighbors); 
  
}

//Up until now max nb count is 20, this is done in order to fix contneib_part
//and to not reallocate dynamically each nb search
__device__ inline void Domain_d::CalcContactNb(const uint *particlenbcount,
const uint *neighborWriteOffsets,
const uint *neighbors){
  
  int i = threadIdx.x + blockDim.x*blockIdx.x;	
  if (i < first_fem_particle_idx ) {
    contneib_count[i] = 0;
    int neibcount = particlenbcount[i];
    const uint writeOffset = neighborWriteOffsets[i];
    
    for (int k=0;k < neibcount;k++) { //Or size
      int j = neighbors[writeOffset + k];
      double3 xij = x[i]-x[j];
      //double h = h[i] + h[j];  //not necessary to take average

      if ( (ID[i] == id_free_surf && ID[j] == contact_surf_id) /*||
           (ID[j] == id_free_surf && ID[i] == contact_surf_id) */) {
        /////if ( norm (Particles[P1]->x - Particles[P2]->x) < ( Particles[P1]->h + Particles[P2]->h ) ){ //2* cutoff being cutoff (h1+h2)/2
        if ( length(xij) < 2.*h[i] ){
          contneib_part[MAX_NB_COUNT*i + contneib_count[i]] = j;
          contneib_count[i]++;
        }
      } //IDs OK
    }//for k neighbours
  }// i < fem index
}

inline void __global__ UpdateContactParticlesKernel(Domain_d *dom){
  dom->UpdateContactParticles();
}

//inline void __device__ Domain_d::UpdateContactParticles(){
  // for (int m=0; m<trimesh.size();m++){
    // for ( int e = 0; e < trimesh[m]->element.Size(); e++ ){
      // Vec3_t v = 0.;
      // for (int en = 0;en<3;en++)
        // v += *trimesh[m] -> node_v[trimesh[m]->element[e] ->node[en]];
      // Particles[first_fem_particle_idx[m] + e] -> v = 
      // Particles[first_fem_particle_idx[m] + e] -> va = 
      // Particles[first_fem_particle_idx[m] + e] -> vb = v/3.;
      // Particles[first_fem_particle_idx[m] + e] -> a = 0.; 
      // Particles[first_fem_particle_idx[m] + e] -> normal  = trimesh[m]->element[e] -> normal;
      // //cout << "v "<< v/3.<<", n "<<Particles[first_fem_particle_idx + e] -> normal<<endl;
    // } 
  // }
//}

//TODO: CHANGE TO SEVERAL CONTACT SURFACES 
inline void __device__ Domain_d::UpdateContactParticles(){

  int e = threadIdx.x + blockDim.x*blockIdx.x;	
  if (e < trimesh->elemcount) {

    //int e = element[i];
    double3 vv = make_double3(0.);
    for (int en = 0; en<3; en++){
      //printf("particle %d \n",i);
      //printf ("node %d \n",trimesh->elnode[3*e+en]);
      // if (trimesh->elnode[3*e+en] < trimesh->nodecount )
      vv += trimesh -> node_v[trimesh->elnode[3*e+en]];
      // else 
        // printf("error \n");
    }
    
    //printf(" particle %d , v %f %f %f \n", i, vv.x, vv.y, vv.z);
    v [first_fem_particle_idx + e] = vv/3.;
    a [first_fem_particle_idx + e] = make_double3(0.);
    normal[first_fem_particle_idx + e] = trimesh -> normal[e];
    //printf("mesh normal, %f %f %f\n", trimesh -> normal[e].x , trimesh -> normal[e].y, trimesh -> normal[e].z);
  }
}


void __global__ CalcContactForcesKernel(Domain_d *dom_d,	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
                                  double3 *cont_forces
																	/*const int &id, const double &totmass*/){
	dom_d->CalcContactForcesWang(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
  cont_forces);

}
//WANG ET AL, FRICTION COEFFICIENT APPROACH    
//Inputs
//max_contact_force
//Neighbours
//vectors v
void __device__ inline Domain_d::CalcContactForcesWang(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors,
                                                  double3 *cont_forces){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  if (i < first_fem_particle_idx ) {  
    // CONTACT OFFSET IS FIX BY NOW
    int neibcount = contneib_count[i];
  
    //printf("i, first fem part, neibcount %d\n",neibcount);
    // printf("Nb indexed,i:%d\n",i);
    // In this Weldform GPU version, is clear than i is SPH particle and 2 is RIGID PARTICLE
    for (int k=0;k < neibcount;k++) { //Or size
      int j = contneib_part[i*MAX_NB_COUNT+k];
      //printf("j neib %d: \n",j);
      //int j = neighbors[writeOffset + k];
      double3 xij;
      double K;
      //int e = element[j]; //Index of mesh Element associated with node
      
      //TODO: MAKE ONE PER MESH
      int e = j - first_fem_particle_idx;
      //printf ("Element e %d\n", e);
      // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
      //IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
      //cout << "Contact pair size: "<<ContPairs[k].Size()<<endl;

      double3 vr = v[i] - v[j];		//Fraser 3-137
      //printf("particle 2 vr %f %f %f\n",v[j].x, v[j].y,v[j].z);
      // if (k==0)
        // printf("particle 2 (%d) x %f %f %f\n",j, x[j].x, x[j].y, x[j].z);
      //cout << "Particle P1v: "<<v[i]<<endl;
      //cout << "Particle P2v: "<<Particles[P2]->v<<endl;
      
      //ok, FEM Particles normals can be calculated by two ways, the one used to
      //calculate SPH ones, and could be given by mesh input
      //delta_ Is the projection of relative velocity 
      // double delta_ = - dot( normal[j] , vr);	//Penetration rate, Fraser 3-138
      
      double3 x_pred = x[i] + v[i] * deltat + a[i] * deltat * deltat/2.;
      //printf("xpred %f %f %f\n", x_pred.x, x_pred.y, x_pred.z);
      //printf ("pplane: %f", trimesh->pplane[e]);
      double dist = dot (normal[j],x_pred)  - trimesh->pplane[e];
      //printf("normal j %d %f %f %f\n", j, normal[j].x, normal[j].y, normal[j].z);
      //printf("dist %f\n", dist);
      if (dist < h[i]) {

        // double deltat_cont = ( h[i] + pplane - dot (normal[j],	x[i]) ) / (-delta_);								//Eq 3-142 
        // double3 Ri = x[i] + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction
        // //printf("delta %f \n",delta_);
        
        contforce[i] = make_double3(0.); //RESET
        //cont_forces[i] = make_double3(0.,0.,0.); //RESET

          //cout << "Inside dt contact" <<endl;
          /////Find point of contact Qj
          /////This is in Fraser algorithm
          double3 Qj = x[i] - dist * normal[j];
          // //Check if it is inside triangular element
          // //Find a vector 
          // //Fraser 3-147
          bool inside = true;
          int l=0,n;		   
          //printf("Entering while \n");
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            // double crit = dot (cross ( *trimesh->node[e -> node[j]] - *trimesh->node[e -> node[i]],
                                                                // Qj  - *trimesh->node[e -> node[i]]),
                              // normal[j]);
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*e+n]] - trimesh->node[trimesh->elnode[3*e+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*e+l]]),
                              normal[j]);
            if (crit < 0.0) inside = false;
            l++;
          }
          //printf("Outside while\n");
          
          if (inside ) { //Contact point inside element, contact proceeds

            // //Calculate penetration depth (Fraser 3-49)
            double delta = h[i] - dist;
            //printf("delta: %f\n", delta);
            // // DAMPING
            // //Calculate SPH and FEM elements stiffness (series)
            // //Since FEM is assumed as rigid, stiffness is simply the SPH one 
            double kij = 2.0 * m[i] / (deltat * deltat);
            //printf("deltat, kij %f %f\n", deltat, kij); 
            //double omega = sqrt (kij/m[i]);
            double psi_cont = kij * delta; // Fraser Eqn 3-158
                      
            // // TANGENTIAL COMPONENNT
            // // Fraser Eqn 3-167
            // // TODO - recalculate vr here too!
            // double3 tgvr = vr + delta_ * normal[j];  // -dot(vr,normal) * normal
            // double3 tgdir = tgvr / length(tgvr);

            contforce[i] = (kij * delta /*- psi_cont * delta_*/) * normal[j]; // NORMAL DIRECTION, Fraser 3-159
            //cont_forces[i] = contforce[i];
            
            // double force2 = dot(contforce[i],contforce[i]);
            
            // // if (force2 > (1.e10))
              // // contforce[i] = 1.e5;
            // double dt_fext = contact_force_factor * (m[i] * 2. * length(v[i]) / length(contforce[i]));////Fraser 3-145

            // if (dt_fext < min_force_ts_){
              // min_force_ts_ = dt_fext;
              // if (dt_fext > 0)
                // this -> min_force_ts = min_force_ts_;
            // }
            //printf("contforce %f %f %f ", contforce[i].x,contforce[i].y,contforce[i].z);
            a[i] += contforce[i] / m[i]; 
            //printf("contforce %f\n",contforce[i].x);
            
            // if (friction_dyn > 0.) {
              // if ( length (vr)  != 0.0 ){
              // // //TG DIRECTION
                // double3 tgforce = friction_dyn * length(contforce[i]) * tgdir;
                // a[i] -= tgforce / m[i]; 
                // //printf("tg force %lf\n", tgforce.z); 
              // }
            // }
            
            // if   (force2 > max_contact_force ) max_contact_force = force2;
            // else if (force2 < min_contact_force ) min_contact_force = force2;
            // inside_pairs++;
        }// if inside

      } //dist <h
    }//neibcount	for (int k=0;k < neibcount;k++) { //Or size
  
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

};//SPH