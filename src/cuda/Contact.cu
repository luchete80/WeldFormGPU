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


__global__ void CalcContactForcesKernel(Domain_d *dom_d,	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors,
                                  double *cont_forces
																	/*const int &id, const double &totmass*/){
	dom_d->CalcContactForces(
	particlenbcount,
	neighborWriteOffsets,
	neighbors,
  cont_forces);

}
               
//Inputs
//max_contact_force
//Neighbours
//vectors v
void __device__ inline Domain_d::CalcContactForces(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors,
                                                  double *cont_forces){
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
      int e = element[j]; //Index of mesh Element associated with node
      // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
      //IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
      //cout << "Contact pair size: "<<ContPairs[k].Size()<<endl;

      double3 vr = v[i] - v[j];		//Fraser 3-137
      //cout << "Particle P1v: "<<v[i]<<endl;
      //cout << "Particle P2v: "<<Particles[P2]->v<<endl;
      //ok, FEM Particles normals can be calculated by two ways, the one used to
      //calculate SPH ones, and could be given by mesh input
      //delta_ Is the projection of relative velocity 
      double delta_ = - dot( normal[j] , vr);	//Penetration rate, Fraser 3-138
      

      //Check if SPH and fem particles are approaching each other
      if (delta_ > 0 ){
        //HERE ELEMENT IS NOT AN OBJECT, THIS IS REPRESENTED AS SOA IN MESH CLASS
        //EACH "RIGID" PARTICLE HAS AN ELEMENT ASOCIATED
        double pplane = trimesh->pplane[e]; 
        //cout<< "contact distance"<<Particles[P1]->h + pplane - dot (normal[j],	x[i])<<endl;

        double deltat_cont = ( h[i] + pplane - dot (normal[j],	x[i]) ) / (-delta_);								//Eq 3-142 
        double3 Ri = x[i] + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction
        //printf("delta %f \n",delta_);
        
        contforce[i] = make_double3(0.); //RESET
        
        // if (dt_fext > deltat)
          // cout << "Time step size ("<<deltat<<" is larger than max allowable contact forces step ("<<dt_fext<<")"<<endl;
        if (deltat_cont < deltat){ //Originaly //	if (deltat_cont < std::min(deltat,dt_fext) 
        
          //cout << "Inside dt contact" <<endl;
          //Find point of contact Qj
          double3 Qj = x[i] + (v[i] * deltat_cont) - ( h[i] * normal[j]); //Fraser 3-146
          //Check if it is inside triangular element
          //Find a vector 
          //Fraser 3-147
          bool inside = true;
          int l=0,n;			
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            // double crit = dot (cross ( *trimesh->node[e -> node[j]] - *trimesh->node[e -> node[i]],
                                                                // Qj  - *trimesh->node[e -> node[i]]),
                              // normal[j]);
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*e+n]] - trimesh->node[trimesh->elnode[3*e+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*e+l]]),
                              normal[j]);
            if (crit < 0.0) inside = false;
            i++;
          }
          
          if (inside ) { //Contact point inside element, contact proceeds

            //Calculate penetration depth (Fraser 3-49)
            double delta = (deltat - deltat_cont) * delta_;
            //cout << "delta: "<<delta<<endl;

            // DAMPING
            //Calculate SPH and FEM elements stiffness (series)
            //Since FEM is assumed as rigid, stiffness is simply the SPH one 
            double kij = PFAC * cont_stiff[i];
            double omega = sqrt (kij/m[i]);
            double psi_cont = 2. * m[i] * omega * DFAC; // Fraser Eqn 3-158
                      
            // TANGENTIAL COMPONENNT
            // Fraser Eqn 3-167
            // TODO - recalculate vr here too!
            double3 tgvr = vr + delta_ * normal[j];  // -dot(vr,normal) * normal
            double3 tgdir = tgvr / length(tgvr);

            contforce[i] = (kij * delta - psi_cont * delta_) * normal[j]; // NORMAL DIRECTION, Fraser 3-159
            double force2 = dot(contforce[i],contforce[i]);
            
            // if (force2 > (1.e10))
              // contforce[i] = 1.e5;
            double dt_fext = contact_force_factor * (m[i] * 2. * length(v[i]) / length(contforce[i]));////Fraser 3-145

            if (dt_fext < min_force_ts_){
              min_force_ts_ = dt_fext;
              if (dt_fext > 0)
                this -> min_force_ts = min_force_ts_;
            }
            a[i] += contforce[i] / m[i]; 
            //printf("contforce %f\n",contforce[i].x);
            
            if (friction_dyn > 0.) {
              if ( length (vr)  != 0.0 ){
              // //TG DIRECTION
                double3 tgforce = friction_dyn * length(contforce[i]) * tgdir;
                a[i] -= tgforce / m[i]; 
                //printf("tg force %lf\n", tgforce.z); 
              }
            }
            
            if   (force2 > max_contact_force ) max_contact_force = force2;
            else if (force2 < min_contact_force ) min_contact_force = force2;
            inside_pairs++;
          }// if inside
        } //deltat <min

      }//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER

    }//neibcount	for (int k=0;k < neibcount;k++) { //Or size
  
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

};//SPH