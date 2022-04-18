//TODO:!!! Pending contact
//contactforce = 0
//In global initialize 
// domain max_contact_force

namespace SPH{




//Inputs
//max_contact_force
//Neighbours
//vectors v
void __device__ inline Domain_d::CalcContactForces(/* int i*/){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
	
	int neibcount = particlenbcount[i];
	const uint writeOffset = neighborWriteOffsets[i];

	// printf("neibcount %d\n",neibcount);
	// printf("Nb indexed,i:%d\n",i);
	// In this Weldform GPU version, is clear than i is SPH particle and 2 is RIGID PARTICLE
	for (int k=0;k < neibcount;k++) { //Or size

    int j = neighbors[writeOffset + k];
		double3 xij;
		double h,K;
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
			Element* e = trimesh-> element[element[j]];
			double pplane = e -> pplane; 
			//cout<< "contact distance"<<Particles[P1]->h + pplane - dot (normal[j],	x[i])<<endl;

			double deltat_cont = ( h[i] + pplane - dot (normal[j],	x[i]) ) / (-delta_);								//Eq 3-142 
			double3 Ri = x[i] + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction

			//Check for contact in this time step 
			//Calculate time step for external forces
			double dt_fext = contact_force_factor * (m[i] * 2. * length(v[i]) / length(contforce[i]) );	//Fraser 3-145
			
			contforce[i] = 0.; //RESET
			
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
				int i=0,j;			
				while (i<3 && inside){
					j = i+1;	if (j>2) j = 0;
					double crit = dot (cross ( *trimesh->node[e -> node[j]] - *trimesh->node[e -> node[i]],
																															Qj  - *trimesh->node[e -> node[i]]),
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
					double kij = PFAC * Particles[P1]-> cont_stiff;
					double omega = sqrt (kij/m[i]);
					double psi_cont = 2. * m[i] * omega * DFAC; // Fraser Eqn 3-158
										
					// TANGENTIAL COMPONENNT
					// Fraser Eqn 3-167
					// TODO - recalculate vr here too!
					double3 tgvr, tgdir;
					if (friction > 0.) {						
						if ( norm (vr)  != 0.0 ) {
							//TODO: THIS VELOCITY SHOULD BE THE CORRECTED ONE 
							//double3 tgvr  = vr - dot(vr,normal[j]) * normal[j];
							// Is Fraser thesis is explained better 
							double3 tgvr = vr + delta_ * normal[j];  // -dot(vr,normal) * normal
							double3 tgdir = tgvr / length(tgvr);
						}
					}

					contforce[i] = (kij * delta - psi_cont * delta_) * normal[j]; // NORMAL DIRECTION, Fraser 3-159
					double force2 = dot(contforce[i],contforce[i]);
					
					// if (force2 > (1.e10))
						// contforce[i] = 1.e5;
					dt_fext = contact_force_factor * (m[i] * 2. * norm(v[i]) / norm (contforce[i]));

					if (dt_fext < min_force_ts_){
						min_force_ts_ = dt_fext;
						if (dt_fext > 0)
							this -> min_force_ts = min_force_ts_;
					}
					a[i] += contforce[i] / m[i]; 
					//cout << "contforce "<<contforce[i]<<endl;
					
					if (friction > 0.) {
						if ( norm (vr)  != 0.0 ){
						// //TG DIRECTION
							double3 tgforce = friction * norm(contforce[i]) * tgdir;
							a[i] += tgforce / m[i]; 
							//cout << "tg force "<< tgforce <<endl;
						}
					}
					
					if   (force2 > max_contact_force ) max_contact_force = force2;
					else if (force2 < min_contact_force ) min_contact_force = force2;
					inside_pairs++;
				}// if inside
			} //deltat <min

		}//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER

	}//neibcount

	max_contact_force = sqrt (max_contact_force);
	min_contact_force = sqrt (min_contact_force);
	if (max_contact_force > 0.){
		//cout << "Min Contact Force"<< min_contact_force<<"Max Contact Force: "<< max_contact_force << "Time: " << Time << ", Pairs"<<inside_pairs<<endl;
		//cout << " Min tstep size: " << min_force_ts << ", current time step: " << deltat <<endl;
		//TEMP
		// if (min_force_ts> 0)
			// deltat = min_force_ts;
	}
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

};//SPH