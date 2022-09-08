#include "Domain_d.cuh"

///////////////////////////////// NO INLINE FUNCTIONS
namespace SPH{

int calcHalfPartCount_c(const double &r, const double &R, const int xinc){
	int ypartcount = -1;
	if ( xinc > 0 ){
		ypartcount = 1;
		double yp = r;
		double xp = r + (double)(xinc - 1 ) *2.*r; 
		double rad = sqrt(yp*yp + xp*xp);
		while( rad <= R -r ){
			yp += 2.*r;
			rad = sqrt(yp*yp + xp*xp);
			ypartcount++;
		}
		ypartcount-=1;
	}
	return ypartcount;
}

void Domain_d::AddCylinderLength(int tag, Vector const & V, double Rxy, double Lz, 
									double r, double Density, double h, bool Fixed) {

//	Util::Stopwatch stopwatch;
    //std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    //size_t PrePS = Particles.size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy=1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount_c(r, Rxy, 1);
	
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	

    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2);

		while (zp <= (V(2)+Lz-r)) {
			j = 0;
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount_c(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//Particles.push_back(new Particle(tag,Vector(xp,yp,zp),Vector(0,0,0),0.0,Density,h,Fixed));
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
			zp = V(2) + (2.0*k+1)*r;
		}
		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.size()-PrePS);
		//double Mass = Vol * Density /Particles.size();

		// for (int i=0; i<Particles.size(); i++)//Like in Domain::Move
		// {
			// Particles[i]->Mass = Mass;
		// }


}

}; //SPH