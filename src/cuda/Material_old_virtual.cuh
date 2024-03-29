#ifndef _MATERIAL_CUH_
#define _MATERIAL_CUH_

class Particle;

#define BILINEAR				0
#define HOLLOMON				1 //POWER LAW
#define JOHNSON_COOK		2

class Elastic_{
	private:
	double E_m, nu_m;	//Poisson and young
	double K_m, G_m;
	
	public:
	__host__ __device__ Elastic_(){}
	__host__ __device__ Elastic_(const double &e, const double &nu):E_m(e),nu_m(nu){}
	__host__ __device__ const double& E()const{return E_m;}
	
};

class Material_{
	
	protected:
	Elastic_ elastic_m;
	double E_m, nu;	//TODO, move to elastic class
	public:
  int			Material_model;	//TODO: Change to enum
	Material_(){}
	Material_(const Elastic_ el):elastic_m(el){}
	virtual __device__ inline double CalcTangentModulus(){};
	virtual __device__ inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp){};
	virtual __device__ inline double CalcTangentModulus(const double &strain){};
	virtual __device__ inline double CalcYieldStress(){};
	virtual __device__ inline double CalcYieldStress(const double &strain){return 0.0;};
	virtual __device__ inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){};
	__host__ __device__ const Elastic_& Elastic()const{return elastic_m;}
};

class _Plastic{
	
	public:
	virtual inline double CalcYieldStress();	
	//virtual inline double CalcYieldStress();
};

class Bilinear:
public Material_{
 	public:
	Bilinear(){ //THIS IS DIFFERENT FROM WELDFORM CPU, IN WHICH BILINEAR IS NOT A MATERIAL
    Material_model = BILINEAR;
  }
};


//TODO: derive johnson cook as plastic material flow
class JohnsonCook:
public Material_{
	double T_t,T_m;	//transition and melting temps
	double A, B, C;
	double n, m;
	double eps_0;
	
	public:
	JohnsonCook(){
    Material_model = JOHNSON_COOK;
  }
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              const double &c, const double &eps_0_,
              const double &m_, const double &T_m_, const double &T_t_):
	Material_(el),A(a),B(b),C(c),
  m(m_),n(n_),eps_0(eps_0_),T_m(T_m_),T_t(T_t_)
  { Material_model = JOHNSON_COOK;
  }
	inline double __device__ CalcYieldStress(){return 0.0;}	
	inline double __device__ CalcYieldStress(const double &plstrain){
     double Et =0.;

    if (plstrain > 0.)
      Et = n * B * pow(plstrain,n-1.);
    else 
      Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
    return Et;
  } //TODO: SEE IF INCLUDE	
	inline double __device__ CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
	inline double __device__ CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp);
  //~JohnsonCook(){}
};

class Hollomon:
public Material_{
	double K, m;
	double eps0, eps1;
  double sy0;
	
	public:
	Hollomon(){}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	//ASSUMING AT FIRST COEFFICIENTS ARE GIVEN TO TOTAL STRAIN-STRESS
	__device__ __host__ Hollomon(const double eps0_, const double &k_, const double &m_):
	K(k_), m(m_){ eps0 = eps0_;}
		__device__ __host__  Hollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_){};
	inline double __device__ CalcTangentModulus(const double &strain);
	inline double __device__ CalcYieldStress(){}	
	inline double __device__ CalcYieldStress(const double &strain);	
};

#include "Material.cu"

#endif