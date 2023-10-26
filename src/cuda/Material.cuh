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
  double Ep;  //If Bilinear this is constant, 
  int			Material_model;	//TODO: Change to enum

  ////// NO //virtual FUNCTIONS, HOLLOMON MATERIAL ///////
  double K, m;
	double eps0, eps1;
  double sy0;
  void InitHollomon(){}

  ////// NO //virtual FUNCTIONS, JOHNSON COOK MATERIAL ///////
	double T_t,T_m;	//transition and melting temps
	double A, B, C;
	double n/*, m*/;
	double eps_0;


  __host__ __device__ void InitHollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_)
  {
    elastic_m = el;
    K = k_;
    m = m_;
  Material_model = HOLLOMON;
  eps0 = sy0_/el.E(); 
  sy0  = sy0_;
  eps1 = pow(sy0_/k_, 1./m);
  printf( "eps_0, %.2e eps_1 %.2e\n",eps0, eps1);
  if (eps0 > eps1){
    printf ("ERROR, Hollomon material bad definition, please correct Yield Stress, Elastic Modulus or Material hardening constants.");
  }
  
  }
  void Init_JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              const double &c, const double &eps_0_,
              const double &m_, const double &T_m_, const double &T_t_)
  //:
	// Material_(el),A(a),B(b),C(c),
  // m(m_),n(n_),eps_0(eps_0_),T_m(T_m_),T_t(T_t_)
  { 
    elastic_m = el;
    Material_model = JOHNSON_COOK;
    A=a;B=b;C=c;
    m=m_;n=n_;eps_0=eps_0_;
    T_m=T_m_;T_t=T_t_;
  }

	Material_(){}
  __host__ __device__ void test(){printf("test\n");}
//  //virtual __host__ __device__ double testret(){return 2.0;}
	Material_(const Elastic_ el):elastic_m(el){}
	// //virtual  __device__ inline double CalcTangentModulus(){};
	// //virtual  __device__ inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp){};
	// //virtual  __device__ inline double CalcTangentModulus(const double &strain){};
	// //virtual  __device__ inline double CalcYieldStress(){};
	// //virtual  __device__ inline double CalcYieldStress(const double &strain){return 0.0;};
	// //virtual  __device__ inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){};
	 __host__ __device__ const Elastic_& Elastic()const{return elastic_m;}
};

class _Plastic{
	
	public:
	//virtual inline double CalcYieldStress();	
	////virtual inline double CalcYieldStress();
};

class Bilinear:
public Material_{

 	public:
	Bilinear(const double &ep){ //THIS IS DIFFERENT FROM WELDFORM CPU, IN WHICH BILINEAR IS NOT A MATERIAL
    Ep = ep;
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
	inline double  __device__ CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
	inline double  __device__ CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp);
  //~JohnsonCook(){}
};

class Hollomon:
public Material_{
	double K, m;
	double eps0, eps1;
  double sy0;
	
	public:
	__host__ __device__ Hollomon(){Material_model = HOLLOMON;}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	//ASSUMING AT FIRST COEFFICIENTS ARE GIVEN TO TOTAL STRAIN-STRESS
	__host__ __device__ Hollomon(const double eps0_, const double &k_, const double &m_):
    K(k_), m(m_){ 
    eps0 = eps0_;
    Material_model = HOLLOMON;}
  __host__ __device__  Hollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_):
  Material_(el),K(k_), m(m_) {
  Material_model = HOLLOMON;
  eps0 = sy0_/el.E(); 
  sy0  = sy0_;
  eps1 = pow(sy0_/k_, 1./m);
  printf( "eps_0, %.2e eps_1 %.2e\n",eps0, eps1);
  if (eps0 > eps1){
    printf ("ERROR, Hollomon material bad definition, please correct Yield Stress, Elastic Modulus or Material hardening constants.");
  }
  
  }
  __host__ __device__ double testret(){printf("hollomon testret\n"); return 2.0;}
	
  inline double  __device__ CalcTangentModulus(const double &strain);
	inline double  __device__ CalcYieldStress(){}	
	inline double  __device__ CalcYieldStress(const double &strain);	
};

__device__ inline double CalcHollomonYieldStress(const double &strain, Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
  double sy = 0.0;
  //printf("K %f eps0 %f , eps1 %f sy0 %f m %f\n",K,eps0,eps1,sy0,m);
  //printf("K %f ", mat->K);
   if (strain + mat->eps0 > mat->eps1) sy = mat->K*pow(strain + mat->eps0,mat->m); //plateau surpassed. If no plateau, eps1=eps0 so 
   else                      sy = mat->sy0; 
	return sy; 
  
}

__device__ inline void ShowProps(Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
  printf("K %f \n eps0 %f \n eps1 \n%f, sy0 \n%f, m %f \n", mat->K,mat->eps0,mat->eps1,mat->sy0,mat->m);
  
}

__device__ inline double CalcJohnsonCookYieldStress(const double &strain, const double &strain_rate, const double &temp, Material_ *mat) //IN CASE OF NOT USING //virtual FUNCTIONS
{
	double T_h = (temp - mat->T_t) / (mat->T_m - mat->T_t);
	double sr = strain_rate;
	if (strain_rate == 0.0)
		sr = 1.e-5;
	
	double sy = (mat->A+mat->B*pow(strain, mat->n))*(1.0 + mat->C * log (sr/ mat->eps_0) ) * (1.0 - pow(T_h,mat->m));
	
	return sy;
}

inline double __device__ CalcHollomonTangentModulus(const double &strain, Material_ *mat) {
	double Et;
  if (strain + mat->eps0 > mat->eps1) Et = mat->K * mat->m * pow(strain + mat->eps0, (mat->m-1.0));
  else                      Et = 0.;
	//cout << "ET: "<<Et<<endl;
	return Et;
}

inline double __device__ CalcJohnsonCookTangentModulus(const double &plstrain, const double &strain_rate, const double &temp, Material_ *mat)	{
	double sy, T_h;
  //cout << "n, B, C, eps_0, T_t, m"<<n<<", "<<B<<", "<<C<<"eps0, "<<eps_0<<", "<<", "<<T_t<<", "<<m<<endl;
	T_h = (temp - mat->T_t) / (mat->T_m - mat->T_t);
	
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

  if (plstrain > 0.)
    Et = mat->n * mat->B * pow(plstrain,mat->n-1.)*(1.0 + mat->C*log(strain_rate/ mat->eps_0)) * (1.0-pow (T_h,mat->m));
  else 
    Et = mat->Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
  return Et;
}	

//#include "Material.cu"

#endif