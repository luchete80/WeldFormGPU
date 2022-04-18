#include <iostream>
#include <cuda_runtime.h>

#define DOM_SIZE 1200

using namespace std;

double step_size = 1e-7;

__device__  void test(double *x, double ts){
  int i = threadIdx.x + blockDim.x*blockIdx.x;
  if ( i < DOM_SIZE)
    x[i] += 3.14159e-3*i*ts;
}

__global__ void test_kernel (double *x, double ts){
  test(x,ts);
}

int main() {
  
  int N = DOM_SIZE;
	int threadsPerBlock = 256; //Or BlockSize
	int blocksPerGrid =				// Or gridsize
	(N + threadsPerBlock - 1) / threadsPerBlock;
  
  
  double *x_d;
  cudaMalloc((void **)&x_d, 	DOM_SIZE * sizeof (double));
  //Init values
  
  double *x =  new double [DOM_SIZE]; //Only to compare
  double *x_h =  new double [DOM_SIZE];
  double3 *vx =  new double3 [DOM_SIZE];
  
  //Init values
	for (int i=0;i<DOM_SIZE;i++){
		//cout <<"i; "<<i<<endl;
		//x[i] = make_double3(dom.Particles[i]->x);
		//x[i] = make_double3(double(dom.Particles[i]->x(0)), double(dom.Particles[i]->x(1)), double(dom.Particles[i]->x(2)));
    x[i] = x_h [i]= 0.;
	}
  cudaMemcpy(x_d, x, sizeof(double) * DOM_SIZE, cudaMemcpyHostToDevice);	
  
  int dtout = 1000;
  int count = 1;
  for (int ts=0; ts < 1e7; ts++){
    
    for (int i=0;i<DOM_SIZE;i++)
      x_h[i]+= 3.14159e-3*i*step_size;
    
    test_kernel <<<blocksPerGrid,threadsPerBlock >>>(x_d, step_size);
    
    if (ts>=dtout*count){
      cudaMemcpy(x, x_d, sizeof(double) * DOM_SIZE, cudaMemcpyDeviceToHost);	
      double error = 0.;
      for (int i=0;i<DOM_SIZE;i++){
        error+= abs(x_h[i] - x[i])/*/x_h[i]*/;
        //cout << "x[i] x_h[i]"<<x[i]<<", "<<x_h[i]<<endl;
      }
      count ++;
      cout <<"Error sum "<<error<<endl;    
    }
  }

  
  
  
  
  cout << "test on CPU"<<endl;
	

}