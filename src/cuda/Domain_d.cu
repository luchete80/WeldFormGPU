#include "Domain_d.cuh"
//Allocating from host
namespace SPH {
void Domain_d::SetDimension(const int &particle_count){
	//Allocae arrays (as Structure of arryays, SOA)

	cudaMalloc((void **)&x, particle_count * sizeof (Vector));
	cudaMalloc((void **)&v, particle_count * sizeof (Vector));
	cudaMalloc((void **)&a, particle_count * sizeof (Vector));

	cudaMalloc((void **)&T, particle_count * sizeof (Vector));
}


Domain_d::~Domain_d(){
	
		cudaFree(a);		
		cudaFree(v);
	
}

};//SPH
    // // Create host pointer to array-like storage of device pointers
    // Obj** h_d_obj = (Obj**)malloc(sizeof(Obj*) * 3); //    <--------- SEE QUESTION 1
    // for (int i = 0; i < 3; i++) {
        // // Allocate space for an Obj and assign
        // cudaMalloc((void**)&h_d_obj[i], sizeof(Obj));
        // // Copy the object to the device (only has single scalar field to keep it simple)
        // cudaMemcpy(h_d_obj[i], &(h_obj[i]), sizeof(Obj), cudaMemcpyHostToDevice);
    // }

    // /**************************************************/
    // /* CREATE DEVICE ARRAY TO PASS POINTERS TO KERNEL */
    // /**************************************************/

    // // Create a pointer which will point to device memory
    // Obj** d_d_obj = NULL;
    // // Allocate space for 3 pointers on device at above location
    // cudaMalloc((void**)&d_d_obj, sizeof(Obj*) * 3);
    // // Copy the pointers from the host memory to the device array
    // cudaMemcpy(d_d_obj, h_d_obj, sizeof(Obj*) * 3, cudaMemcpyHostToDevice);
