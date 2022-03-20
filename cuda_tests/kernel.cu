//https://stackoverflow.com/questions/65325842/how-do-i-properly-implement-classes-whose-members-are-called-both-from-host-and

#include <iostream>
#include <cuda_runtime.h>

#define CUDA_CHECK cuda_check(__FILE__,__LINE__)
inline void cuda_check(std::string file, int line)
{
    cudaError_t e = cudaGetLastError();
    if (e != cudaSuccess) {
        std::cout << std::endl
                  << file << ", line " << line << ": "
                  << cudaGetErrorString(e) << " (" << e << ")" << std::endl;
        exit(1);
    }
}

class TestTable {

    float* vector_;
    int num_cells_;

public:

    void Init() {
        num_cells_ = 1e4;
        cudaMallocManaged(&vector_, num_cells_*sizeof(float));
        CUDA_CHECK;
    }

    void Free() {
        cudaFree(vector_);
    }

    __device__
    bool UpdateValue(int global_index, float val) {
        int index = global_index % num_cells_;
        vector_[index] = val;
        return false;
    }

};

class TestClass {

private:

    float value_;
    TestTable* test_table_;

public:

    TestClass() : value_(1.) {
        // test_table_ = new TestTable;
        cudaMallocManaged(&test_table_, sizeof(TestTable));
        test_table_->Init();
        CUDA_CHECK;
    }

    ~TestClass() {
        test_table_->Free();
        cudaFree(test_table_);
        CUDA_CHECK;
    }

    __host__ __device__
    float GetValue() {
        return value_;
    }

    __host__
    void RunKernel();

};

__global__
void test_kernel(TestClass* test_class, TestTable* test_table) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < 1e6; i += stride) {
        const float val = test_class->GetValue();
        test_table->UpdateValue(i, val);
    }
}
  ////
__host__
void TestClass::RunKernel() {
    test_kernel<<<1,1>>>(this, test_table_);
    cudaDeviceSynchronize(); CUDA_CHECK;
}

int main(int argc, char *argv[]) {

    //TestClass* test_class = new TestClass();
    TestClass* test_class;
		cudaMallocManaged(&test_class, sizeof(TestClass));
		new(test_class) TestClass();

		std::cout << "TestClass successfully constructed" << std::endl;

    test_class->RunKernel();
    std::cout << "Kernel successfully run" << std::endl;

    delete test_class;
    std::cout << "TestClass successfully destroyed" << std::endl;
		
		test_class->~TestClass();
		cudaFree(test_class);

    return 0;
}