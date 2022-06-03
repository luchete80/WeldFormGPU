WeldFormGPU is the GPU version of [WeldForm](https://github.com/luchete80/WeldForm)

Work in progress.

- replaces Vec3_t to Vector for CPU work
- implements (via vector_math.h) 

tensor is defined on CUDA, 


https://docs.nvidia.com/cuda/pdf/CUDA_C_Best_Practices_Guide.pdf
https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html


set CUDA_PATH = C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4

"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"

cmake ..\WeldFormGPU -G "NMake Makefiles"


