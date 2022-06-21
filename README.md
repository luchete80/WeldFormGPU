WeldFormGPU is the GPU version of [WeldForm](https://github.com/luchete80/WeldForm)

Not tested on Linux yet.

## Features
Has been exclusively adapted to solid mechaincs, and it includes:

- Mechanic Solver
- Thermal Solver
- Coupled ThermoMechanical Solver (in progress)
- Contact formulation (in progress)
- Adaptive search only in case of plastic strain threshold (in progress)

![alt text](https://github.com/luchete80/WeldForm/blob/master/compression.PNG)

Is hevaily based on: 

1) [PersianSPH](https://github.com/mghkorzani/persiansph) - Maziar Gholami Korzani and Sergio Galindo Torres
2) Kirk Fraser [Thesis](https://constellation.uqac.ca/4246/1/Fraser_uqac_0862D_10345.pdf) and [works](https://pdfs.semanticscholar.org/b09e/8c8023d56b130cc6fa5314cb66bce364df8e.pdf) on SPH model of FSW


## Building Instructions

1) Install Visual Studio Community [2019](https://visualstudio.microsoft.com/es/vs/older-downloads/) (Not tested on 2022 yet) 
2) Install [cmake](https://cmake.org/download/)
3) Download and install CUDA [compilers](https://developer.nvidia.com/cuda-downloads)
4) Set CUDA compiler path on CMD: set CUDA_PATH="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4"
5) Create a direcory for built binaries
6) Run make.bat (located here on root directory), to set MSVC env vars
7) Clone this repo (for example to c:\WeldFormGPU\src)
8) Run cmake c:\WeldFormGPU\src -G "NMake Makefiles"
9) Run nmake.exe 

