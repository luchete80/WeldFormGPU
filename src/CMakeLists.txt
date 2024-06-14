SET(BUILD_GPU ON)

if (BUILD_GPU)
project (gpuweldform LANGUAGES CXX CUDA)
project (testcuns LANGUAGES CXX CUDA)
include_directories(${PROJECT_SOURCE_DIR})

include_directories(cuNSearch/include)
include_directories(cuNSearch/src)
include_directories(cuNSearch/Utils)

link_directories(${CUNSEARCH_LIB_DIR})
link_directories(${cuNSearch_BINARY_DIR})


# target_link_libraries (gpuweldform cuNSearch )

#set (CMAKE_CXX_FLAGS "-fopenmp")


set(CMAKE_CXX_STANDARD 14) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --fmad=true and --fmad=false enables and disables the contraction respectively. 
 # This switch is supported only when the --gpu-architecture option is set with     
 # compute_20, sm_20, or higher. For other architecture classes, the contraction is     
  # always enabled. 
 # The --use_fast_math option implies --fmad=true, and enables the contraction.
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false --prec-sqrt=true --ftz=false --prec-div=true")
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math=false")


# FROM --use_fast_math implies --ftz=true --prec-div=false --prec-sqrt=false --fmad=true.
# https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html


set_source_files_properties(Domain_d.C Matrices.C Mechanical.C Solver.C Tensor.C main.C  main_1_elem_3d.C
                            WeldFormFEM.C include/*.h 
                            PROPERTIES LANGUAGE CUDA)
                            
                            
else ()

endif()


SET(EXES

	main

)

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

##------------THIS WORKS OK -------------------
# add_executable(weldform 
                          # WeldForm.cpp
													# Domain.cpp 
													# Output.cpp
													# Vector.cpp
                          # Mesh.cpp
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
													# cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													                          
                          # ) 
# add_executable(main_mech 
													# main_mech.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# Output.cpp
													# Vector.cpp
                          # Mesh.cpp
                          # cuda/thermal.cu
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
                          # cuda/Geom.cu
													# cuda/Domain_d.cu
                          # cuda/cudautils.cu
													# #cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# )                          
                 

# add_executable(main_contact 
													# main_contact.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# NastranReader.cpp
													# Output.cpp
													# Vector.cpp
                          # Mesh.cpp
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
                          # cuda/cudautils.cu
													# cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# )   

add_executable(WeldFormGPU
                        WeldFormGPU.cu
                        Domain.cpp 
													Output.cpp
													Vector.cpp
                          Mesh.cpp
                          NastranReader.cpp
													cuda/Domain_d.cu
                          cuda/cudautils.cu
                          cuda/Thermal.cu
                          cuda/Mesh_reader.cu #TO AVOID CONFLICT WITH INLINE FUNCTIONS
                          
													cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													cuNSearch/Utils/IDFactory.cpp			
              )

# add_executable(main_thermal 
													# main_thermal.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# Output.cpp
													# Vector.cpp
                          # cuda/thermal.cu
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
													# #cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu cuNSearch/src/cuNSearchDeviceData.cu cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# )
                          

# add_executable(main_pairs 
													# #main_thermal.cu
													# main_thermal_pairs.cu
													# Domain.cpp 
													# #NbSearch.cpp 
													# Output.cpp
													# Vector.cpp
													# #TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													# cuda/Domain_d.cu
													# cuda/Thermal.cu
													# #cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													# cuNSearch/src/PointSet.cu
													# cuNSearch/src/PointSetImplementation.cu
													# cuNSearch/src/cuNSearch.cu
													# cuNSearch/src/cuNSearchDeviceData.cu
													# cuNSearch/src/cuNSearchKernels.cu
													# cuNSearch/Utils/cuda_helper.cu
													# cuNSearch/Utils/Timing.cpp
													# cuNSearch/Utils/IDFactory.cpp													
													# ) 


#tensor and all functions defined in separate modulus
# CASE NO INLINE VERSION
#set_property(TARGET gpuweldform PROPERTY CUDA_SEPARABLE_COMPILATION ON)

#add_subdirectory (cuda)

TARGET_LINK_LIBRARIES (WeldFormGPU cuda  lsdynareader )

set_property(TARGET WeldFormGPU PROPERTY CUDA_ARCHITECTURES 75)

#FOREACH(var ${EXES})
#    ADD_EXECUTABLE        (${var} "${var}.cpp" #
#
#	#"${CMAKE_SOURCE_DIR}/Source/CompactNSearch/CompactNSearch.cpp"
#	)
#    TARGET_LINK_LIBRARIES (${var} ${LIBS} )
#    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
#ENDFOREACH(var)

# if (NOT BUILD_GPU)
# else()
# TARGET_LINK_LIBRARIES (${var} cuda)
# endif()
    # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")
    # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)
    # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
