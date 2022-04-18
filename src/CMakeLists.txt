SET(BUILD_GPU ON)


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


set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --fmad=true and --fmad=false enables and disables the contraction respectively. 
 # This switch is supported only when the --gpu-architecture option is set with     
 # compute_20, sm_20, or higher. For other architecture classes, the contraction is     
  # always enabled. 
 # The --use_fast_math option implies --fmad=true, and enables the contraction.
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fast_math=false")

SET(EXES
	#Ale
	#Ale_modelo6
	#SolidBlock_Heat_conv3faces_Flux_n30
	#Ale_Modelo6_placa
	# 1-Poiseuille
	# 2-Couette
	# 3-DamBreak
	# 4-DamBreak3D
	#Cylinder_Test
	#Plate-Yield
	#Compression
	#Compression_10
	#Compression_12_solver_orig
	#Compression_smallD
	#Compression_smallL
	#Compression_smallL_1side
	#Compression_smallD_CNS
	#Compression_10ts
	#Compression_Elastic_CNS_10ts
	#Compression_Elastic_CNS
	#Cantilever
	#Cantilever2
	#Cantilever2_pl
	#Elastic
	main
	#Heat_2D
	# Traction
	# 4-DamBreak3D
	#Plate-Yield_hcorr
	#Plate-Yield_CNS
	#Plate-Yield_Hard
	#Plate-Yield_CNS_Hard
	#Plate-Yield_CNS_v2
	#Plate-Yield_3d_al
	#Plate-Yield_3d_al_CNS
	#SolidBlock_Heat
	#SolidBlock_Heat_CNS
	#SolidBlock_Heat_conv3faces
	#SolidBlock_Heat_source
	#SolidBlock_Heat_conv3faces_Flux
	#SolidBlock_Heat_source2
	#SolidBlock_Heat_source2_boundpart
	#Neighbour_Test
	#Neighbour_Test_CNS
	#"Source/WeldForm"
)

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

add_executable(gpuweldform 
													#main_thermal.cu
													main_mech.cu
													Domain.cpp 
													#NbSearch.cpp 
													Output.cpp
													Vector.cpp
													#TODO: ALLOW BUILD CUDA DIR WITH SEPARABLE COMPILATION
													cuda/Domain_d.cu
													#cuda/tensor.cu
													#cuda/Interaction.cu
													#cuda/Mechanical.cu
													cuda/Thermal.cu
													#cuda/Interaction_global.cu	#TEMPORARY, ALL IS EXTERNAL
													cuNSearch/src/PointSet.cu
													cuNSearch/src/PointSetImplementation.cu
													cuNSearch/src/cuNSearch.cu
													cuNSearch/src/cuNSearchDeviceData.cu
													cuNSearch/src/cuNSearchKernels.cu
													cuNSearch/Utils/cuda_helper.cu
													cuNSearch/Utils/Timing.cpp
													cuNSearch/Utils/IDFactory.cpp													
													) 

# add_executable(test_cuns 
													# main_mech_cuns.cu
													# Domain.cpp 
													# Output.cpp
													# Vector.cpp
													# cuda/Domain_d.cu
													# cuda/Thermal.cu
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

#TARGET_LINK_LIBRARIES (gpuweldform cuda)

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
