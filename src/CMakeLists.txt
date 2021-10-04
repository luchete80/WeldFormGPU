SET(BUILD_GPU ON)

if (NOT BUILD_GPU)

else ()

project (gpuweldform LANGUAGES CXX CUDA)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(Source_cuda/CuNSearch/include)
include_directories(${CUNSEARCH_SRC_DIR}/Utils)
link_directories(${CUNSEARCH_LIB_DIR})
# add_executable(gpuweldform example_gpu.cu) 
# target_link_libraries (gpuweldform cuNSearch )

endif()

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

add_subdirectory (cuda)

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp" 

	#"${CMAKE_SOURCE_DIR}/Source/CompactNSearch/CompactNSearch.cpp"
	)
    TARGET_LINK_LIBRARIES (${var} ${LIBS} )
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

if (NOT BUILD_GPU)
else()
TARGET_LINK_LIBRARIES (${var} cuda)
endif()
    # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")
    # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)
    # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
