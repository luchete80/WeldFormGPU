
project (test LANGUAGES CXX CUDA)


set(CMAKE_CXX_STANDARD 11) #LIKE IN CUNSearch, this allows to include string (conflicts with CXX17)!!!
set(CMAKE_CXX_STANDARD_REQUIRED ON)



include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

add_executable(roundoff roundoff_test.cu
													#Domain.cpp NbSearch.cpp 
													#Vector.cpp
													) 
 ##
#add_subdirectory (cuda)
