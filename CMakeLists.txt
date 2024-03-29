#####################################################################################
# PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids         #
#             and soils) using Smoothed Particle Hydrodynamics method               #
# Copyright (C) 2016 Maziar Gholami Korzani and Sergio Galindo-Torres               #
#                                                                                   #
# This file is part of PersianSPH                                                   #
#                                                                                   #
# This is free software; you can redistribute it and/or modify it under the         #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 3 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY   #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# PersianSPH; if not, see <http://www.gnu.org/licenses/>                            #
#####################################################################################

CMAKE_MINIMUM_REQUIRED (VERSION 2.8.12)

PROJECT                (New_Test)
cmake_policy(SET CMP0054 NEW) #elseif
add_definitions(-D_CRT_NO_VA_START_VALIDATION)

add_definitions(-D_USE_MATH_DEFINES) #COMPACTNSEARCH??

include_directories("lib/cuNSearch/include")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
MESSAGE("Using GNU compiler")
include_directories("D:/Luciano/Numerico/blitz/include")
link_directories("D:/Luciano/Numerico/blitz/lib")
  # using GCC
link_directories("D:/Luciano/Numerico/CompactNSearch_bin_mingw")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
MESSAGE("Using MSVC compilers")
# using Visual Studio C++
include_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/gsl-2.5_install_mswin/include")
include_directories("D:/Luciano/Numerico/Libs/hdf5-1.12.0_install_mswin/include")
link_directories("D:/Luciano/Numerico/Libs/blitz-0.9_install_mswin/lib")
link_directories("D:/Luciano/Numerico/CompactNSearch_bin_nmake_release")

include_directories("src")
include_directories("src/LSDynaReader")
include_directories("External")

endif()

include_directories("LSDynaReader/src")
#BEFORE WELDFORM
add_subdirectory("LSDynaReader/src") 
add_subdirectory("src")
add_subdirectory("cuda_tests")




