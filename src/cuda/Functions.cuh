/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *
*             and soils) using Smoothed Particle Hydrodynamics method              *
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#ifndef SPH_SPECIAL_FUNCTIONS_H
#define SPH_SPECIAL_FUNCTIONS_H

//#include "matvec.h"
#include "tensor.cuh"

//enum Kernels_Type_d { Qubic_Spline=0, Quintic=1, Quintic_Spline=2 ,Hyperbolic_Spline=3};
	
namespace SPH {
	
	__device__ double Kernel								(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	__device__ inline double GradKernel						(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	__device__ double LaplaceKernel				(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	__device__ double SecDerivativeKernel	(size_t const & Dim, size_t const & KT, double const & q, double const & h);

	__device__ __forceinline__ double EOS									(size_t const & EQ, double const & Cs0, double const & P00, double const & Density, double const & Density0);

	__device__ __host__ __forceinline__ double SoundSpeed						(size_t const & EQ, double const & Cs0, double const & Density, double const & Density0);

	__device__ __host__ __forceinline__ double DensitySolid					(size_t const & EQ, double const & Cs0, double const & P00, double const & Pressure, double const & Density0);

	//__device__ void   Rotation							(tensor3 Input, tensor3 & Vectors, tensor3 & VectorsT, tensor3 & Values);

	//tensor3 abab									(tensor3 const & A, tensor3 const & B);
	
	
	//NEW
	// class iKernel{
		// public: //Ok TODO: Move members to private
		// iKernel(size_t const & Dim,double const & h);
		// iKernel(){}
		// virtual ~iKernel(){}
		// double gradW(double const & q);
		// double W(double const & q);
		
		// double m_w,m_gradw,m_lapw,m_inv_h;	//Precomputed Kernel constants

		
	// };

}; // namespace SPH

#include "Functions.cu"

#endif // SPH_SPECIAL_FUNCTIONS_H
