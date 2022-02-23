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

#include "Functions.cuh"
#include "vector_math.h"
#include <stdio.h>

namespace SPH {

	inline double Kernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
	{
		double C;

		switch (KT)
		{
			case 0:	// Qubic Spline
				Dim == 2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);

				if 		(q<1.0)	return C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
				else if (q<2.0)	return C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q));
				else						return 0.0;
				break;

			case 1:	// Quintic
				Dim ==2 ? C = 7.0/(4.0*h*h*M_PI) : C = 7.0/(8.0*h*h*h*M_PI);

				if			(q<2.0)	return C*pow((1.0-q/2.0),4.0)*(2.0*q+1.0);
				else						return 0.0;
				break;

			case 2:	// Quintic Spline
				Dim ==2 ? C = 7.0/(478.0*h*h*M_PI) : C = 1.0/(120.0*h*h*h*M_PI);

				if			(q<1.0)	return C*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0)+15.0*pow((1.0-q),5.0));
				else if (q<2.0)	return C*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0));
				else if (q<3.0)	return C*(pow((3.0-q),5.0));
				else						return 0.0;
				break;

			case 3:	// Hyperbolic Spline, Fraser, Eq 3-27
				Dim ==2 ? C = 1.0/(3.0*h*h*M_PI) : C = 15.0/(62.0*h*h*h*M_PI);

				if		(q<1.0)	return C*( pow(q,3.0)- 6.0 *q  + 6. );
				else if (q<2.0)	return C*( pow((2.0-q),3.0) );
				else			return 0.0;
				break;
				
			default:
				//std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
				//std::cout << "0 => Qubic Spline" << std::endl;
				//std::cout << "1 => Quintic" << std::endl;
				//std::cout << "2 => Quintic Spline" << std::endl;
				//std::cout << "3 => Hyperbolic Spline" << std::endl;
				abort();
				break;
		}
	}

	inline double GradKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
	{
		double C;

		switch (KT)
		{
			case 0:	// Qubic Spline
				Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);
				if 		(q==0.0)	return C/h    *(-3.0+(9.0/2.0)*q);
				else if (q<1.0)		return C/(q*h)*(-3.0*q+(9.0/4.0)*q*q);
				else if (q<2.0)		return C/(q*h)*((-3.0/4.0)*(2.0-q)*(2.0-q));
				else							return 0.0;
				break;

			case 1:	// Quintic
				Dim ==2 ? C = 7.0/(4.0*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*M_PI);

				if 		(q==0.0)	return C*-5.0/h*(pow((1.0-q/2.0),3.0)-3.0*q/2.0*pow((1.0-q/2.0),2.0));
				else if (q<2.0)		return C/(q*h)*-5.0*q*pow((1.0-q/2.0),3.0);
				else							return 0.0;
				break;

			case 2:	// Quintic Spline
				Dim ==2 ? C = 7.0/(478.0*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*M_PI);

				if		(q==0.0)	return C/h*    (20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0)+300.0*pow((1.0-q),3.0));
				else if (q<1.0)		return C/(q*h)*(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0)-75.0*pow((1.0-q),4.0));
				else if (q<2.0)		return C/(q*h)*(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0));
				else if (q<3.0)		return C/(q*h)*(-5.0*pow((3.0-q),4.0));
				else							return 0.0;
				break;

			case 3:	// Hyperbolic Spline, Fraser, Eq 3-27 & Eqn 3-12, is dW/dR * 1/h *1/r
				Dim ==2 ? C = 1.0/(3.0*h*h*h*M_PI) : C = 15.0/(62.0*h*h*h*h*M_PI);

				if		(q<1.0)	return C/(q*h)*( 3.0 * pow(q,2.0)- 6.0 );
				else if (q<2.0)	return C/(q*h)*( 3.0 * pow((2.0-q),2.0) - 1.0 );
				else			return 0.0;
				break;
				
			default:
				//std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
				//std::cout << "0 => Qubic Spline" << std::endl;
				//std::cout << "1 => Quintic" << std::endl;
				//std::cout << "2 => Quintic Spline" << std::endl;
				//std::cout << "3 => Hyperbolic Spline" << std::endl;
				abort();
				break;
		}
	}

	inline double LaplaceKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
	{
		double C;

		switch (KT)
		{
			case 0:	// Qubic Spline
				Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

				if		(q<1.0)	return C*(-3.0+(9.0/2.0)*q)  + C*(Dim-1.0)/q * (-3.0*q+(9.0/4.0)*q*q);
				else if (q<2.0) return C*((3.0/2.0)*(2.0-q)) + C*(Dim-1.0)/q * ((-3.0/4.0)*(2.0-q)*(2.0-q));
				else						return 0.0;
				break;

			case 1:	// Quintic
				Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

				if 			(q<2.0)	return C*pow((1.0-q/2.0),2.0)*(10.0*q-5.0) + C*(Dim-1.0)/q * -5.0*q*pow((1.0-q/2.0),3.0);
				else						return 0.0;
				break;

			case 2:	// Quintic Spline
				Dim ==2 ? C = 7.0/(478.0*h*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*h*M_PI);

				if			(q<1.0)	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2-q),3.0)+300.0*pow((1-q),3.0))	+ C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0)-75.0*pow((1.0-q),4.0));
				else if (q<2.0)	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2-q),3.0))												+ C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0));
				else if (q<3.0)	return C*(20.0*pow((3.0-q),3.0))																						+ C*(Dim-1.0)/q * (-5.0*pow((3.0-q),4.0));
				else						return 0.0;
				break;

			default:
				//std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
				//std::cout << "0 => Qubic Spline" << std::endl;
				//std::cout << "1 => Quintic" << std::endl;
				//std::cout << "2 => Quintic Spline" << std::endl;
				abort();
				break;
		}
	}

	inline double SecDerivativeKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
	{
		double C;

		switch (KT)
		{
			case 0:	// Qubic Spline
				Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

				if 			(q<1.0)	return C*(-3.0+(9.0/2.0)*q);
				else if (q<2.0)	return C*((3.0/2.0)*(2.0-q));
				else						return 0.0;
				break;

			case 1:	// Quintic
				Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

				if			(q<2.0)	return C*pow((1.0-q/2.0),2.0)*(10.0*q-5.0);
				else						return 0.0;
				break;

			case 2:	// Quintic Spline
				Dim ==2 ? C = 7.0/(478.0*h*h*h*h*M_PI) : C = 1.0/(120.0*h*h*h*h*h*M_PI);

				if		(q<1.0)	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0)+300.0*pow((1.0-q),3.0));
				else if (q<2.0)	return C*(20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0));
				else if (q<3.0)	return C*(20.0*pow((3.0-q),3.0));
				else						return 0.0;
				break;

			default:
				//std::cout << "Kernel Type No is out of range. Please correct it and run again" << std::endl;
				//std::cout << "0 => Qubic Spline" << std::endl;
				//std::cout << "1 => Quintic" << std::endl;
				//std::cout << "2 => Quintic Spline" << std::endl;
				abort();
				break;
		}
	}

	__device__ __forceinline__ double EOS(size_t const & EQ, double const & Cs0, double const & P00, double const & Density, double const & Density0)
	{
		switch (EQ)
		{
			case 0:
				return P00+(Cs0*Cs0)*(Density-Density0);
				break;

			case 1:
				return P00+(Density0*Cs0*Cs0/7.0)*(pow(Density/Density0,7.0)-1);
				break;

			case 2:
				return (Cs0*Cs0)*Density;
				break;

			default:
				//std::cout << "Please correct Pressure Equation No and run again" << std::endl;
				//std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
				//std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
				//std::cout << "2 => (Cs*Cs)*Density" << std::endl;
				abort();
				break;
		}
	}

	inline double SoundSpeed(size_t const & EQ, double const & Cs0, double const & Density, double const & Density0)
	{
		switch (EQ)
		{
			case 0:
				return Cs0;
				break;

			case 1:
				return sqrt((Cs0*Cs0)*pow(Density/Density0,6.0));
				break;

			case 2:
				return Cs0;
				break;

			default:
				//std::cout << "Please correct Pressure Equation No and run again" << std::endl;
				//std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
				//std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
				//std::cout << "2 => (Cs*Cs)*Density" << std::endl;
				abort();
				break;
		}
	}

	inline double DensitySolid(size_t const & EQ, double const & Cs0, double const & P00, double const & Pressure, double const & Density0)
	{
		switch (EQ)
		{
			case 0:
				return (Pressure-P00)/(Cs0*Cs0) + Density0;
				break;

			case 1:
				return pow( ((Pressure-P00)*(7.0/(Density0*Cs0*Cs0))+1) , 1.0/7.0 ) * Density0;
				break;

			case 2:
				return Pressure/(Cs0*Cs0);
				break;

			default:
				//std::cout << "Please correct Pressure Equation No and run again" << std::endl;
				//std::cout << "0 => P0+(Cs*Cs)*(Density-Density0)" << std::endl;
				//std::cout << "1 => P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1)" << std::endl;
				//std::cout << "2 => (Cs*Cs)*Density" << std::endl;
				abort();
				break;
		}
	}

	// //NEW
	// iKernel::iKernel(size_t const & Dim,double const & h){
	
	// m_inv_h=1./h;
	
	// Dim ==2 ? m_w = 7.0/(478.0*h*h*M_PI) 		: m_w = 1.0/(120.0*h*h*h*M_PI);				//m_w
	// Dim ==2 ? m_gradw = 7.0/(478.0*h*h*h*M_PI) 	: m_gradw = 1.0/(120.0*h*h*h*h*M_PI);		//m_gradw
	// Dim ==2 ? m_lapw = 7.0/(478.0*h*h*h*h*M_PI) : m_lapw = 1.0/(120.0*h*h*h*h*h*M_PI);		//m_lapw
		
	// }
	
	// inline double iKernel::W( double const & q ) {
		// if		(q<1.0)	return m_w*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0)+15.0*pow((1.0-q),5.0));
		// else if (q<2.0)	return m_w*(pow((3.0-q),5.0)-6.0*pow((2.0-q),5.0));
		// else if (q<3.0)	return m_w*(pow((3.0-q),5.0));
		// else			return 0.0;
		
	// }
	
	// inline double iKernel::gradW( double const & q ) {
		// if		(q==0.0)	return m_gradw*m_inv_h*    (20.0*pow((3.0-q),3.0)-120.0*pow((2.0-q),3.0)+300.0*pow((1.0-q),3.0));
		// else if (q<1.0)		return m_gradw*m_inv_h/q *(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0)-75.0*pow((1.0-q),4.0));
		// else if (q<2.0)		return m_gradw*m_inv_h/q *(-5.0*pow((3.0-q),4.0)+30.0*pow((2.0-q),4.0));
		// else if (q<3.0)		return m_gradw*m_inv_h/q *(-5.0*pow((3.0-q),4.0));
		// else							return 0.0;		
		
		
	// }
	
	// THIS IS FROM GPUSPH, euler_kernel.cu

	/// Apply rotation to a given vector
	/*! Apply the rotation given by the matrix rot to the vector relPos.
	 *  The change in the relPos vector due to the rotation is computed
	 *  and added to the pos vector.
	 *
	 *	\param[in] rot : rotation matrix
	 *	\param[in] relPos: position with respect to center of gravity
	 *	\param[in] pos: position with respect to the local cell center
	 *
	 *	\return local postion rotated according to rot
	 */
	__device__ __forceinline__ void
	applyrot(const float* rot, const float3 & relPos, float4 & pos)
	{
		// Applying rotation
		pos.x += (rot[0] - 1.0f)*relPos.x + rot[1]*relPos.y + rot[2]*relPos.z;
		pos.y += rot[3]*relPos.x + (rot[4] - 1.0f)*relPos.y + rot[5]*relPos.z;
		pos.z += rot[6]*relPos.x + rot[7]*relPos.y + (rot[8] - 1.0f)*relPos.z;
	}

	/// Apply counter rotation to a given vector
	/*! Apply the inverse rotation given by the matrix rot to the vector relPos.
	 *  The change in the relPos vector due to the rotation is computed
	 *  and added to the pos vector.
	 *
	 *	\param[in] rot : rotation matrix
	 *	\param[in] relPos: position with respect to center of gravity
	 *	\param[in] pos: position with respect to the local cell center
	 *
	 *	\return local postion rotated according to rot^{-1}
	 */
	__device__ __forceinline__ void
	applycounterrot(const float* rot, const float3 & relPos, float4 & pos)
	{
		// Applying counter rotation (using R^{-1} = R^T)
		pos.x += (rot[0] - 1.0f)*relPos.x + rot[3]*relPos.y + rot[6]*relPos.z;
		pos.y += rot[1]*relPos.x + (rot[4] - 1.0f)*relPos.y + rot[7]*relPos.z;
		pos.z += rot[2]*relPos.x + rot[5]*relPos.y + (rot[8] - 1.0f)*relPos.z;
	}

	__device__ __forceinline__ void
	applyrot2(float* rot, float3 & pos, const float3 & cg)
	{
		float3 relpos = pos - cg;
		float3 new_relpos;

		// Applying rotation
		new_relpos.x = rot[0]*relpos.x + rot[1]*relpos.y + rot[2]*relpos.z;
		new_relpos.y = rot[3]*relpos.x + rot[4]*relpos.y + rot[5]*relpos.z;
		new_relpos.z = rot[6]*relpos.x + rot[7]*relpos.y + rot[8]*relpos.z;

		pos.x = new_relpos.x + cg.x;
		pos.y = new_relpos.y + cg.y;
		pos.z = new_relpos.z + cg.z;
	}

	__device__ inline void Rotation (float* Input, float* & Vectors, float* & VectorsT, float* & Values)
	{
	//	float3 Val,V0,V1,V2;
	//	Eig(Input,Val,V0,V1,V2,true,false);
	//
	//	float* V;
	//	Vectors	=	V0(0), V1(0), V2(0),
	//						V0(1), V1(1), V2(1),
	//						V0(2), V1(2), V2(2);

	//	Trans(Vectors,VectorsT);
	//	Values	= 	Val(0), 0.0   , 0.0   ,
	//							0.0   , Val(1), 0.0   ,
	//							0.0   , 0.0   , Val(2);
	}



}; // namespace SPH
