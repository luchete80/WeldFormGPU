/*  Copyright (c) 2011-2017 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  Vector.h
 *  NNS
 *
 *  Created by Alexis Herault on 27/07/06.
 *
 */

#ifndef _GEOMVECTOR_H
#define _GEOMVECTOR_H

#include "cuda_runtime.h"
#include <math.h>
//#include "vector_math.h"
//#include "Point.h"

//class Point;

/// 3D vector class
/*! 3D vector class provide :
	- standard operators for vectors
	- norm of vector
	- normal vector
	- access to coordinates values
*/
class Vector {
	private:
		double	x[4];	///< coordinates of vector

	public:
		//Vector(const Point &, const Point &);
		__host__ __device__ Vector(const Vector &);
		__host__ __device__ Vector(double xx = 0, double yy = 0, double zz = 0);
		__host__ __device__ Vector(const float3 &);
		__host__ __device__ Vector(const double3 &);
		__host__ __device__ Vector(const float4 &);
		__host__ __device__ Vector(const double4 &);
		__host__ __device__ Vector(const float *);
		__host__ __device__ Vector(const double *);
		__host__ __device__ ~Vector(void) {};

		/*! Return the norm of vector */
		__host__ __device__ double norm(void) const;
		void normalize(void);
		double normSquared(void) const;
		/*! Return a normal vector of vector */
		Vector Normal(void) const;
		/*! Return the vector rotated by the given angle (in radians)
			around the given vector.
		 */
		Vector rotated(const double &angle, const Vector &normal) const;
		Vector cross(const Vector &) const;

		/*! \name
			Overloaded operators
		*/
		//\{
		__host__ __device__ Vector &operator+=(const Vector &);
		__host__ __device__ Vector &operator-=(const Vector &);
		__host__ __device__ Vector &operator*=(double);
		__host__ __device__ Vector &operator/=(double);
		__host__ __device__ Vector &operator=(const Vector &);
		__host__ __device__ Vector &operator=(const float &);
		__host__ __device__ double &operator()(int);
		__host__ __device__ double operator()(int) const;
		//\}

		/*! \name
			Overloaded friend operators
		*/
		//\{
		__host__ __device__ friend Vector operator+(const Vector &, const Vector &);
		__host__ __device__ friend Vector operator-(const Vector &, const Vector &);
		__host__ __device__ friend Vector operator*(double, const Vector &);
		__host__ __device__ friend Vector operator*(const Vector &, double);
		__host__ __device__ friend double operator*(const Vector &, const Vector &);
		__host__ __device__ friend Vector operator/(const Vector &, double);
		__host__ __device__ friend Vector operator-(const Vector &);
		//\}

		// DEBUG
		void print(void);
};

inline double dot (const Vector &v1, const Vector &v2){
  return (v1(0)*v2(0)+v1(1)*v2(1)*+v1(2)*v2(2));
}

inline Vector cross(const Vector &a, const Vector &b){
	Vector ret;
	ret(0) = a(1)*b(2)-a(2)*b(1);
	ret(1) = a(2)*b(0)-a(0)*b(2);
	ret(2) = a(0)*b(1)-a(1)*b(0);
	return ret;
}

float3 make_float3(const Vector &);
double3 make_double3(const Vector &);
float4 make_float4(const Vector &);
double4 make_double4(const Vector &);

#endif
