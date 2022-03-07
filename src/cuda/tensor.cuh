/*  Copyright (c) 2013-2015 INGV, EDF, UniCT, JHU

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
#ifndef _TENSOR_H_
#define _TENSOR_H_

typedef struct {
	double xx;
	double xy;
	double xz;
	double yy;
	double yz;
	double zz;
} symtensor3 ;

typedef struct {
	double xx;
	double xy;
	double xz;
	double xw;
	double yy;
	double yz;
	double yw;
	double zz;
	double zw;
	double ww;
} symtensor4 ;

#define __spec __device__ __forceinline__

class tensor3
{
private:
		double m_data[4][4]{};    
public:
    __device__ tensor3();
    __device__ tensor3(double flat[]);	//Six components
    __device__ void FromFlatSym(double flat[]);	//Six components
		__device__ void ToFlatSymPtr(double *flat, int initial);
    __device__ void FromFlatAntiSym(double flat[]);	//Six components
    __device__ void FromFlatSymPtr(double *flat);	//Six components
    __device__ double& operator()(int row, int col);
    __device__ double& operator[](const int &i);	//0 to 8
    __device__ double operator()(int row, int col) const;
    __device__ void operator()();
		__device__ void Set(const int &i, const int &j);
		__device__ tensor3 operator+ (const tensor3 &b);
		__device__ tensor3 operator- (const tensor3 &b);
		__device__ tensor3 operator* (const tensor3 &b);
		__device__ tensor3 operator- (const double &f);
		__device__ tensor3 operator= (const double &f);
		
		__spec tensor3 Identity();
		//tensor3 operator* (const double &f);
		__device__ tensor3 Trans ();
		__device__ ~tensor3(){};
};

__device__ tensor3 operator* (const double &f, const tensor3 &b);
__device__ tensor3 operator/ (const tensor3 &b, const double &f);

__device__ double3 dot(tensor3 const& T, double3 const& v);

__device__ /*__forceinline__*/ tensor3 Identity();

#endif
