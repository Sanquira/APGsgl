#ifndef OBJECTS_H_

#define OBJECTS_H_

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stack>
#include <sgl.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <climits>

using namespace std;

class Color {
public:
	float red=0, green=0, blue=0;

	Color(){}

	Color(float R, float G, float B){
		red = R;
		green = G;
		blue = B;
	}

	bool compare(Color clr){
		if (clr.red == red&&clr.green == green&&clr.blue == blue){
			return true;
		}
		return false;
	}
	
	void toTerminal(){
		cout << red << ", " << green << ", " << blue << endl;
	}

};

class Vector4f{
public:
	float x = 0;
	float y = 0;
	float z = 0;
	float w = 1;

	Vector4f(){
	}

	Vector4f(float x, float y, float z, float w){
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	Vector4f minus(Vector4f &vector){
		Vector4f ret;
		ret.x = x - vector.x;
		ret.y = y - vector.y;
		ret.z = z - vector.z;
		ret.w = w;
		return ret;
	}
	
	float dotNoHomo(Vector4f &vector){
		float ret=0;
		ret += x*vector.x;
		ret += y*vector.y;
		ret += z*vector.z;
		return ret;
	}
	
	Vector4f cross(Vector4f &vector){
		Vector4f ret;
		ret.x = y*vector.z-z*vector.y;
		ret.y = z*vector.x-x*vector.z;
		ret.z = x*vector.y-y*vector.x;
		return ret;
	}
	
	Vector4f mulByConst(float num){
		Vector4f ret;
		ret.x = x * num;
		ret.y = y * num;
		ret.z = z * num;
		return ret;
	}

	void homoNorm(){
		x = x / w;
		y = y / w;
		z = z / w;
		w = w / w;
	}
	
	void normalize(){
		float norm = sqrt(x*x+y*y+z*z);
		x = x/norm;
		y = y/norm;
		z = z/norm;
	}
	
	void removeHomo(){
		w = 0;
	}

	void toTerminal(){
		cout << x << ", " << y << ", " << z << ", " << w << endl;
	}
};

class Matrix4f {
public:
	float m[4][4];
	Matrix4f(){

		// identity matrix
		for (int i = 0; i < 4; i++)
		{
			for (int j = i + 1; j < 4; j++)
			{
				m[j][i] = 0;
				m[i][j] = 0;
			}
			m[i][i] = 1;
		}

	}

	Matrix4f(const float *matrix){
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				m[j][i] = matrix[i * 4 + j];
			}
		}
	}

	Matrix4f mulByMatrix(Matrix4f &right){
		Matrix4f ret;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret.m[i][j] = 0;
				for (int k = 0; k < 4; k++)
				{
					ret.m[i][j] += m[i][k] * right.m[k][j];
				}
			}
		}
		return ret;
	}

	void mulByMatrixToItself(Matrix4f *right){
		Matrix4f tmp = mulByMatrix(*right);
		memcpy(m, tmp.m, sizeof(m));
	}

	Vector4f mulByVec(Vector4f &vec){
		float temp1[4] = { vec.x, vec.y, vec.z, vec.w };
		float temp2[4] = {};
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				temp2[i] += temp1[j] * m[i][j];
			}
		}
		Vector4f ret;
		ret.x = temp2[0];
		ret.y = temp2[1];
		ret.z = temp2[2];
		ret.w = temp2[3];
		return ret;
	}

	Matrix4f inverse(){
		float x[4][4];
		memcpy(x,m,sizeof(m));
		int indxc[4],indxr[4],ipiv[4];
		int i,icol,irow,j,k,l,ll,n;
		float big,dum,pivinv,temp;
		// satisfy the compiler
		icol = irow = 0;
		
		// the size of the matrix
		n = 4;
		
		for ( j = 0 ; j < n ; j++){ /* zero pivots */
			ipiv[j] = 0;
		}
		
		for ( i = 0; i < n; i++){
			big = 0.0;
			for (j = 0 ; j < n ; j++){
				if (ipiv[j] != 1){
					for ( k = 0 ; k<n ; k++){
						if (ipiv[k] == 0){
							if (fabs(x[k][j]) >= big){
								big = fabs(x[k][j]);
								irow = j;
								icol = k;
							}
						}else if (ipiv[k] > 1){
							cerr << "ERROR MATRIX IS SINGULAR!!!" << endl;
//							return 1; /* singular matrix */
						}
					}
				}
			}
			++(ipiv[icol]);
			if (irow != icol){
				for ( l = 0 ; l<n ; l++){
					temp = x[l][icol];
					x[l][icol] = x[l][irow];
					x[l][irow] = temp;
				}
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (x[icol][icol] == 0.0){
				cerr << "ERROR MATRIX IS SINGULAR!!!" << endl;
//				return 1; /* singular matrix */
			}
				
			pivinv = 1.0 / x[icol][icol];
			x[icol][icol] = 1.0 ;
			for ( l = 0 ; l<n ; l++){
				x[l][icol] = x[l][icol] * pivinv ;
			}
				
			for (ll = 0 ; ll < n ; ll++){
				if (ll != icol){
					dum = x[icol][ll];
					x[icol][ll] = 0.0;
					for ( l = 0 ; l<n ; l++){
						x[l][ll] = x[l][ll] - x[l][icol] * dum ;
					}
				}
			}
		}
		for ( l = n; l--; ){
			if (indxr[l] != indxc[l]){
				for ( k = 0; k<n ; k++){
					temp = x[indxr[l]][k];
					x[indxr[l]][k] = x[indxc[l]][k];
					x[indxc[l]][k] = temp;
				}
			}
	  	}
	Matrix4f ret;
	memcpy(ret.m,x,sizeof(x));
	return ret;
	}

	void toTerminal(){
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				std::cout << m[i][j] << ", ";
			}
			std::cout << std::endl;
		}
	}

};

class Material {
public:

	Color color = Color(1,1,1);
	float kd = 0;
	float ks = 0;
	float shine = 0;
	float T = 0;
	float ior = 0;

	Material(){
	}

	Material(float r, float g, float b, float kd, float ks, float shine, float T, float ior){
		this->color = Color(r,g,b);
		this->kd = kd;
		this->ks = ks;
		this->shine = shine;
		this->T = T;
		this->ior = ior;
	}
};

bool compareVector4fX(Vector4f i, Vector4f j) {
	return (i.x < j.x);
}


#endif

