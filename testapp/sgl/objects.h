
#include <cstdlib>
#include <iostream>
#include <vector>
#include <stack>

class Color {
public:
	float red, green, blue;
};

class Vector4f{
public:
	float x, y, z, w;
	Vector4f(float x, float y, float z, float w){
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
};

class Matrix4f {
private:
	float m[4][4];



public:
	Matrix4f(){

		// TODO - matice identity
		/*m = (float*)calloc(16,sizeof(float));
		m[1][1] = 1;
		m[2][2] = 1;
		m[3][3] = 1;
		m[4][4] = 1;*/
	}

	Matrix4f(float *matrix){
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				m[j][i] = matrix[i * 4 + j];
			}
		}
	}

	Matrix4f mulByMatrix(Matrix4f right){
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
	}

	float* mulByVec(float *rightVector){
	}

};

class Context {
private:
	Color *currColor;
	Color *bcgColor;
	float *colorBuffer;
	std::vector<Vector4f> vertexBuffer;
	std::stack<Matrix4f> modelViewStack, projectionStack;
	Matrix4f *currMatrix;

public:
	// constructor
	Context(int width, int height){
		colorBuffer = (float*)calloc(width*height * 3, sizeof(float));
		vertexBuffer.reserve(8);
		//currMatrix = modelViewStack.top(); //????
		//vertexBuffer.push_back(new Vector4f(1., 2., 3.,4.)); //????
	}

	// destructor
	~Context(){
		free(colorBuffer);
		vertexBuffer.shrink_to_fit();
	}




};

