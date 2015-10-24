
#include <cstdlib>
#include <iostream>
#include <vector>
#include <stack>

class Color {
public:
	float red, green, blue;

	Color(float R, float G, float B){
		red = R;
		green = G;
		blue = B;
	}

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
	float **m;



public:
	Matrix4f(){

		// identity matrix
		m = (float**)calloc(4, sizeof(float*));
		for (int i = 0; i < 4; i++)
		{
			m[i] = (float*)calloc(4, sizeof(float));
			m[i][i] = 1;
		}
		
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

	float** getMatrix(){
		return m;
	}

	void toTerminal(){
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				std::cout << m[j][i]<<", ";
			}
			std::cout << std::endl;
		}
	}

};

class Context {
private:
	Color *currColor;
	Color *bcgColor;
	float *colorBuffer;
	std::vector<Vector4f> vertexBuffer;
	std::stack<Matrix4f> modelViewStack, projectionStack;
	std::stack<Matrix4f> *currMatrixStack;

public:
	// constructor
	Context(int width, int height){
		colorBuffer = (float*)calloc(width*height * 3, sizeof(float));
		// vertex buffer
		vertexBuffer.reserve(8);

		//matrix stacks
		modelViewStack.push(Matrix4f());
		projectionStack.push(Matrix4f());
		currMatrixStack = &modelViewStack;
		
		//colors
		bcgColor = new Color(0,0,0);
		currColor = new Color(1, 1, 1);

	}

	// destructor
	~Context(){
		free(colorBuffer);
		vertexBuffer.shrink_to_fit();
	}

	Matrix4f getMatrix(){
		return (*currMatrixStack).top();
	}

	Color* getDrawColor(){

	}

	float* getColorBuffer(){
		return colorBuffer;
	}


};

