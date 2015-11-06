
#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <vector>
#include <stack>
#include <sgl.h>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

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
	float vec[4];

	Vector4f(){
		Vector4f(0, 0, 0, 1);
	}

	Vector4f(float x, float y, float z, float w){
		vec[0] = x;
		vec[1] = y;
		vec[2] = z;
		vec[3] = w;
	}

	void toTerminal(){
		for (size_t i = 0; i < 4; i++)
		{
			std::cout << vec[i] << ", ";
		}
		std::cout << std::endl;
	}
};

class Matrix4f {
public:
	float m[4][4];
	Matrix4f(){

		// identity matrix
		for (int i = 0; i < 4; i++)
		{
			for (int j = i+1; j < 4; j++)
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

	Matrix4f* mulByMatrix(Matrix4f* right){
		Matrix4f* ret = new Matrix4f();
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				ret->m[i][j] = 0;
				for (int k = 0; k < 4; k++)
				{
					ret->m[i][j] += m[i][k] * (*right).m[k][j];
				}
			}
		}
		return ret;
	}

	void mulByMatrixToItself(Matrix4f * right){
		Matrix4f* tmp = mulByMatrix(right);
		memcpy(m, tmp->m, sizeof(m));
		delete tmp;
	}

	Vector4f* mulByVec(Vector4f * vec){
		Vector4f* ret = new Vector4f();
		for (int i = 0; i < 4; i++)
		{
			ret->vec[i] = 0;
			for (int j = 0; j < 4; j++)
			{
				ret->vec[i] += vec->vec[j] * m[i][j];
			}
		}
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

class Context {
private:
	Color *currColor;
	Color *bcgColor;
	float *colorBuffer;
	int width, height;
	std::vector<Vector4f> vertexBuffer;
	std::stack<Matrix4f> modelViewStack, projectionStack;
	std::stack<Matrix4f> *currMatrixStack;
	int areaDrawMode, vertexDrawMode;
	float pointSize;
	bool depthTestEnable;
	Matrix4f viewport;

	void symetricPoints(int x, int y, Vector4f * center){
		int xs = (int)center->vec[0];
		int ys = (int)center->vec[1];
		setPixel(y + xs, x + ys, currColor);
		setPixel(x + xs, y + ys, currColor);
		setPixel(-x + xs, y + ys, currColor);
		setPixel(-y + xs, x + ys, currColor);
		setPixel(-y + xs, -x + ys, currColor);
		setPixel(-x + xs, -y + ys, currColor);
		setPixel(x + xs, -y + ys, currColor);
		setPixel(y + xs, -x + ys, currColor);
	}

public:
	// constructor
	Context(int width, int height){
		this->width = width;
		this->height = height;
		colorBuffer = (float*)calloc(width*height * 3, sizeof(float));
		// vertex buffer
		vertexBuffer.reserve(8);

		//matrix stacks
		modelViewStack.push(Matrix4f());
		projectionStack.push(Matrix4f());
		currMatrixStack = &modelViewStack;

		//colors
		bcgColor = new Color(0, 0, 0);
		currColor = new Color(1, 1, 1);

		//draw modes
		areaDrawMode = 2; //SGL_FILL

	}

	// destructor
	~Context(){
		//buffers
		vertexBuffer.clear();
		free(colorBuffer);
		while (!modelViewStack.empty()){
			modelViewStack.pop();
		}
		delete(&modelViewStack);
		while (!projectionStack.empty()){
			projectionStack.pop();
		}
		delete(&projectionStack);

		//others
		delete(&viewport);
		delete(bcgColor);
		delete(currColor);
	}

	void setMatrixMode(int mode){
		if (mode == 0){
			currMatrixStack = &modelViewStack;
		}
		if (mode == 1){
			currMatrixStack = &projectionStack;
		}
	}

	void pushMatrix(){
		Matrix4f * mat = new Matrix4f();
		mat->mulByMatrixToItself(&currMatrixStack->top());
		currMatrixStack->push(*mat);
	}

	void popMatrix(){
		currMatrixStack->pop();
	}

	std::stack<Matrix4f>* getMatrixStack(){
		return currMatrixStack;
	}

	void setDrawColor(float r, float g, float b){
		delete(currColor);
		currColor = new Color(r, g, b);
	}

	Color* getBcgColor(){
		return bcgColor;
	}

	void setBcgColor(float r, float g, float b){
		delete(bcgColor);
		bcgColor = new Color(r, g, b);
	}

	void clearVertexBuffer(){
		vertexBuffer.clear();
	}

	std::vector<Vector4f> getVertexBuffer(){
		return vertexBuffer;
	}

	void addToVertexBuffer(Vector4f vec){
		vertexBuffer.push_back(vec);
	}

	float* getColorBuffer(){
		return colorBuffer;
	}

	int getWidth(){
		return width;
	}

	int getHeight(){
		return height;
	}

	void setAreaDrawMode(int mode){
		areaDrawMode = mode;
	}

	int getAreaDrawMode(){
		return areaDrawMode;
	}

	void setVertexDrawMode(int mode){
		vertexDrawMode = mode;
	}

	int getVertexDrawMode(){
		return vertexDrawMode;
	}

	void setPointSize(float size){
		pointSize = size;
	}

	float getPointSize(){
		return pointSize;
	}

	void setDepthTest(bool ena){
		depthTestEnable = ena;
	}

	bool getDepthTest(){
		return depthTestEnable;
	}

	Matrix4f getViewport(){
		return viewport;
	}

	void setViewport(Matrix4f mat){
		memcpy(&viewport, &mat, sizeof(Matrix4f));
	}

	Matrix4f* computeTransformation(){
	//	projectionStack.top().toTerminal();
	//	modelViewStack.top().toTerminal();
	//	viewport.toTerminal();
	//	cout << "=====" << endl;
		Matrix4f * tmp = projectionStack.top().mulByMatrix(&modelViewStack.top());
		return viewport.mulByMatrix(tmp);
	}

	//------------------------------------------------------------------
	//RENDERING METHODS
	//------------------------------------------------------------------

	void setPixel(int x, int y, Color * clr){
		if (x < 0 || y < 0 || x >= width || y >= height)
			return;
		colorBuffer[((y)*width + x) * 3 + 0] = clr->red;
		colorBuffer[((y)*width + x) * 3 + 1] = clr->green;
		colorBuffer[((y)*width + x) * 3 + 2] = clr->blue;
	}

	void renderCircle(Vector4f * vec, float radii){
		Matrix4f* mat = computeTransformation();
	//	mat->toTerminal();
	//	cout << "*******************" << endl;
		float scale = std::sqrt(mat->m[0][0] * mat->m[1][1] - mat->m[1][0] * mat->m[0][1]);
		Vector4f* center = mat->mulByVec(vec);
		int x, y, p;
		x = 0;
		y = (int)(radii*scale);
		p = (int)(1 - radii*scale); // p pro [0,r]
		while (x < y) {
			symetricPoints(x, y, center);
			if (p < 0) {
				p += 2 * x + 1;
			}
			else {
				p += 2 * (x - y) + 1;
				y -= 1;
			}
			x += 1;
		}
		if (x == y) {// 45° pixely, jen 4
			symetricPoints(x, y, center);
		}
		delete(mat);
		delete(center);
	}

	void renderEllipse(Vector4f * center, float a, float b){
		Matrix4f* mat = computeTransformation();

		for (int i = 0; i < 40; i++)
		{
			float uhel = (float)(M_PI*i / 20.);
			float uhelDalsi = (float)(M_PI*(i + 1) / 20.);
			Vector4f *p1 = new Vector4f(center->vec[0] + a*cos(uhel), center->vec[1] + b*sin(uhel), center->vec[2], 1);
			Vector4f *p2 = new Vector4f(center->vec[0] + a*cos(uhelDalsi), center->vec[1] + b*sin(uhelDalsi), center->vec[2], 1);
			renderLine(mat->mulByVec(p1), mat->mulByVec(p2));
			delete(p1);
			delete(p2);
		}
		delete(mat);
	}

	void renderArc(Vector4f * center, float radius, float from, float to){
		Matrix4f* mat = computeTransformation();

		int numOfVert = (int)(40 * fabs(to - from) / (2 * M_PI));
		float step = (to - from) / (numOfVert - 1);

		for (int i = 0; i < numOfVert - 1; i++){
			Vector4f *p1 = new Vector4f(center->vec[0] + cos(from + i*step)*radius, center->vec[1] + sin(from + i*step)*radius, center->vec[2], 1);
			Vector4f *p2 = new Vector4f(center->vec[0] + cos(from + (i + 1)*step)*radius, center->vec[1] + sin(from + (i + 1)*step)*radius, center->vec[2], 1);
			renderLine(mat->mulByVec(p1), mat->mulByVec(p2));
			delete(p1);
			delete(p2);
		}
		delete(mat);
	}

	void renderLine(Vector4f * p1, Vector4f * p2){
		int x1 = (int)round(p1->vec[0]);
		int y1 = (int)round(p1->vec[1]);
		int x2 = (int)round(p2->vec[0]);
		int y2 = (int)round(p2->vec[1]);

		bool uhel = (abs(y2 - y1) > abs(x2 - x1)); //>45°
		if (uhel){
			std::swap(x1, y1);
			std::swap(x2, y2);
		}

		if (x1 > x2){ //leva polorovina
			std::swap(x1, x2);
			std::swap(y1, y2);
		}

		int dx = x2 - x1;
		int dy = abs(y2 - y1);

		float error = dx / 2.0f;
		int ystep = (y1 < y2) ? 1 : -1;
		int y = y1;

		int maxX = x2;

		for (int x = x1; x <= maxX; x++){
			if (uhel){
				setPixel(y, x, currColor);
			}
			else{
				setPixel(x, y, currColor);
			}
			error -= dy;
			if (error < 0){
				y += ystep;
				error += dx;
			}
		}
	}

	void renderPoint(Vector4f * vec){
		int tmp = (int)(pointSize / 2);
		for (int i = -tmp; i <= tmp; i++){
			for (int j = -tmp; j <= tmp; j++){
				setPixel((int)(vec->vec[0] + i), (int)(vec->vec[1] + j), currColor);
			}
		}
	}

	//------------------------------------------------------------------
	//VERTEX BUFFER DRAW POINTS
	//------------------------------------------------------------------

	void drawPoints(){
		Matrix4f* mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size(); i++)
		{
			Vector4f * vec = mat->mulByVec(&vertexBuffer[i]);
			renderPoint(vec);
			delete(vec);
		}
		delete(mat);
	}

	void drawLines(){
		Matrix4f* mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size(); i += 2)
		{
			Vector4f * vec1 = mat->mulByVec(&vertexBuffer[i]);
			Vector4f * vec2 = mat->mulByVec(&vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
	}

	void drawStrip(){
		Matrix4f* mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size() - 1; i++)
		{
			Vector4f * vec1 = mat->mulByVec(&vertexBuffer[i]);
			Vector4f * vec2 = mat->mulByVec(&vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
			delete(vec1);
			delete(vec2);
		}
		delete(mat);
	}

	void drawLoop(){
		Matrix4f* mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size() - 1; i++)
		{
			Vector4f * vec1 = mat->mulByVec(&vertexBuffer[i]);
			Vector4f * vec2 = mat->mulByVec(&vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
			delete(vec1);
			delete(vec2);
		}
		Vector4f * vec1 = mat->mulByVec(&vertexBuffer.back());
		Vector4f * vec2 = mat->mulByVec(&vertexBuffer.front());
		renderLine(vec1, vec2);
		delete(vec1);
		delete(vec2);
		delete(mat);
	}

};

