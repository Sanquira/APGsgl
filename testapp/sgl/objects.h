
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

	bool compare(Color clr){
		if (clr.red == red&&clr.green == green&&clr.blue == blue){
			return true;
		}
		return false;
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

	Vector4f minus(Vector4f &vector){
		Vector4f ret;
		ret.vec[0] = vec[0] - vector.vec[0];
		ret.vec[1] = vec[1] - vector.vec[1];
		ret.vec[2] = vec[2] - vector.vec[2];
		ret.vec[3] = vec[3];
		return ret;
	}

	void homoNorm(){
		vec[0] = vec[0] / vec[3];
		vec[1] = vec[1] / vec[3];
		vec[2] = vec[2] / vec[3];
		vec[3] = vec[3] / vec[3];
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
		//delete tmp;
	}

	Vector4f mulByVec(Vector4f &vec){
		Vector4f ret;
		for (int i = 0; i < 4; i++)
		{
			ret.vec[i] = 0;
			for (int j = 0; j < 4; j++)
			{
				ret.vec[i] += vec.vec[j] * m[i][j];
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

bool compareVector4fX(Vector4f i, Vector4f j) { return (i.vec[0]<j.vec[0]); }

class Context {
private:
	Color *currColor;
	Color *bcgColor;
	float *colorBuffer;
	float *depthBuffer;
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
		float z = center->vec[2];
		setPixel(y + xs, x + ys, z, currColor);
		setPixel(x + xs, y + ys, z, currColor);
		setPixel(-x + xs, y + ys, z, currColor);
		setPixel(-y + xs, x + ys, z, currColor);
		setPixel(-y + xs, -x + ys, z, currColor);
		setPixel(-x + xs, -y + ys, z, currColor);
		setPixel(x + xs, -y + ys, z, currColor);
		setPixel(y + xs, -x + ys, z, currColor);
	}

public:
	// constructor
	Context(int width, int height){
		this->width = width;
		this->height = height;
		colorBuffer = (float*)calloc(width*height * 3, sizeof(float));
		depthBuffer = (float*)calloc(width*height, sizeof(float));
		fill_n(depthBuffer, width*height, std::numeric_limits<float>::infinity());
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
		free(depthBuffer);
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
		Matrix4f mat;
		mat.mulByMatrixToItself(&currMatrixStack->top());
		currMatrixStack->push(mat);
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

	void cleanDepthBuffer(){
		free(depthBuffer);
		depthBuffer = (float*)calloc(width*height, sizeof(float));
		fill_n(depthBuffer, width*height, std::numeric_limits<float>::infinity());
	}

	Matrix4f computeTransformation(){
		Matrix4f tmp = projectionStack.top().mulByMatrix(modelViewStack.top());
		return viewport.mulByMatrix(tmp);
	}

	//------------------------------------------------------------------
	//RENDERING METHODS
	//------------------------------------------------------------------

	void setPixel(int x, int y, float z, Color * clr){
		if (x < 0 || y < 0 || x >= width || y >= height)
			return;
		if (depthBuffer[((y)*width + x)] >= z || !depthTestEnable){
			colorBuffer[((y)*width + x) * 3 + 0] = clr->red;
			colorBuffer[((y)*width + x) * 3 + 1] = clr->green;
			colorBuffer[((y)*width + x) * 3 + 2] = clr->blue;
			depthBuffer[((y)*width + x)] = z;
		}
	}

	Color getPixelColor(int x, int y){
		return  Color(colorBuffer[((y)*width + x) * 3 + 0], colorBuffer[((y)*width + x) * 3 + 1], colorBuffer[((y)*width + x) * 3 + 2]);
	}

	void renderCircle(Vector4f &vec, float radii){
		Matrix4f mat = computeTransformation();
		//	mat->toTerminal();
		//	cout << "*******************" << endl;
		float scale = std::sqrt(mat.m[0][0] * mat.m[1][1] - mat.m[1][0] * mat.m[0][1]);
		Vector4f center = mat.mulByVec(vec);
		int x, y, p;
		x = 0;
		y = (int)(radii*scale);
		p = (int)(1 - radii*scale); // p pro [0,r]
		while (x < y) {
			symetricPoints(x, y, &center);
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
			symetricPoints(x, y, &center);
		}
		//delete(mat);
		//delete(center);
	}

	void renderEllipse(Vector4f &center, float a, float b){
		Matrix4f mat = computeTransformation();

		for (int i = 0; i < 40; i++)
		{
			float uhel = (float)(M_PI*i / 20.);
			float uhelDalsi = (float)(M_PI*(i + 1) / 20.);
			Vector4f p1, p2;
			p1.vec[0] = center.vec[0] + a*cos(uhel);
			p1.vec[1] = center.vec[1] + b*sin(uhel);
			p1.vec[2] = center.vec[2];
			p1.vec[3] = 1;
			p2.vec[0] = center.vec[0] + a*cos(uhelDalsi);
			p2.vec[1] = center.vec[1] + b*sin(uhelDalsi);
			p2.vec[2] = center.vec[2];
			p2.vec[3] = 1;
			renderLine(mat.mulByVec(p1), mat.mulByVec(p2));
			//delete(p1);
			//delete(p2);
		}
		//delete(mat);
	}

	void renderArc(Vector4f &center, float radius, float from, float to){
		Matrix4f mat = computeTransformation();

		int numOfVert = (int)ceil(40 * fabs(to - from) / (2 * M_PI)) + 1;
		float step = (to - from) / (numOfVert - 1);

		for (int i = 0; i < numOfVert - 1; i++){
			Vector4f p1, p2;
			p1.vec[0] = center.vec[0] + cos(from + i*step)*radius;
			p1.vec[1] = center.vec[1] + sin(from + i*step)*radius;
			p1.vec[2] = center.vec[2];
			p1.vec[3] = 1;
			p2.vec[0] = center.vec[0] + cos(from + (i + 1)*step)*radius;
			p2.vec[1] = center.vec[1] + sin(from + (i + 1)*step)*radius;
			p2.vec[2] = center.vec[2];
			p2.vec[3] = 1;
			renderLine(mat.mulByVec(p1), mat.mulByVec(p2));
		}
	}

	void renderLine(Vector4f &p1, Vector4f &p2){
		int x1 = (int)round(p1.vec[0]);
		int y1 = (int)round(p1.vec[1]);
		int x2 = (int)round(p2.vec[0]);
		int y2 = (int)round(p2.vec[1]);
		float z1 = p1.vec[2];
		float z2 = p2.vec[2];

		float depth;

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

		float zdd = (z2 - z1) / (maxX - x1);
		depth = z1;

		for (int x = x1; x <= maxX; x++){
			if (uhel){
				setPixel(y, x, depth, currColor);
			}
			else{
				setPixel(x, y, depth, currColor);
			}
			depth += zdd;
			error -= dy;
			if (error < 0){
				y += ystep;
				error += dx;
			}
		}
	}

	void renderPoint(Vector4f &vec){
		int tmp = (int)(pointSize / 2);
		for (int i = -tmp; i <= tmp; i++){
			for (int j = -tmp; j <= tmp; j++){
				setPixel((int)(vec.vec[0] + i), (int)(vec.vec[1] + j), vec.vec[2], currColor);
			}
		}
	}

	//------------------------------------------------------------------
	//VERTEX BUFFER DRAW POINTS
	//------------------------------------------------------------------

	void drawPoints(){
		Matrix4f mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size(); i++)
		{
			Vector4f  vec = mat.mulByVec(vertexBuffer[i]);
			renderPoint(vec);
		}
	}

	void drawLines(){
		Matrix4f mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size(); i += 2)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			Vector4f vec2 = mat.mulByVec(vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
	}

	void drawStrip(){
		Matrix4f mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size() - 1; i++)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			Vector4f vec2 = mat.mulByVec(vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
	}

	void drawLoop(){
		Matrix4f mat = computeTransformation();
		for (size_t i = 0; i < vertexBuffer.size() - 1; i++)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			Vector4f vec2 = mat.mulByVec(vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
		Vector4f vec1 = mat.mulByVec(vertexBuffer.back());
		Vector4f vec2 = mat.mulByVec(vertexBuffer.front());
		renderLine(vec1, vec2);
	}

	//------------------------------------------------------------------
	//FILLED OBJECTS
	//------------------------------------------------------------------

	//TODO - zaokrouhlovani
	void drawFilledLoop(){
		std::vector<Vector4f> prus;
		Matrix4f mat = computeTransformation();
		int maxY = -INT_MIN, minY = INT_MAX;

		float depth;

		for (size_t i = 0; i < vertexBuffer.size(); i++)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			vec1.homoNorm();
			int y0 = (int)ceil(vec1.vec[1]); //TODO
			(maxY <= y0) ? maxY = y0 : y0;
			(minY >= y0) ? minY = y0 : y0;
		}

		for (int y = maxY; y >= minY; y--)
		{
			prus.clear();
			for (size_t i = 0; i < vertexBuffer.size(); i++)
			{
				Vector4f vec0 = mat.mulByVec(vertexBuffer[(i==0)?vertexBuffer.size()-1:i-1]);
				Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
				Vector4f vec2 = mat.mulByVec(vertexBuffer[(i == vertexBuffer.size() - 1) ? 0 : i + 1]);
				vec0.homoNorm();
				vec1.homoNorm();
				vec2.homoNorm();
				Vector4f vec = vec2.minus(vec1);
				float t = (y - vec1.vec[1]) / vec.vec[1];
				if((vec1.vec[1]<vec0.vec[1]&&vec1.vec[1]<vec2.vec[0])||(vec1.vec[1]>vec0.vec[1]&&vec1.vec[1]>vec2.vec[0])){
					if (t > 0 && t < 1){	//TODO - meze paramteru usecky
						float x = vec1.vec[0] + vec.vec[0] * t;
						float z = vec1.vec[2] + vec.vec[2] * t;
						prus.push_back(Vector4f(ceil(x), (float)y, z, 1)); //TODO
					}
				}else{
					if (t >= 0 && t < 1){	//TODO - meze paramteru usecky
						float x = vec1.vec[0] + vec.vec[0] * t;
						float z = vec1.vec[2] + vec.vec[2] * t;
						prus.push_back(Vector4f(ceil(x), (float)y, z, 1)); //TODO
					}
				}
			}
			std::sort(prus.begin(), prus.end(), compareVector4fX);

			if (prus.size() != 0){
				
				for (size_t i = 0; i < prus.size() - 1; i += 2)
				{
					float zdd = (prus[i + 1].vec[2] - prus[i].vec[2]) / (prus[i + 1].vec[0] - prus[i].vec[0]);
					depth = prus[i].vec[2];
					for (int x = (int)prus[i].vec[0]; x < (int)prus[i + 1].vec[0]; x++)
					{
						setPixel(x, y, depth, currColor);
						depth += zdd;
					}
				}
			}
		}
		prus.clear();
	}

	void drawArcFilled(Vector4f &center, float radius, float from, float to){
		int numOfVert = (int)ceil(40 * fabs(to - from) / (2 * M_PI)) + 1;
		float step = (to - from) / (numOfVert - 1);
		vertexBuffer.clear();
		for (int i = 0; i < numOfVert; i++){
			Vector4f p1, p2;
			p1.vec[0] = center.vec[0] + cos(from + i*step)*radius;
			p1.vec[1] = center.vec[1] + sin(from + i*step)*radius;
			p1.vec[2] = center.vec[2];
			p1.vec[3] = 1;
			vertexBuffer.push_back(p1);
		}
		vertexBuffer.push_back(center);
		drawFilledLoop();
		vertexBuffer.clear();
	}

	void drawEllipseFilled(Vector4f &center, float a, float b){
		vertexBuffer.clear();
		for (int i = 0; i < 40; i++)
		{
			float uhel = (float)(M_PI*i / 20.);
			Vector4f p1;
			p1.vec[0] = center.vec[0] + a*cos(uhel);
			p1.vec[1] = center.vec[1] + b*sin(uhel);
			p1.vec[2] = center.vec[2];
			p1.vec[3] = 1;
			vertexBuffer.push_back(p1);
		}
		drawFilledLoop();
		vertexBuffer.clear();
	}

	// TODO - dodelat Scan line seed fill pokud bude chut a potreba nejaky nastrel je v commitu 951b2f86 
	void drawCircleFilled(Vector4f &vec, float radii){
		float tmpred = currColor->red;
		currColor->red = -1;
		renderCircle(vec, radii);
		currColor->red = tmpred;
		Matrix4f mat = computeTransformation();
		Vector4f center = mat.mulByVec(vec);
		std::vector<Vector4f> seeds;
		seeds.push_back(center);
		Color clr = Color(-1, -1, -1);

		do{
			Vector4f seed = seeds.back();
			seeds.pop_back();
			// right
			for (int i = (int)seed.vec[0]; i < width - 1; i++)
			{
				setPixel(i, (int)seed.vec[1], center.vec[2], currColor);
				//end
				clr = getPixelColor(i + 1, (int)seed.vec[1]);
				if (clr.compare(Color(-1, currColor->green, currColor->blue))){
					break;
				}
			}
			// left
			for (size_t i = (int)seed.vec[0]; i >= 0; i--)
			{
				setPixel(i, (int)seed.vec[1], center.vec[2], currColor);
				//end
				clr = getPixelColor(i - 1, (int)seed.vec[1]);
				if (clr.compare(Color(-1, currColor->green, currColor->blue))){
					break;
				}
			}

			clr = getPixelColor((int)seed.vec[0], (int)seed.vec[1] + 1);
			if (!clr.compare(Color(-1, currColor->green, currColor->blue)) && !clr.compare(*currColor)){
				seeds.push_back(Vector4f(seed.vec[0], seed.vec[1] + 1, 0, 1));
			}

			clr = getPixelColor((int)seed.vec[0], (int)seed.vec[1] - 1);
			if (!clr.compare(Color(-1, currColor->green, currColor->blue)) && !clr.compare(*currColor)){
				seeds.push_back(Vector4f(seed.vec[0], seed.vec[1] - 1, 0, 1));
			}

		} while (!seeds.empty());
	}
};

