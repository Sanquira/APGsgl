#ifndef CONTEXT_H_

#define CONTEXT_H_

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
#include <memory>

#include "objects.h"
#include "graphic_primitives.h"

#define FINFINITY std::numeric_limits<float>::infinity()

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class Context {
private:
	Color currColor;
	Color bcgColor;
	float *colorBuffer;
	float *depthBuffer;
	int width, height;
	std::vector<Vector4f> vertexBuffer;
	vector<std::unique_ptr<AbstractPrimitivum>> scenePrimitives;
	vector<std::unique_ptr<PointLight>> lights;
	std::stack<Matrix4f> modelViewStack, projectionStack;
	std::stack<Matrix4f> *currMatrixStack;
	int areaDrawMode, vertexDrawMode;
	float pointSize;
	bool depthTestEnable;
	Matrix4f viewport;
	Material material;

	void symetricPoints(int x, int y, Vector4f * center){
		int xs = (int)center->x;
		int ys = (int)center->y;
		float z = center->z;
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
		pointSize = 1.0f;

		// buffers
		colorBuffer = (float*)calloc(width*height * 3, sizeof(float));
		depthBuffer = (float*)calloc(width*height, sizeof(float));
		fill_n(depthBuffer, width*height, FINFINITY);
		vertexBuffer.reserve(8);
		scenePrimitives.clear();
		lights.clear();

		// matrix stacks
		modelViewStack.push(Matrix4f());
		projectionStack.push(Matrix4f());
		currMatrixStack = &modelViewStack;

		// colors
		bcgColor = Color(0, 0, 0);
		currColor = Color(1, 1, 1);

		// draw modes
		areaDrawMode = 2; //SGL_FILL
	}

	// destructor
	~Context(){
		// buffers
		vertexBuffer.clear();
		scenePrimitives.clear();
		lights.clear();
		free(colorBuffer);
		free(depthBuffer);

		// matrix stacks
		while (!modelViewStack.empty()){
			modelViewStack.pop();
		}
		while (!projectionStack.empty()){
			projectionStack.pop();
		}

		// others
//		delete(bcgColor);
//		delete(currColor);
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
//		delete(currColor);
		currColor = Color(r, g, b);
	}

	Color getBcgColor(){
		return bcgColor;
	}

	void setBcgColor(float r, float g, float b){
//		delete(bcgColor);
		bcgColor = Color(r, g, b);
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
		fill_n(depthBuffer, width*height, FINFINITY);
	}

	Matrix4f computeTransformation(){
		Matrix4f tmp = projectionStack.top().mulByMatrix(modelViewStack.top());
		return viewport.mulByMatrix(tmp);
	}

	//------------------------------------------------------------------
	//RENDERING METHODS
	//------------------------------------------------------------------

	void setPixel(int x, int y, float z, Color clr){
		if (x < 0 || y < 0 || x >= width || y >= height)
			return;
		if (depthBuffer[((y)*width + x)] >= z || !depthTestEnable){
			colorBuffer[((y)*width + x) * 3 + 0] = clr.red;
			colorBuffer[((y)*width + x) * 3 + 1] = clr.green;
			colorBuffer[((y)*width + x) * 3 + 2] = clr.blue;
			depthBuffer[((y)*width + x)] = z;
		}
	}

	Color getPixelColor(int x, int y){
		return  Color(colorBuffer[((y)*width + x) * 3 + 0], colorBuffer[((y)*width + x) * 3 + 1], colorBuffer[((y)*width + x) * 3 + 2]);
	}

	void renderCircle(Vector4f vec, float radii){
		Matrix4f mat = computeTransformation();
		//	mat->toTerminal();
		//	cout << "*******************" << endl;
		float scale = std::sqrt(mat.m[0][0] * mat.m[1][1] - mat.m[1][0] * mat.m[0][1]);
		Vector4f center = mat.mulByVec(vec);
		int x, y, p;
		x = 0;
		y = (int)round(radii*scale);
		p = (int)round(3 - 2*radii*scale); // p pro [0,r]
		while (x < y) {
			symetricPoints(x, y, &center);
			if (p < 0) {
				p += 4 * x + 6;
			}
			else {
				p += 4 * (x - y) + 10;
				y -= 1;
			}
			x += 1;
		}
		if (x == y) {// 45° pixely, jen 4
			symetricPoints(x, y, &center);
		}
	}

	void renderEllipse(Vector4f center, float a, float b){
		Matrix4f mat = computeTransformation();

		for (int i = 0; i < 40; i++)
		{
			float uhel = (float)(M_PI*i / 20.);
			float uhelDalsi = (float)(M_PI*(i + 1) / 20.);
			Vector4f p1, p2;
			p1.x = center.x + a*cos(uhel);
			p1.y = center.y + b*sin(uhel);
			p1.z = center.z;
			p1.w = 1;
			p2.x = center.x + a*cos(uhelDalsi);
			p2.y = center.y + b*sin(uhelDalsi);
			p2.z = center.z;
			p2.w = 1;
			renderLine(mat.mulByVec(p1), mat.mulByVec(p2));
		}
	}

	void renderArc(Vector4f center, float radius, float from, float to){
		Matrix4f mat = computeTransformation();
		
		int numOfVert = (int)ceil(40 * fabs(to - from) / (2 * M_PI)) + 1;
		if((to-from)<0){
			to+=(float)(2*M_PI);
		}
		float step = (to - from) / (numOfVert - 1);

		for (int i = 0; i < numOfVert - 1; i++){
			Vector4f p1, p2;
			p1.x = center.x + cos(from + i*step)*radius;
			p1.y = center.y + sin(from + i*step)*radius;
			p1.z = center.z;
			p1.w = 1;
			p2.x = center.x + cos(from + (i + 1)*step)*radius;
			p2.y = center.y + sin(from + (i + 1)*step)*radius;
			p2.z = center.z;
			p2.w = 1;
			renderLine(mat.mulByVec(p1), mat.mulByVec(p2));
		}
	}

	void renderLine(Vector4f p1, Vector4f p2){
		int x1 = (int)round(p1.x);
		int y1 = (int)round(p1.y);
		int x2 = (int)round(p2.x);
		int y2 = (int)round(p2.y);
		float z1 = p1.z;
		float z2 = p2.z;

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

		int dx =2*( x2 - x1);
		int dy = 2*abs(y2 - y1);

		int error = dx;
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
				setPixel((int)(vec.x + i), (int)(vec.y + j), vec.z, currColor);
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

	//TODO - osetrit vodorovnou usecku
	void drawFilledLoop(){
		std::vector<Vector4f> prus;
		Matrix4f mat = computeTransformation();
		int y0;
		int maxY = INT_MIN;
		int minY = INT_MAX;

		float depth;

		for (size_t i = 0; i < vertexBuffer.size(); i++)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			vec1.homoNorm();
			y0 = (int)ceil(vec1.y);
			if (maxY < y0)
				maxY = y0;
			if (minY > y0)
				minY = y0;
		}

		for (int y = maxY; y >= minY; y--)
		{
			prus.clear();
			for (size_t i = 0; i < vertexBuffer.size(); i++)
			{
				Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
				Vector4f vec2 = mat.mulByVec(vertexBuffer[(i == vertexBuffer.size() - 1) ? 0 : i + 1]);
				vec1.homoNorm();
				vec2.homoNorm();
				if (vec1.y > vec2.y) {
					Vector4f tmp = vec1;
					vec1 = vec2;
					vec2 = tmp;
				}
				Vector4f vec = vec2.minus(vec1);
				float t = (y - vec1.y) / vec.y;
				if (t >= 0 && t < 1){
					float x = vec1.x + vec.x * t;
					float z = vec1.z + vec.z * t;
					prus.push_back(Vector4f(ceil(x), (float)y, z, 1));
				}
			}
			std::sort(prus.begin(), prus.end(), compareVector4fX);

			if (prus.size() != 0){
				
				for (size_t i = 0; i < prus.size() - 1; i += 2)
				{
					float zdd = (prus[i + 1].z - prus[i].z) / (prus[i + 1].x - prus[i].x);
					depth = prus[i].z;
					for (int x = (int)prus[i].x; x < (int)prus[i + 1].x; x++)
					{
						setPixel(x, y, depth, currColor);
						depth += zdd;
					}
				}
			}
		}
		prus.clear();
	}

	void drawArcFilled(Vector4f center, float radius, float from, float to){
		int numOfVert = (int)ceil(40 * fabs(to - from) / (2 * M_PI)) + 1;
		float step = (to - from) / (numOfVert - 1);
		vertexBuffer.clear();
		for (int i = 0; i < numOfVert; i++){
			Vector4f p1, p2;
			p1.x = center.x + cos(from + i*step)*radius;
			p1.y = center.y + sin(from + i*step)*radius;
			p1.z = center.z;
			p1.w = 1;
			vertexBuffer.push_back(p1);
		}
		vertexBuffer.push_back(center);
		drawFilledLoop();
		vertexBuffer.clear();
	}

	void drawEllipseFilled(Vector4f center, float a, float b){
		vertexBuffer.clear();
		for (int i = 0; i < 40; i++)
		{
			float uhel = (float)(M_PI*i / 20.);
			Vector4f p1;
			p1.x = center.x + a*cos(uhel);
			p1.y = center.y + b*sin(uhel);
			p1.z = center.z;
			p1.w = 1;
			vertexBuffer.push_back(p1);
		}
		drawFilledLoop();
		vertexBuffer.clear();
	}

	// TODO - dodelat Scan line seed fill pokud bude chut a potreba nejaky nastrel je v commitu 951b2f86 
	void drawCircleFilled(Vector4f vec, float radii){
		float tmpred = currColor.red;
		currColor.red = -1;
		renderCircle(vec, radii);
		currColor.red = tmpred;
		Matrix4f mat = computeTransformation();
		Vector4f center = mat.mulByVec(vec);
		std::vector<Vector4f> seeds;
		seeds.push_back(center);
		Color clr = Color(-1, -1, -1);

		do{
			Vector4f seed = seeds.back();
			seeds.pop_back();
			// right
			for (int i = (int)seed.x; i < width - 1; i++)
			{
				setPixel(i, (int)seed.y, center.z, currColor);
				//end
				clr = getPixelColor(i + 1, (int)seed.y);
				if (clr.compare(Color(-1, currColor.green, currColor.blue))){
					break;
				}
			}
			// left
			for (size_t i = (int)seed.x; i >= 0; i--)
			{
				setPixel(i, (int)seed.y, center.z, currColor);
				//end
				clr = getPixelColor(i - 1, (int)seed.y);
				if (clr.compare(Color(-1, currColor.green, currColor.blue))){
					break;
				}
			}

			clr = getPixelColor((int)seed.x, (int)seed.y + 1);
			if (!clr.compare(Color(-1, currColor.green, currColor.blue)) && !clr.compare(currColor)){
				seeds.push_back(Vector4f(seed.x, seed.y + 1, 0, 1));
			}

			clr = getPixelColor((int)seed.x, (int)seed.y - 1);
			if (!clr.compare(Color(-1, currColor.green, currColor.blue)) && !clr.compare(currColor)){
				seeds.push_back(Vector4f(seed.x, seed.y - 1, 0, 1));
			}

		} while (!seeds.empty());
	}
	
	//------------------------------------------------------------------
	//SCENE
	//------------------------------------------------------------------
	
	
	
	void initializeScene(){
		scenePrimitives.clear();
		lights.clear();
	}
	
	void addSphere(Vector4f center, float radius){
		scenePrimitives.push_back(::make_unique<SpherePrimitivum>(center,radius,material));
	}

	void addTriangle(){
		if(vertexBuffer.size()==3){
			scenePrimitives.push_back(::make_unique<TrianglePrivitivum>(vertexBuffer[0],vertexBuffer[1],vertexBuffer[2],material));
			vertexBuffer.clear();
		}else{
			cerr << "ERROR!!! addTriangle is NOT triangle!!!" << endl;
		}
	}
	
	void addLight(Vector4f position, Color clr){
		lights.push_back(::make_unique<PointLight>(position,clr));
	}
	
	void setMaterial(Material mat){
		this->material = mat;
	}
	
	void renderRayTrace(){
		Vector4f camera = Vector4f(0,0,0,1);	//eye space
		camera = modelViewStack.top().inverse().mulByVec(camera); //world space
		Matrix4f iv = viewport.inverse();
		Matrix4f ips = projectionStack.top().inverse();
		Matrix4f imv = modelViewStack.top().inverse();
		for (int y=0;y<height;y++){
			for(int x=0;x<width;x++){
				Vector4f pixel = Vector4f((float)x,(float)y,-1,1);
				pixel = iv.mulByVec(pixel);	//normalized space
				pixel = ips.mulByVec(pixel);
				pixel.z = -1;
				pixel.w = 0;	//eye space
				pixel = imv.mulByVec(pixel); //world space
				float minDist = FINFINITY;
				int idxMin = -1;
				for(int i=0;i<(int)scenePrimitives.size();i++){
					float dist = (*scenePrimitives[i]).intersect(camera,pixel);
					if(dist!=FINFINITY && dist<=minDist){
						minDist=dist;
						idxMin = i;
					}
				}
				if(minDist!=FINFINITY){
					Color clr = computePixelColor(scenePrimitives[idxMin],camera,pixel,&lights);
					setPixel(x,y,0,clr);
				}
			}
		}
		
		
	}
	
	Color computePixelColor(std::unique_ptr<AbstractPrimitivum>& primitivum, Vector4f origin, Vector4f ray, vector<std::unique_ptr<PointLight>> *lights){
		Color color;
		Vector4f normalVec = (*primitivum).getNormal();
		Vector4f toCamVec = (*primitivum).getToOrigin(origin);
		for(auto & l : *lights){
			Vector4f toLightVec = (*primitivum).getToLight(l->position);
			bool shade = false;
			float minDist = primitivum->intersect(primitivum->intPoint, toLightVec.reverse())-0.001;
			for (int i = 0; i < scenePrimitives.size(); i++){
				float dist = (*scenePrimitives[i]).intersect(primitivum->intPoint, toLightVec.reverse());
				if (dist < minDist){
					shade = true;
					break;
				}
			}
			float cosa = toLightVec.dotNoHomo(normalVec);
			if (!shade) {
				color.red +=  l->color.red * (*primitivum).material.kd * (*primitivum).material.color.red * cosa;
				color.green +=  l->color.green * (*primitivum).material.kd * (*primitivum).material.color.green * cosa;
				color.blue +=  l->color.blue * (*primitivum).material.kd * (*primitivum).material.color.blue * cosa;
			}
			
			toLightVec = normalVec.mulByConst(2*cosa).minus(toLightVec);
			float cosb = toCamVec.dotNoHomo(toLightVec);
			if(cosb<0){
				cosb = 0;
			}
			color.red +=  l->color.red * (*primitivum).material.ks * pow(cosb,(*primitivum).material.shine);
			color.green +=  l->color.green * (*primitivum).material.ks * pow(cosb,(*primitivum).material.shine);
			color.blue +=  l->color.blue * (*primitivum).material.ks * pow(cosb,(*primitivum).material.shine);
		}
		return color;
	}
	
	
};

#endif

