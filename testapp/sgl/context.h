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

//#define FINFINITY std::numeric_limits<float>::infinity()
#define MAX_RECURSION 8

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class Context {
private:
	int width, height;
	float pointSize;
	Color currColor, bcgColor;
	float *colorBuffer, *depthBuffer;

	Matrix4f viewport;

	std::stack<Matrix4f> *currMatrixStack;
	std::stack<Matrix4f> modelViewStack, projectionStack;

	std::vector<Vector4f> vertexBuffer;
	int areaDrawMode, vertexDrawMode;

	bool depthTestEnable;

	vector<std::unique_ptr<AbstractPrimitivum>> scenePrimitives;
	vector<std::unique_ptr<AbstractLight>> lights;

	Material material;
	EmissiveMaterial emissiveMaterial;

	int envMapWidth, envMapHeight;
	float *envMap;
	bool envMapSet;

	/*
		8-fold symmetry for Bresenham's algorithm
	*/
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
	Context(int width, int height) {
		this->width = width;
		this->height = height;
		pointSize = 1.0f;

		// buffers
		colorBuffer = (float*) calloc(width*height * 3, sizeof(float));
		depthBuffer = (float*) calloc(width*height, sizeof(float));
		fill_n(depthBuffer, width*height, FINFINITY);
		vertexBuffer.reserve(8);
		scenePrimitives.clear();
		lights.clear();

		envMapSet = false;
		envMap = (float*) calloc(1, sizeof(float));

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
	~Context() {
		// buffers
		vertexBuffer.clear();
		scenePrimitives.clear();
		lights.clear();
		free(colorBuffer);
		free(depthBuffer);
		free(envMap);

		// matrix stacks
		while (!modelViewStack.empty()) {
			modelViewStack.pop();
		}
		while (!projectionStack.empty()) {
			projectionStack.pop();
		}
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
		currColor = Color(r, g, b);
	}

	Color getBcgColor(){
		return bcgColor;
	}

	void setBcgColor(float r, float g, float b){
		bcgColor = Color(r, g, b);
	}

	void clearVertexBuffer(){
		vertexBuffer.clear();
	}

	std::vector<Vector4f> getVertexBuffer(){
		return vertexBuffer;
	}

	void addToVertexBuffer(Vector4f vec) {
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

	void setAreaDrawMode(int mode) {
		areaDrawMode = mode;
	}

	int getAreaDrawMode() {
		return areaDrawMode;
	}

	void setVertexDrawMode(int mode) {
		vertexDrawMode = mode;
	}

	int getVertexDrawMode() {
		return vertexDrawMode;
	}

	void setPointSize(float size) {
		pointSize = size;
	}

	float getPointSize() {
		return pointSize;
	}

	void setDepthTest(bool ena) {
		depthTestEnable = ena;
	}

	bool getDepthTest() {
		return depthTestEnable;
	}

	Matrix4f getViewport() {
		return viewport;
	}

	void setViewport(Matrix4f mat) {
		memcpy(&viewport, &mat, sizeof(Matrix4f));
	}

	void cleanDepthBuffer() {
		//free(depthBuffer);
		//depthBuffer = (float*)calloc(width*height, sizeof(float));
		fill_n(depthBuffer, width*height, FINFINITY);
	}

	Matrix4f computeTransformation() {
		Matrix4f tmp = projectionStack.top().mulByMatrix(modelViewStack.top());
		return viewport.mulByMatrix(tmp);
	}

	bool isEnvMapSet() {
		return envMapSet;
	}



	//------------------------------------------------------------------
	//RENDERING METHODS
	//------------------------------------------------------------------

	void setPixel(int x, int y, float z, Color clr) {
		if (x < 0 || y < 0 || x >= width || y >= height)
			return;

		int dPixel = y * width + x;
		int cPixel = dPixel * 3;

		if (depthBuffer[dPixel] >= z || !depthTestEnable) {
			colorBuffer[cPixel + 0] = clr.red;
			colorBuffer[cPixel + 1] = clr.green;
			colorBuffer[cPixel + 2] = clr.blue;
			depthBuffer[dPixel] = z;
		}
	}

	Color getPixelColor(int x, int y) {
		int pixel = (y * width + x) * 3;

		return Color(colorBuffer[pixel + 0], colorBuffer[pixel + 1], colorBuffer[pixel + 2]);
	}

	void renderCircle(Vector4f vec, float radii) {
		Matrix4f mat = computeTransformation();
		Vector4f center = mat.mulByVec(vec);
		float scale = std::sqrt(mat.m[0][0] * mat.m[1][1] - mat.m[1][0] * mat.m[0][1]);

		int x, y, p;
		x = 0;
		y = (int) round(radii*scale);
		p = (int) round(3 - 2*radii*scale); // p pro [0,r]

		while (x < y) {
			symetricPoints(x, y, &center);

			if (p < 0) {
				p += 4 * x + 6;
			} else {
				p += 4 * (x - y) + 10;
				y -= 1;
			}
			
			x += 1;
		}

		if (x == y) {// 45° pixely, jen 4
			symetricPoints(x, y, &center);
		}
	}

	void renderEllipse(Vector4f center, float a, float b) {
		Matrix4f mat = computeTransformation();

		for (int i = 0; i < 40; i++) {
			float uhel = (float) (M_PI*i / 20.);
			float uhelDalsi = (float) (M_PI*(i + 1) / 20.);

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

	void renderArc(Vector4f center, float radius, float from, float to) {
		Matrix4f mat = computeTransformation();
		float len = to - from;

		if (len < 0) {
			to += (float) (2*M_PI);
		}

		int numOfVert = (int) ceil(40 * fabs(len) / (2 * M_PI));
		float step = len / numOfVert;

		for (int i = 0; i < numOfVert; i++) {
			Vector4f p1, p2;
			p1.x = center.x + cos(from + i*step) * radius;
			p1.y = center.y + sin(from + i*step) * radius;
			p1.z = center.z;
			p1.w = 1;

			p2.x = center.x + cos(from + (i + 1)*step) * radius;
			p2.y = center.y + sin(from + (i + 1)*step) * radius;
			p2.z = center.z;
			p2.w = 1;

			renderLine(mat.mulByVec(p1), mat.mulByVec(p2));
		}
	}

	void renderLine(Vector4f p1, Vector4f p2) {
		int x1 = (int) round(p1.x);
		int y1 = (int) round(p1.y);
		int x2 = (int) round(p2.x);
		int y2 = (int) round(p2.y);
		float z1 = p1.z;
		float z2 = p2.z;

		bool uhel = (abs(y2 - y1) > abs(x2 - x1)); //>45°
		if (uhel) {
			std::swap(x1, y1);
			std::swap(x2, y2);
		}

		if (x1 > x2) { //leva polorovina
			std::swap(x1, x2);
			std::swap(y1, y2);
		}

		int maxX = x2;

		int dx = 2 *    (x2 - x1);
		int dy = 2 * abs(y2 - y1);

		int error = dx;
		int ystep = (y1 < y2) ? 1 : -1;
		int y = y1;

		float depth = z1;
		float zdd = (z2 - z1) / (maxX - x1);

		for (int x = x1; x <= maxX; x++) {
			if (uhel) {
				setPixel(y, x, depth, currColor);
			} else {
				setPixel(x, y, depth, currColor);
			}

			depth += zdd;
			error -= dy;

			if (error < 0) {
				y += ystep;
				error += dx;
			}
		}
	}

	void renderPoint(Vector4f &vec) {
		int tmp = (int) (pointSize / 2);

		for (int i = -tmp; i <= tmp; i++) {
			for (int j = -tmp; j <= tmp; j++) {
				setPixel((int)(vec.x + i), (int)(vec.y + j), vec.z, currColor);
			}
		}
	}



	//------------------------------------------------------------------
	//VERTEX BUFFER DRAW POINTS
	//------------------------------------------------------------------

	void drawPoints() {
		Matrix4f mat = computeTransformation();

		for (size_t i = 0; i < vertexBuffer.size(); i++)
		{
			Vector4f vec = mat.mulByVec(vertexBuffer[i]);
			renderPoint(vec);
		}
	}

	void drawLines() {
		Matrix4f mat = computeTransformation();

		for (size_t i = 0; i < vertexBuffer.size(); i += 2)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			Vector4f vec2 = mat.mulByVec(vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
	}

	void drawStrip() {
		Matrix4f mat = computeTransformation();

		for (size_t i = 0; i < vertexBuffer.size() - 1; i++)
		{
			Vector4f vec1 = mat.mulByVec(vertexBuffer[i]);
			Vector4f vec2 = mat.mulByVec(vertexBuffer[i + 1]);
			renderLine(vec1, vec2);
		}
	}

	void drawLoop() {
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

	void drawFilledLoop() {
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

			y0 = (int) ceil(vec1.y);
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

				if (t >= 0 && t < 1) {
					float x = vec1.x + vec.x * t;
					float z = vec1.z + vec.z * t;
					prus.push_back(Vector4f(ceil(x), (float) y, z, 1));
				}
			}
			std::sort(prus.begin(), prus.end(), compareVector4fX);

			if (prus.size() != 0) {
				for (size_t i = 0; i < prus.size() - 1; i += 2)
				{
					float zdd = (prus[i + 1].z - prus[i].z) / (prus[i + 1].x - prus[i].x);
					depth = prus[i].z;
					for (int x = (int) prus[i].x; x < (int) prus[i + 1].x; x++)
					{
						setPixel(x, y, depth, currColor);
						depth += zdd;
					}
				}
			}
		}

		prus.clear();
	}

	void drawArcFilled(Vector4f center, float radius, float from, float to) {
		float len = to - from;
		int numOfVert = (int) ceil(40 * fabs(len) / (2 * M_PI));
		float step = len / numOfVert;
		vertexBuffer.clear();

		for (int i = 0; i <= numOfVert; i++)
		{
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

	void drawEllipseFilled(Vector4f center, float a, float b) {
		vertexBuffer.clear();

		for (int i = 0; i < 40; i++)
		{
			float uhel = (float) (M_PI*i / 20.);
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

	void drawCircleFilled(Vector4f vec, float radii) {
		float tmpred = currColor.red;
		currColor.red = -1;
		renderCircle(vec, radii);
		currColor.red = tmpred;

		Matrix4f mat = computeTransformation();
		Vector4f center = mat.mulByVec(vec);
		std::vector<Vector4f> seeds;
		seeds.push_back(center);
		Color clr = Color(-1, -1, -1);

		do {
			Vector4f seed = seeds.back();
			seeds.pop_back();

			// right
			for (int i = (int) seed.x; i < width - 1; i++)
			{
				setPixel(i, (int) seed.y, center.z, currColor);
				//end
				clr = getPixelColor(i + 1, (int) seed.y);
				if (clr.compare(Color(-1, currColor.green, currColor.blue))) {
					break;
				}
			}

			// left
			for (size_t i = (int)seed.x; i >= 0; i--)
			{
				setPixel(i, (int)seed.y, center.z, currColor);
				//end
				clr = getPixelColor(i - 1, (int)seed.y);
				if (clr.compare(Color(-1, currColor.green, currColor.blue))) {
					break;
				}
			}

			clr = getPixelColor((int)seed.x, (int)seed.y + 1);
			if (!clr.compare(Color(-1, currColor.green, currColor.blue)) && !clr.compare(currColor)) {
				seeds.push_back(Vector4f(seed.x, seed.y + 1, 0, 1));
			}

			clr = getPixelColor((int)seed.x, (int)seed.y - 1);
			if (!clr.compare(Color(-1, currColor.green, currColor.blue)) && !clr.compare(currColor)) {
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
	
	void addSphere(Vector4f center, float radius) {
		scenePrimitives.push_back(::make_unique<SpherePrimitivum>(center,radius,material));
	}

	void addTriangle() {
		if (vertexBuffer.size() == 3) {
			scenePrimitives.push_back(::make_unique<TrianglePrivitivum>(vertexBuffer[0],vertexBuffer[1],vertexBuffer[2],material));

			if (emissiveMaterial.c0 != -1) {
				lights.push_back(::make_unique<AreaLight>(vertexBuffer[0],vertexBuffer[1],vertexBuffer[2],emissiveMaterial));
				scenePrimitives.back()->material.color = emissiveMaterial.color;
				scenePrimitives.back()->isLight = true;
			}

			vertexBuffer.clear();
		} else {
			cerr << "ERROR!!! addTriangle is NOT triangle!!!" << endl;
		}
	}
	
	void addLight(Vector4f position, Color clr) {
		lights.push_back(::make_unique<PointLight>(position,clr));
	}
	
	void setMaterial(Material mat) {
		this->material = mat;
	}
	
	void setEmissiveMaterial(EmissiveMaterial mat) {
		this->emissiveMaterial = mat;
	}
	
	void setEnviromentMap(int width, int height, float *texels) {
		free(envMap);
		envMapWidth = width;
		envMapHeight = height;
		envMap = texels;
		envMapSet = true;
	}
	
	void renderRayTrace() {
		Vector4f camera = Vector4f(0,0,0,1);	//eye space
		camera = modelViewStack.top().inverse().mulByVec(camera); //world space
		Matrix4f iv = viewport.inverse();
		Matrix4f ips = projectionStack.top().inverse();
		Matrix4f imv = modelViewStack.top().inverse();

		for (int y = 0; y < height; y++) {
			for(int x = 0; x < width; x++) {
				Vector4f pixel = Vector4f((float) x, (float) y, 0, 1);
				pixel = iv.mulByVec(pixel);	//normalized space
				pixel = ips.mulByVec(pixel);
				pixel.z = -1;
				pixel.w = 0;	//eye space
				pixel = imv.mulByVec(pixel); //world space
				pixel.normalize();

				Color clr = computeRecursionColor(0, camera, pixel);

				setPixel(x, y, 0, clr);
			}
		}


	}

private:
	Color computeRecursionColor(int depth, Vector4f origin, Vector4f ray) {
		Color clr;
		Vector4f normalVec;
		Vector4f intPoint;

		// find intersection point
		float intDist = FINFINITY;
		int idxMin = -1;
		for (size_t i = 0; i < scenePrimitives.size(); i++) {
			Vector4f tmp = scenePrimitives[i]->intersect(origin, ray);
			float dist = origin.minus(tmp).getSize();

			if (tmp.w != -1 && dist <= intDist) {
				intPoint = tmp;
				intDist = dist;
				idxMin = i;
			}
		}

		// no intersection
		if (intDist == FINFINITY) {
			if (isEnvMapSet()) {
				clr = getEnvMapColor(ray);
			}

			return clr;
		}

		// is light
		if (scenePrimitives[idxMin]->isLight) {
			clr = scenePrimitives[idxMin]->material.color;
			return clr;
		}

		// compute phong
		clr = computePixelColor(scenePrimitives[idxMin], intPoint, origin, ray, &lights);

		// steps exceeded
		if (depth >= MAX_RECURSION) {
			return clr;
		}

		normalVec = scenePrimitives[idxMin]->getNormal(intPoint);

		// compute reflected ray
		if (scenePrimitives[idxMin]->material.ks != 0) {
			Vector4f refRay = reflectRay(ray, normalVec);

			Color ret = computeRecursionColor(++depth, intPoint, refRay);
			clr.add(ret, scenePrimitives[idxMin]->material.ks);
		}

		// compute refracted ray
		if (scenePrimitives[idxMin]->material.T != 0) {
			Vector4f refRay = refractRay(ray, normalVec, scenePrimitives[idxMin]->material.ior);

			// find refract out 
			Vector4f tmp = scenePrimitives[idxMin]->intersect(intPoint, refRay);
			if (tmp.w != -1) {
				refRay = refractRay(refRay, scenePrimitives[idxMin]->getNormal(tmp), scenePrimitives[idxMin]->material.ior);
				if (refRay.w != -1) {
					Color ret = computeRecursionColor(++depth, tmp, refRay);
					clr.add(ret, scenePrimitives[idxMin]->material.T*scenePrimitives[idxMin]->material.T);
				}
			}
		}

		return clr;
	}

	Vector4f reflectRay(Vector4f ray, Vector4f normal) {
		float cth1 = ray.dotNoHomo(normal);

		Vector4f ret = normal.mulByConst(2 * cth1);
		ret = ray.minus(ret);
		ret.normalize();

		return ret;
	}

	Vector4f refractRay(Vector4f ray, Vector4f normal, float ior) {
		float pior;
		float cth1 = ray.dotNoHomo(normal);

		if (cth1 < 0) {
			// from outside
			pior = 1.f / ior;
		} else {
			// from inside
			pior = ior;
			cth1 = -cth1;
			normal = normal.reverse();
		}

		float cth2 = (float) (1. - (pior*pior) * (1. - (cth1*cth1)));
		if (cth2 < 0) {
			// end
			return Vector4f(0, 0, 0, -1);
		}

		cth2 = pior*cth1 + sqrt(cth2);
		Vector4f ret = normal.mulByConst(-cth2);
		Vector4f tmp = ray.mulByConst(pior).reverse();
		ret = ret.minus(tmp);
		ret.normalize();

		return ret;
	}
	
	Color getEnvMapColor(Vector4f ray) {
		float d = sqrt(ray.x*ray.x + ray.y*ray.y);
		float r = (d > 0) ? (float) (acos(ray.z) / (2*M_PI*d)) : 0;
		float s = 0.5f + r * ray.x;
		float t = 0.5f - r * ray.y;
		
		int x = (int) (s * envMapWidth);
		int y = (int) (t * envMapHeight);
		int pixel = (y * envMapWidth + x) * 3;
		
		return Color(envMap[pixel + 0], envMap[pixel + 1], envMap[pixel + 2]);
	}
	
	Color computePixelColor(std::unique_ptr<AbstractPrimitivum>& primitivum, Vector4f intPoint, Vector4f origin, Vector4f ray, vector<std::unique_ptr<AbstractLight>> *lights) {
		Color color;
		Vector4f normalVec = primitivum->getNormal(intPoint);
		Vector4f toCamVec = origin.minus(intPoint);
		toCamVec.normalize();
		float prColorRate = 1;

		for (auto & l : *lights) {
			vector<Vector4f> positions = l->getLightPositions();
			Vector4f lightNormal = l->getLightNormal();

			if (positions.size() == 1) {
				lightNormal = intPoint.minus(positions[0]);
				lightNormal.normalize();
			}

			lightNormal = lightNormal.reverse();
			float lightAreaRatio = l->getLightArea() / positions.size();

			for (size_t lp = 0; lp < positions.size(); lp++) {
		
				Vector4f toLightVec = positions[lp].minus(intPoint);
				
				if (toLightVec.getSize() <= FEPSILON) {
					cout << "t" << endl;
				}
				
				float distToLight = toLightVec.getSize();
				toLightVec.normalize();
				bool skipLight = false;

				for (size_t i = 0; i < scenePrimitives.size(); i++) {
					Vector4f tmp = scenePrimitives[i]->intersect(intPoint, toLightVec);
					float dist = intPoint.minus(tmp).getSize();

					if ((tmp.w != -1) && (fabs(dist - distToLight) >= FEPSILON) && (dist < distToLight)) {
						skipLight=true;
						break;
					}
				}

				if (!skipLight) {
					Color clrTmp;
					float cosa = toLightVec.dotNoHomo(normalVec);
					clrTmp.red   += prColorRate * l->emat.color.red   * primitivum->material.kd * primitivum->material.color.red * cosa;
					clrTmp.green += prColorRate * l->emat.color.green * primitivum->material.kd * primitivum->material.color.green * cosa;
					clrTmp.blue  += prColorRate * l->emat.color.blue  * primitivum->material.kd * primitivum->material.color.blue * cosa;
				
					Vector4f tmp = normalVec.mulByConst(2*cosa).minus(toLightVec);
					float cosb = toCamVec.dotNoHomo(tmp);

					if (cosb < 0) {
						cosb = 0;
					}
					
					clrTmp.red   +=  l->emat.color.red   * primitivum->material.ks * pow(cosb,primitivum->material.shine);
					clrTmp.green +=  l->emat.color.green * primitivum->material.ks * pow(cosb,primitivum->material.shine);
					clrTmp.blue  +=  l->emat.color.blue  * primitivum->material.ks * pow(cosb,primitivum->material.shine);
					
					float cosPhi = lightNormal.dotNoHomo(toLightVec);
					float denum = (l->emat.c0) + (l->emat.c1) * distToLight + (l->emat.c2) * distToLight * distToLight;
					
					color.red   += clrTmp.red   * cosPhi * lightAreaRatio / denum;
					color.green += clrTmp.green * cosPhi * lightAreaRatio / denum;
					color.blue  += clrTmp.blue  * cosPhi * lightAreaRatio / denum;
				}
			}
		}

		return color;
	}
	
};

#endif

