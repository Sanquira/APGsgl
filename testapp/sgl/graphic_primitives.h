#ifndef GRAPHIC_PRIMITIVES_H_

#define GRAPHIC_PRIMITIVES_H_

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

#define FINFINITY std::numeric_limits<float>::infinity()

using namespace std;

class PointLight {
private:

public:
	
	Vector4f position;
	Color color;
	
	PointLight(Vector4f position, Color clr){
		this->position = position;
		this->color = clr;
	}


};

class AbstractPrimitivum {
	
protected:
	
	Vector4f toCamVec;

public:

	Material material;
	Vector4f intPoint;

	virtual Vector4f getToOrigin(Vector4f origin) = 0;
	virtual Vector4f getNormal() = 0;
	virtual Vector4f getToLight(Vector4f lightPos) = 0;
	
	AbstractPrimitivum(){}

	virtual ~AbstractPrimitivum(){}

	virtual float intersect(Vector4f &origin, Vector4f &ray) = 0;

	void setMaterial(Material mat){
		material = mat;
	}	
	
	

};

class SpherePrimitivum : public AbstractPrimitivum {
private:

	Vector4f center;
	float radius;

	Matrix4f invTran;
	Matrix4f tran;
	
	void setTransformationMatrix(){
		invTran.m[0][0] = 1/radius;
		invTran.m[1][1] = 1/radius;
		invTran.m[2][2] = 1/radius;
		invTran.m[0][3] = -center.x/radius;
		invTran.m[1][3] = -center.y/radius;
		invTran.m[2][3] = -center.z/radius;
		tran = invTran.inverse();
	}
	
	Vector4f transformVector(Vector4f &vec){
		return invTran.mulByVec(vec);
	}
	
public:
	
	virtual Vector4f getToOrigin(Vector4f origin){
		Vector4f ret = origin.minus(intPoint);
		ret.normalize();
		return ret;
	}
	
	virtual Vector4f getNormal(){
		Vector4f ret = intPoint.minus(center);
		ret.normalize();
		return ret;
	}
	
	virtual Vector4f getToLight(Vector4f lightPos){
		Vector4f ret = lightPos.minus(intPoint);
		ret.normalize();
		return ret;
	}
	
	
	SpherePrimitivum(Vector4f center, float radius, Material material){
		this->center = center;
		this->radius = radius;
		setMaterial(material);
		setTransformationMatrix();
	}
	
	virtual float intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = transformVector(origin);
		ray.removeHomo();
		Vector4f tr = transformVector(ray);
		
		float a = tr.dotNoHomo(tr);
		float b = to.mulByConst(2).dotNoHomo(tr);
		float c = to.dotNoHomo(to) - 1;

		float D = (b * b) - (4 * a * c);
		if (D < 0) {	// No intersection
			return FINFINITY;
		}
		float d1 = (-b + sqrt(D)) / (2 * a);
		float d2 = (-b - sqrt(D)) / (2 * a);

		float ret = min(d1, d2);
		if(ret>=0){	// intersect
			Vector4f tmp = tr.mulByConst(-ret);
			tmp = to.minus(tmp);
			intPoint = tran.mulByVec(tmp);
			return ret;
		}
		return FINFINITY; // no positive intersect
	}

};

class TrianglePrivitivum : public AbstractPrimitivum {
private:
	Vector4f points[3];
	
	bool barycentricInside(Vector4f &point){
		Vector4f u = points[1].minus(points[0]);
		Vector4f v = points[2].minus(points[0]);
		Vector4f w = point.minus(points[0]);
		
		Vector4f vw = v.cross(w);
		Vector4f vu = v.cross(u);
		
		if(vw.dotNoHomo(vu)<0){
			return false;
		}
		
		Vector4f uw = u.cross(w);
		vu = vu.mulByConst(-1);
		
		if(uw.dotNoHomo(vu)<0){
			return false;
		}
		
		float denom = sqrt(vu.dotNoHomo(vu));
		float r = sqrt(vw.dotNoHomo(vw))/denom;
		float t = sqrt(uw.dotNoHomo(uw))/denom;
		
		return (r+t <= 1);
	}

public:

	virtual Vector4f getToOrigin(Vector4f origin){
		Vector4f ret = origin.minus(intPoint);
		ret.normalize();
		return ret;
	}
	
	virtual Vector4f getNormal(){
		Vector4f u = points[1].minus(points[0]);
		Vector4f v = points[2].minus(points[0]);
		Vector4f ret = u.cross(v);
		ret.normalize();
		return ret;
	}
	
	virtual Vector4f getToLight(Vector4f lightPos){
		Vector4f ret = lightPos.minus(intPoint);
		ret.normalize();
		return ret;
	}

	TrianglePrivitivum(Vector4f v1, Vector4f v2, Vector4f v3, Material material){
		setMaterial(material);
		points[0] = v1;
		points[1] = v2;
		points[2] = v3;
	}

	virtual float intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = origin;
		ray.removeHomo();
		Vector4f tr = ray;
		
		Vector4f tmp = points[0].minus(points[1]);
		Vector4f normalPlane = points[2].minus(points[1]).cross(tmp);
		
		tmp = points[0].minus(to);
		float t = normalPlane.dotNoHomo(tmp)/normalPlane.dotNoHomo(tr);
		
		tmp = tr.mulByConst(-t);
		tmp = to.minus(tmp);
		
		if(barycentricInside(tmp)){
			Vector4f tmp = tr.mulByConst(-t);
			intPoint = to.minus(tmp);
			return t;
		}
		return FINFINITY;
	}

};

#endif

