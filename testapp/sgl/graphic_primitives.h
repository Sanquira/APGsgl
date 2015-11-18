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

using namespace std;

class AbstractPrimitivum {
	
protected:
	
	Matrix4f invTran;
	
	
	virtual void upgradeTransformationMatrix() = 0;

public:

	void setTransformationMatrix(Matrix4f mat){
		Matrix4f tmp = mat.inverse();
		invTran = invTran.mulByMatrix(tmp);
	}
	
	Matrix4f getInverseTransformationMatrix(){
		return invTran;
	}
	
	Vector4f transformVector(Vector4f &vec){
		return invTran.mulByVec(vec);
	}
	
	virtual ~AbstractPrimitivum(){}

	virtual void intersect(Vector4f &origin, Vector4f &ray) = 0;
	

};

class SpherePrimitivum : public AbstractPrimitivum {
private:
	Vector4f center;
	float radius;
	
	void upgradeTransformationMatrix(){
		invTran.m[0][0] = 1/radius;
		invTran.m[1][1] = 1/radius;
		invTran.m[2][2] = 1/radius;
		invTran.m[0][3] = -center.x;
		invTran.m[1][3] = -center.y;
		invTran.m[2][3] = -center.z;
	}
	
public:
	
	SpherePrimitivum(Vector4f center, float radius){
		this->center = center;
		this->radius = radius;
		upgradeTransformationMatrix();
	}
	
	virtual void intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = transformVector(origin);
		ray.removeHomo();
		Vector4f tr = transformVector(ray);
		
		float a = tr.dotNoHomo(tr);
		float b = to.mulByConst(2).dotNoHomo(tr);
		float c = to.dotNoHomo(to) - 0.25;

		float D = (b * b) - (4 * a * c);
		if (D < 0) {
			cout << "No intersection" << endl;
		}
		float d1 = (-b + sqrt(D)) / (2 * a);
		float d2 = (-b - sqrt(D)) / (2 * a);

		float ret = min(d1, d2);
		if(ret>=0){			
			cout << "intersect in " << ret << endl;
		}else{
			cout << "No positive intersection " << d1 << ", " << d2 <<endl;
		}
	}

};

class TrianglePrivitivum : public AbstractPrimitivum {
private:
	Vector4f points[3];

	void upgradeTransformationMatrix(){}
	
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

	TrianglePrivitivum(Vector4f v1, Vector4f v2, Vector4f v3){
		points[0] = v1;
		points[1] = v2;
		points[2] = v3;
	}

	virtual void intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = transformVector(origin);
		ray.removeHomo();
		Vector4f tr = transformVector(ray);
		
		Vector4f tmp = points[0].minus(points[1]);
		Vector4f normalPlane = points[2].minus(points[1]).cross(tmp);
		
		tmp = points[0].minus(to);
		float t = normalPlane.dotNoHomo(tmp)/normalPlane.dotNoHomo(tr);
		
		tmp = tr.mulByConst(-t);
		tmp = to.minus(tmp);
		
		if(barycentricInside(tmp)){
			cout << "inside" << endl;
		}else{
			cout << "not inside" << endl;
		}
		
	}

};

#endif

