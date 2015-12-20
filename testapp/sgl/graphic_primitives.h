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
#define NULLVECTOR Vector4f(0,0,0,-1);
#define FEPSILON 1e-2
#define LIGHTAREARAYS 16

using namespace std;

class AbstractLight{

public:
	
	EmissiveMaterial emat;
	
	AbstractLight(){}
	
	virtual ~AbstractLight(){}
	
	virtual vector<Vector4f> getLightPositions() = 0; 
	
	virtual float getLightArea() = 0;
	
	virtual Vector4f getLightNormal() = 0;

};

class PointLight : public AbstractLight {

public:
	
	Vector4f position;
	
	PointLight(Vector4f position, Color clr){
		this->position = position;
		this->emat = EmissiveMaterial(clr.red,clr.green,clr.blue,1,0,0);
	}

	virtual vector<Vector4f> getLightPositions(){
		vector<Vector4f> ret;
		ret.push_back(position);
		return ret;
	}
	
	virtual float getLightArea(){
		return 1;
	}
	
	virtual Vector4f getLightNormal(){
		return NULLVECTOR;
	}

};

class AreaLight : public AbstractLight {
public:
	Vector4f u;
	Vector4f v;
	Vector4f p1;
	
	AreaLight(Vector4f p1, Vector4f p2, Vector4f p3, EmissiveMaterial emat){
		u = p2.minus(p1);
		v = p3.minus(p1);
		this->p1 = p1;
		this->emat = emat;
	}
	
	virtual vector<Vector4f> getLightPositions(){
		vector<Vector4f> ret;
		for(int i=0; i<LIGHTAREARAYS; i++){
			float r1 = (float)rand()/(float)RAND_MAX;
			float r2 = (float)rand()/(float)RAND_MAX;
			float t,s;
			if(r1+r2>1){
				t=1-r1;
				s=1-r2;
			}else{
				t=r1;
				s=r2;
			}
			Vector4f tmp0 = u.mulByConst(-t);
			Vector4f tmp1 = v.mulByConst(-s);
			Vector4f tmp = p1.minus(tmp0).minus(tmp1);
			ret.push_back(tmp);
		}
		return ret;
	}
	
	virtual float getLightArea(){
		return (float) 0.5*u.cross(v).getSize();
	}
	
	virtual Vector4f getLightNormal(){
		Vector4f tmp = u.cross(v);
		tmp.normalize();
		return tmp;
	}

};

class AbstractPrimitivum {

public:

	Material material;
	bool isLight = false;
	
	AbstractPrimitivum(){}

	virtual ~AbstractPrimitivum(){}

	virtual Vector4f getNormal(Vector4f intPoint) = 0;

	virtual Vector4f intersect(Vector4f &origin, Vector4f &ray) = 0;

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
	
	virtual Vector4f getNormal(Vector4f intPoint){
		Vector4f ret = intPoint.minus(center);
		ret.removeHomo();
		ret.normalize();
		return ret;
	}
	
	
	SpherePrimitivum(Vector4f center, float radius, Material material){
		this->center = center;
		this->radius = radius;
		setMaterial(material);
		setTransformationMatrix();
	}
	
	virtual Vector4f intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = transformVector(origin);
		ray.removeHomo();
		Vector4f tr = transformVector(ray);
		
		float a = tr.dotNoHomo(tr);
		float b = to.mulByConst(2).dotNoHomo(tr);
		float c = to.dotNoHomo(to) - 1;

		float D = (b * b) - (4 * a * c);
		if (D < 0) {	// No intersection
			return NULLVECTOR;
		}
		float d1 = (-b + sqrt(D)) / (2 * a);
		float d2 = (-b - sqrt(D)) / (2 * a);

		float ret = d1>=FEPSILON ? (d2>=FEPSILON ? (d1>d2 ? d2 : d1) : d1) : -1;

		if( ret >= FEPSILON ){	// intersect
			Vector4f tmp = tr.mulByConst(-ret);
			tmp = to.minus(tmp);
			return tran.mulByVec(tmp);
		}
		return NULLVECTOR; // no positive intersect
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
	
	virtual Vector4f getNormal(Vector4f intPoint){
		Vector4f u = points[1].minus(points[0]);
		Vector4f v = points[2].minus(points[0]);
		Vector4f ret = u.cross(v);
		ret.normalize();
		ret.removeHomo();
		return ret;
	}

	TrianglePrivitivum(Vector4f v1, Vector4f v2, Vector4f v3, Material material){
		setMaterial(material);
		points[0] = v1;
		points[1] = v2;
		points[2] = v3;
	}

	virtual Vector4f intersect(Vector4f &origin, Vector4f &ray){
		Vector4f to = origin;
		Vector4f tr = ray;
		tr.removeHomo();
		
		Vector4f u = points[1].minus(points[0]);
		Vector4f v = points[2].minus(points[0]);
		Vector4f normalPlane = u.cross(v);
		
		Vector4f tmp = points[0].minus(to);
		float cth1 = normalPlane.dotNoHomo(tr);
		float t = normalPlane.dotNoHomo(tmp)/cth1;

		if(cth1 > 0 || t < 0 || fabs(t) < FEPSILON){
			return NULLVECTOR;
		}
		
		tmp = tr.mulByConst(-t);
		tmp = to.minus(tmp);
		
		if(barycentricInside(tmp)){
			return tmp;
		}
		return NULLVECTOR;
	}

};

#endif

