//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
using namespace std;


/// Current error code.
static sglEErrorCode _libStatus = SGL_NO_ERROR;

static inline void setErrCode(sglEErrorCode c)
{
	if (_libStatus == SGL_NO_ERROR)
		_libStatus = c;
}

//---------------------------------------------------------------------------
// sglGetError()
//---------------------------------------------------------------------------
sglEErrorCode sglGetError(void)
{
	sglEErrorCode ret = _libStatus;
	_libStatus = SGL_NO_ERROR;
	return ret;
}

//---------------------------------------------------------------------------
// sglGetErrorString()
//---------------------------------------------------------------------------
const char* sglGetErrorString(sglEErrorCode error)
{
	static const char *errStrigTable[] =
	{
		"Operation succeeded",
		"Invalid argument(s) to a call",
		"Invalid enumeration argument(s) to a call",
		"Invalid call",
		"Quota of internal resources exceeded",
		"Internal library error",
		"Matrix stack overflow",
		"Matrix stack underflow",
		"Insufficient memory to finish the requested operation"
	};

	if ((int)error<(int)SGL_NO_ERROR || (int)error>(int)SGL_OUT_OF_MEMORY) {
		return "Invalid value passed to sglGetErrorString()";
	}

	return errStrigTable[(int)error];
}

//---------------------------------------------------------------------------
// Variables
//---------------------------------------------------------------------------

#define MAXCONTEXT 32

Context **contextBuffer = NULL;
int currContext = -1;

//---------------------------------------------------------------------------
// Initialization functions
//---------------------------------------------------------------------------

void sglInit(void) {
	contextBuffer = new Context*[MAXCONTEXT];
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		contextBuffer[i] = NULL;
	}
}

void sglFinish(void) {
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		contextBuffer[i]->~Context();
	}
	delete contextBuffer;
	cout << "sgl Finished" << endl;
}

// TODO - out of memory exception
int sglCreateContext(int width, int height) {
	int i = 0;
	for (; i < MAXCONTEXT; i++)
	{
		if (contextBuffer[i] == NULL){
			contextBuffer[i] = new Context(width, height);
			return i;
		}
	}
	if (i == MAXCONTEXT){
		throw SGL_OUT_OF_RESOURCES;
	}
	return i;
}


void sglDestroyContext(int id) {
	if (id == currContext){
		throw SGL_INVALID_OPERATION;
	}
	if (id < 0 || id >= MAXCONTEXT || contextBuffer[id] == NULL){
		throw  SGL_INVALID_VALUE;
	}
	contextBuffer[id]->~Context();
	contextBuffer[id] = NULL;
}

void sglSetContext(int id) {
	if (id < 0 || id >= MAXCONTEXT || contextBuffer[id] == NULL){
		throw  SGL_INVALID_VALUE;
	}
	currContext = id;
}

int sglGetContext(void) {
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		if (contextBuffer[i] == NULL){
			throw SGL_INVALID_OPERATION;
			break;
		}
	}
	if (currContext == -1){
		throw SGL_INVALID_OPERATION;
	}
	return currContext;
}

float *sglGetColorBufferPointer(void) {
	if (contextBuffer[currContext] == NULL){
		return NULL;
	}
	return contextBuffer[currContext]->getColorBuffer();
}

//---------------------------------------------------------------------------
// Drawing functions
//---------------------------------------------------------------------------

void sglClearColor(float r, float g, float b, float alpha) {

}

void sglClear(unsigned what) {}

void sglBegin(sglEElementType mode) {}

void sglEnd(void) {}

void sglVertex4f(float x, float y, float z, float w) {}

void sglVertex3f(float x, float y, float z) {}

void sglVertex2f(float x, float y) {}

void sglCircle(float x, float y, float z, float radius) {}

void sglEllipse(float x, float y, float z, float a, float b) {}

void sglArc(float x, float y, float z, float radius, float from, float to) {}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode(sglEMatrixMode mode) {}

void sglPushMatrix(void) {}

void sglPopMatrix(void) {}

void sglLoadIdentity(void) {}

void sglLoadMatrix(const float *matrix) {}

void sglMultMatrix(const float *matrix) {}

void sglTranslate(float x, float y, float z) {}

void sglScale(float scalex, float scaley, float scalez) {}

void sglRotate2D(float angle, float centerx, float centery) {}

void sglRotateY(float angle) {}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {}

void sglViewport(int x, int y, int width, int height) {}

//---------------------------------------------------------------------------
// Attribute functions
//---------------------------------------------------------------------------

void sglColor3f(float r, float g, float b) {}

void sglAreaMode(sglEAreaMode mode) {}

void sglPointSize(float size) {}

void sglEnable(sglEEnableFlags cap) {}

void sglDisable(sglEEnableFlags cap) {}

//---------------------------------------------------------------------------
// RayTracing oriented functions
//---------------------------------------------------------------------------

void sglBeginScene() {}

void sglEndScene() {}

void sglSphere(const float x,
	const float y,
	const float z,
	const float radius) {}

void sglMaterial(const float r,
	const float g,
	const float b,
	const float kd,
	const float ks,
	const float shine,
	const float T,
	const float ior) {}

void sglPointLight(const float x,
	const float y,
	const float z,
	const float r,
	const float g,
	const float b) {}

void sglRayTraceScene() {}

void sglRasterizeScene() {}

void sglEnvironmentMap(const int width,
	const int height,
	float *texels)
{}

void sglEmissiveMaterial(
	const float r,
	const float g,
	const float b,
	const float c0,
	const float c1,
	const float c2
	)
{}

