//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
//#include <vld.h>
#include <objects.h>
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

bool transactionEnabled = false;


//---------------------------------------------------------------------------
// Initialization functions
//---------------------------------------------------------------------------

//TODO - mam pocit ze tohle nestaci
void sglInit(void) {
	contextBuffer = new Context*[MAXCONTEXT];
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		contextBuffer[i] = NULL;
	}
}

//TODO - ti same zde
void sglFinish(void) {
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		if (contextBuffer[i] != NULL)
			delete(contextBuffer[i]);
	}
	delete [] contextBuffer;
	cout << "sgl Finished" << endl;
}

int sglCreateContext(int width, int height) {
	int i = 0;
	for (; i < MAXCONTEXT; i++)
	{
		if (contextBuffer[i] == NULL){
			contextBuffer[i] = new Context(width, height);
			if (contextBuffer[i] == NULL)
				throw SGL_OUT_OF_MEMORY;

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
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->setBcgColor(r, g, b);
}

void sglClear(unsigned what) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((what & ~(SGL_COLOR_BUFFER_BIT | SGL_DEPTH_BUFFER_BIT)) != 0){
		throw SGL_INVALID_VALUE;
	}
	Context *con = contextBuffer[currContext];
	if ((what&SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT){
		for (int i = 0; i < con->getWidth(); i++)
		{
			for (int j = 0; j < con->getHeight(); j++)
			{
				con->setPixel(i, j, *(con->getBcgColor()));
			}
		}
	}
	if ((what&SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT){
		// TODO - implementovat az bude implementovany depth buffer
	}
}

void sglBegin(sglEElementType mode) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((mode & ~(SGL_POINTS | SGL_LINES | SGL_LINE_STRIP | SGL_LINE_LOOP | SGL_TRIANGLES | SGL_POLYGON | SGL_AREA_LIGHT | SGL_LAST_ELEMENT_TYPE)) != 0){
		throw SGL_INVALID_ENUM;
	}
	contextBuffer[currContext]->setVertexDrawMode(mode);
	transactionEnabled = true;
}

void sglEnd(void) {
	if (transactionEnabled != true){
		throw SGL_INVALID_OPERATION;
	}
	transactionEnabled = false;
	Context* con = contextBuffer[currContext];
	if (con->getVertexDrawMode() == SGL_POINTS){
		con->drawPoints();
	}
	if (con->getVertexDrawMode() == SGL_LINES){
		con->drawLines();
	}
	if (con->getVertexDrawMode() == SGL_LINE_STRIP){
		con->drawStrip();
	}
	if (con->getVertexDrawMode() == SGL_LINE_LOOP){
		con->drawLoop();
	}
	con->clearVertexBuffer();
}

void sglVertex4f(float x, float y, float z, float w) {
	contextBuffer[currContext]->addToVertexBuffer(Vector4f(x, y, z, w));
}

void sglVertex3f(float x, float y, float z) {
	sglVertex4f(x, y, z, 1);
}

void sglVertex2f(float x, float y) {
	sglVertex4f(x, y, 0, 1);
}

void sglCircle(float x, float y, float z, float radius) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (radius < 0){
		throw SGL_INVALID_VALUE;
	}
	contextBuffer[currContext]->renderCircle(Vector4f(x, y, z, 1), radius);
}

void sglEllipse(float x, float y, float z, float a, float b) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (a < 0 || b < 0){
		throw SGL_INVALID_VALUE;
	}
	contextBuffer[currContext]->renderEllipse(Vector4f(x, y, z, 1), a, b);
	if ((contextBuffer[currContext]->getVertexDrawMode() & SGL_POINT) == SGL_POINT){
		contextBuffer[currContext]->renderPoint(Vector4f(x, y, z, 1));
	}
}

void sglArc(float x, float y, float z, float radius, float from, float to) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (radius < 0){
		throw SGL_INVALID_VALUE;
	}
	contextBuffer[currContext]->renderArc(Vector4f(x, y, z, 1), radius, from, to);
}

//---------------------------------------------------------------------------
// Transform functions
//---------------------------------------------------------------------------

void sglMatrixMode(sglEMatrixMode mode) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((mode & ~(SGL_MODELVIEW | SGL_PROJECTION)) != 0){
		throw SGL_INVALID_ENUM;
	}
	if ((mode&SGL_MODELVIEW) == SGL_MODELVIEW){
		contextBuffer[currContext]->setMatrixMode(0);
	}
	if ((mode&SGL_PROJECTION) == SGL_PROJECTION){
		contextBuffer[currContext]->setMatrixMode(1);
	}

}

void sglPushMatrix(void) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->pushMatrix();
}

void sglPopMatrix(void) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (contextBuffer[currContext]->getMatrixStack()->size() == 1){
		throw SGL_STACK_UNDERFLOW;
	}
	contextBuffer[currContext]->popMatrix();
}

void sglLoadIdentity(void) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->getMatrixStack()->top() = Matrix4f();
}

void sglLoadMatrix(const float *matrix) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->getMatrixStack()->top() = Matrix4f(matrix);
}

//TODO - podivne instrukce
void sglMultMatrix(const float *matrix) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(&Matrix4f(matrix));
}

void sglTranslate(float x, float y, float z) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	Matrix4f* tr = new Matrix4f();
	tr->getMatrix()[0][3] = x;
	tr->getMatrix()[1][3] = y;
	tr->getMatrix()[2][3] = z;
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(tr);
	delete(tr);
}

void sglScale(float scalex, float scaley, float scalez) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	Matrix4f* sc = new Matrix4f();
	sc->getMatrix()[0][0] = scalex;
	sc->getMatrix()[1][1] = scaley;
	sc->getMatrix()[2][2] = scalez;
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(sc);
	delete(sc);
}

void sglRotate2D(float angle, float centerx, float centery) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	sglTranslate(centerx, centery, 0);

	Matrix4f* rotZ = new Matrix4f();
	rotZ->getMatrix()[0][0] = cos(angle);
	rotZ->getMatrix()[1][0] = sin(angle);
	rotZ->getMatrix()[0][1] = -sin(angle);
	rotZ->getMatrix()[1][1] = cos(angle);
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(rotZ);
	delete(rotZ);

	sglTranslate(-centerx, -centery, 0);
}

void sglRotateY(float angle) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}

	Matrix4f* rotY = new Matrix4f();
	rotY->getMatrix()[0][0] = cos(angle);
	rotY->getMatrix()[0][2] = -sin(angle);
	rotY->getMatrix()[0][2] = sin(angle);
	rotY->getMatrix()[2][2] = cos(angle);
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(rotY);
	delete(rotY);
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (right == left || top == bottom || far == near){
		throw SGL_INVALID_VALUE;
	}
	Matrix4f* mat = new Matrix4f();
	mat->getMatrix()[0][0] = 2 / (right - left);
	mat->getMatrix()[1][1] = 2 / (top - bottom);
	mat->getMatrix()[2][2] = 2 / (far - near);
	mat->getMatrix()[0][3] = -(right + left) / (right - left);
	mat->getMatrix()[1][3] = -(top + bottom) / (top - bottom);
	mat->getMatrix()[2][3] = -(far + near) / (far - near);
	contextBuffer[currContext]->getMatrix().mulByMatrixToItself(mat);
	delete(mat);
}

//TODO - zatim neni potreba
void sglFrustum(float left, float right, float bottom, float top, float near, float far) {
	cout << "sglFrustum need to be implemented" << endl;
}

void sglViewport(int x, int y, int width, int height) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (width < 0 || height < 0){
		throw SGL_INVALID_VALUE;
	}
	Matrix4f* mat = new Matrix4f();
	mat->getMatrix()[0][0] = (float)width / 2;
	mat->getMatrix()[1][1] = (float)height / 2;
	mat->getMatrix()[0][3] = x + (float)width / 2;
	mat->getMatrix()[1][3] = y + (float)height / 2;
	contextBuffer[currContext]->setViewport(*mat);
	delete(mat);
}

//---------------------------------------------------------------------------
// Attribute functions
//---------------------------------------------------------------------------

void sglColor3f(float r, float g, float b) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	contextBuffer[currContext]->setDrawColor(r, g, b);
}

//Sets the current drawing mode of !!!CLOSED!!! areas for subsequent operations.
void sglAreaMode(sglEAreaMode mode) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((mode & ~(SGL_POINT | SGL_LINE | SGL_FILL)) != 0){
		throw SGL_INVALID_ENUM;
	}
	if ((mode&SGL_POINT) == SGL_POINT){
		contextBuffer[currContext]->setAreaDrawMode(SGL_POINT);
		return;
	}
	if ((mode&SGL_LINE) == SGL_LINE){
		contextBuffer[currContext]->setAreaDrawMode(SGL_LINE);
		return;
	}
	contextBuffer[currContext]->setAreaDrawMode(SGL_FILL);
}

void sglPointSize(float size) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (size < 0){
		throw SGL_INVALID_VALUE;
	}
	contextBuffer[currContext]->setPointSize(size);
}

void sglEnable(sglEEnableFlags cap) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((cap & ~(SGL_DEPTH_TEST)) != 0){
		throw SGL_INVALID_ENUM;
	}
	if ((cap&SGL_DEPTH_TEST) == SGL_DEPTH_TEST){
		contextBuffer[currContext]->setDepthTest(true);
	}
}

void sglDisable(sglEEnableFlags cap) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if ((cap & ~(SGL_DEPTH_TEST)) != 0){
		throw SGL_INVALID_ENUM;
	}
	if ((cap&SGL_DEPTH_TEST) == SGL_DEPTH_TEST){
		contextBuffer[currContext]->setDepthTest(false);
	}
}

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

