//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "objects.h"
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
		if (contextBuffer[i] != NULL) {
			delete(contextBuffer[i]);
			contextBuffer[i] = NULL;
		}
	}
	delete[] contextBuffer;
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
	delete(contextBuffer[id]);
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
				con->setPixel(i, j, -std::numeric_limits<float>::infinity(), con->getBcgColor());
			}
		}
	}
	if ((what&SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT){
		con->cleanDepthBuffer();
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
	if (con->getVertexDrawMode() == SGL_POLYGON){
		if (con->getAreaDrawMode() == SGL_FILL){
			con->drawFilledLoop();
		}
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
	Vector4f v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = 1;
	if (contextBuffer[currContext]->getAreaDrawMode() == SGL_FILL){
		contextBuffer[currContext]->drawCircleFilled(v, radius);
	}
	else{
		contextBuffer[currContext]->renderCircle(v, radius);
	}
}

void sglEllipse(float x, float y, float z, float a, float b) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (a < 0 || b < 0){
		throw SGL_INVALID_VALUE;
	}
	Vector4f v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = 1;
	if (contextBuffer[currContext]->getAreaDrawMode() == SGL_FILL){
		contextBuffer[currContext]->drawEllipseFilled(v, a, b);
	}
	contextBuffer[currContext]->renderEllipse(v, a, b);
	if ((contextBuffer[currContext]->getVertexDrawMode() & SGL_POINT) == SGL_POINT){
		contextBuffer[currContext]->renderPoint(v);
	}
}

void sglArc(float x, float y, float z, float radius, float from, float to) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (radius < 0){
		throw SGL_INVALID_VALUE;
	}
	Vector4f v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = 1;
	if (contextBuffer[currContext]->getAreaDrawMode() == SGL_FILL){
		contextBuffer[currContext]->drawArcFilled(v, radius, from, to);
	}

	contextBuffer[currContext]->renderArc(v, radius, from, to);
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

void sglMultMatrix(const float *matrix) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	Matrix4f mul(matrix);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&mul);
}

void sglTranslate(float x, float y, float z) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	Matrix4f tr;
	tr.m[0][3] = x;
	tr.m[1][3] = y;
	tr.m[2][3] = z;
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&tr);
}

void sglScale(float scalex, float scaley, float scalez) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	Matrix4f sc;
	sc.m[0][0] = scalex;
	sc.m[1][1] = scaley;
	sc.m[2][2] = scalez;
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&sc);
}

void sglRotate2D(float angle, float centerx, float centery) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	sglTranslate(centerx, centery, 0);

	Matrix4f rotZ;
	rotZ.m[0][0] = cos(angle);
	rotZ.m[1][0] = sin(angle);
	rotZ.m[0][1] = -sin(angle);
	rotZ.m[1][1] = cos(angle);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&rotZ);

	sglTranslate(-centerx, -centery, 0);
}

void sglRotateY(float angle) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}

	Matrix4f rotY;
	rotY.m[0][0] = cos(angle);
	rotY.m[0][2] = -sin(angle);
	rotY.m[2][0] = sin(angle);
	rotY.m[2][2] = cos(angle);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&rotY);
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (right == left || top == bottom || far == near){
		throw SGL_INVALID_VALUE;
	}
	Matrix4f mat;
	mat.m[0][0] = 2 / (right - left);
	mat.m[1][1] = 2 / (top - bottom);
	mat.m[2][2] = 2 / (far - near);
	mat.m[0][3] = -(right + left) / (right - left);
	mat.m[1][3] = -(top + bottom) / (top - bottom);
	mat.m[2][3] = -(far + near) / (far - near);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&mat);
}

void sglFrustum(float left, float right, float bottom, float top, float near, float far) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (near <= 0 || far <= 0 || left == right || top == bottom || far == near){
		throw SGL_INVALID_VALUE;
	}
	Matrix4f mat;
	mat.m[0][0] = (2 * near) / (right - left);
	mat.m[0][2] = (right + left) / (right - left);
	mat.m[1][1] = (2 * near) / (top - bottom);
	mat.m[1][2] = (top + bottom) / (top - bottom);
	mat.m[2][2] = -(far + near) / (far - near);
	mat.m[2][3] = -(2*far * near) / (far - near);
	mat.m[3][2] = -1;
	mat.m[3][3] = 0;
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&mat);
}

void sglViewport(int x, int y, int width, int height) {
	if (transactionEnabled || contextBuffer[currContext] == NULL){
		throw SGL_INVALID_OPERATION;
	}
	if (width < 0 || height < 0){
		throw SGL_INVALID_VALUE;
	}
	Matrix4f mat;
	mat.m[0][0] = (float)width / 2;
	mat.m[1][1] = (float)height / 2;
	mat.m[0][3] = x + (float)width / 2;
	mat.m[1][3] = y + (float)height / 2;
	contextBuffer[currContext]->setViewport(mat);
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

	if ((mode & SGL_LINE) == SGL_LINE){
		contextBuffer[currContext]->setAreaDrawMode(SGL_LINE);
		return;
	}
	if ((mode & SGL_FILL) == SGL_FILL){
		contextBuffer[currContext]->setAreaDrawMode(SGL_FILL);
		return;
	}
	contextBuffer[currContext]->setAreaDrawMode(SGL_POINT);
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

