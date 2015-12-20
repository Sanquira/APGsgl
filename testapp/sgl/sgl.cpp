//---------------------------------------------------------------------------
// sgl.cpp
// Empty implementation of the SGL (Simple Graphics Library)
// Date:  2011/11/1
// Author: Jaroslav Krivanek, Jiri Bittner CTU Prague
//---------------------------------------------------------------------------

#include "sgl.h"
#include "objects.h"
#include "graphic_primitives.h"
#include "context.h"

//#define FINFINITY std::numeric_limits<float>::infinity()

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

bool transactionFragmentEnabled = false;
bool transactionSceneEnabled = false;


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
			if (contextBuffer[i] == NULL){
				setErrCode(SGL_OUT_OF_MEMORY);
				return -1;
				}

			return i;
		}
	}
	if (i == MAXCONTEXT){
		setErrCode(SGL_OUT_OF_RESOURCES);
		return -1;
	}
	return i;
}

void sglDestroyContext(int id) {
	if (id == currContext){
		setErrCode( SGL_INVALID_OPERATION);
		return;
	}
	if (id < 0 || id >= MAXCONTEXT || contextBuffer[id] == NULL){
		setErrCode(  SGL_INVALID_VALUE);
		return;
	}
	delete(contextBuffer[id]);
	contextBuffer[id] = NULL;
}

void sglSetContext(int id) {
	if (id < 0 || id >= MAXCONTEXT || contextBuffer[id] == NULL){
		setErrCode(SGL_INVALID_VALUE);
		return;
	}
	currContext = id;
}

int sglGetContext(void) {
	for (int i = 0; i < MAXCONTEXT; i++)
	{
		if (contextBuffer[i] == NULL){
			setErrCode( SGL_INVALID_OPERATION);
			break;
		}
	}
	if (currContext == -1){
		setErrCode( SGL_INVALID_OPERATION);
		return -1;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode( SGL_INVALID_OPERATION);
	}
	contextBuffer[currContext]->setBcgColor(r, g, b);
}

void sglClear(unsigned what) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((what & ~(SGL_COLOR_BUFFER_BIT | SGL_DEPTH_BUFFER_BIT)) != 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
	}
	Context *con = contextBuffer[currContext];
	if ((what&SGL_COLOR_BUFFER_BIT) == SGL_COLOR_BUFFER_BIT){
		for (int i = 0; i < con->getWidth(); i++)
		{
			for (int j = 0; j < con->getHeight(); j++)
			{
				con->setPixel(i, j, -FINFINITY, con->getBcgColor());
			}
		}
	}
	if ((what&SGL_DEPTH_BUFFER_BIT) == SGL_DEPTH_BUFFER_BIT){
		con->cleanDepthBuffer();
	}
}

void sglBegin(sglEElementType mode) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((mode & ~(SGL_POINTS | SGL_LINES | SGL_LINE_STRIP | SGL_LINE_LOOP | SGL_TRIANGLES | SGL_POLYGON | SGL_AREA_LIGHT | SGL_LAST_ELEMENT_TYPE)) != 0){
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	contextBuffer[currContext]->setVertexDrawMode(mode);
	transactionFragmentEnabled = true;
}

void sglEnd(void) {
	if (transactionFragmentEnabled != true){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	transactionFragmentEnabled = false;
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
		if(transactionSceneEnabled){
			con->addTriangle();
		}else{
			if (con->getAreaDrawMode() == SGL_FILL){
				con->drawFilledLoop();
			}
			con->drawLoop();
		}
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (radius < 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (a < 0 || b < 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (radius < 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((mode & ~(SGL_MODELVIEW | SGL_PROJECTION)) != 0){
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	if ((mode&SGL_MODELVIEW) == SGL_MODELVIEW){
		contextBuffer[currContext]->setMatrixMode(0);
	}
	if ((mode&SGL_PROJECTION) == SGL_PROJECTION){
		contextBuffer[currContext]->setMatrixMode(1);
	}

}

void sglPushMatrix(void) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->pushMatrix();
}

void sglPopMatrix(void) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (contextBuffer[currContext]->getMatrixStack()->size() == 1){
		setErrCode(SGL_STACK_UNDERFLOW);
		return;
	}
	contextBuffer[currContext]->popMatrix();
}

void sglLoadIdentity(void) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->getMatrixStack()->top() = Matrix4f();
}

void sglLoadMatrix(const float *matrix) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->getMatrixStack()->top() = Matrix4f(matrix);
}

void sglMultMatrix(const float *matrix) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	Matrix4f mul(matrix);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&mul);
}

void sglTranslate(float x, float y, float z) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	Matrix4f tr;
	tr.m[0][3] = x;
	tr.m[1][3] = y;
	tr.m[2][3] = z;
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&tr);
}

void sglScale(float scalex, float scaley, float scalez) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	Matrix4f sc;
	sc.m[0][0] = scalex;
	sc.m[1][1] = scaley;
	sc.m[2][2] = scalez;
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&sc);
}

void sglRotate2D(float angle, float centerx, float centery) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}

	Matrix4f rotY;
	rotY.m[0][0] = cos(angle);
	rotY.m[0][2] = -sin(angle);
	rotY.m[2][0] = sin(angle);
	rotY.m[2][2] = cos(angle);
	contextBuffer[currContext]->getMatrixStack()->top().mulByMatrixToItself(&rotY);
}

void sglOrtho(float left, float right, float bottom, float top, float near, float far) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (right == left || top == bottom || far == near){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (near <= 0 || far <= 0 || left == right || top == bottom || far == near){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (width < 0 || height < 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->setDrawColor(r, g, b);
}

//Sets the current drawing mode of !!!CLOSED!!! areas for subsequent operations.
void sglAreaMode(sglEAreaMode mode) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((mode & ~(SGL_POINT | SGL_LINE | SGL_FILL)) != 0){
		setErrCode(SGL_INVALID_ENUM);
		return;
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
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if (size < 0){
		setErrCode(SGL_INVALID_VALUE);
		return;
	}
	contextBuffer[currContext]->setPointSize(size);
}

void sglEnable(sglEEnableFlags cap) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((cap & ~(SGL_DEPTH_TEST)) != 0){
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	if ((cap&SGL_DEPTH_TEST) == SGL_DEPTH_TEST){
		contextBuffer[currContext]->setDepthTest(true);
	}
}

void sglDisable(sglEEnableFlags cap) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	if ((cap & ~(SGL_DEPTH_TEST)) != 0){
		setErrCode(SGL_INVALID_ENUM);
		return;
	}
	if ((cap&SGL_DEPTH_TEST) == SGL_DEPTH_TEST){
		contextBuffer[currContext]->setDepthTest(false);
	}
}

//---------------------------------------------------------------------------
// RayTracing oriented functions
//---------------------------------------------------------------------------

void sglBeginScene() {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->initializeScene();
	transactionSceneEnabled = true;
}

void sglEndScene() {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	transactionSceneEnabled = false;
}

void sglSphere(const float x, const float y, const float z, const float radius) {
	if (!transactionSceneEnabled || transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	
	contextBuffer[currContext]->addSphere(Vector4f(x,y,z,1), radius);
}

void sglMaterial(const float r,
	const float g,
	const float b,
	const float kd,
	const float ks,
	const float shine,
	const float T,
	const float ior) {
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	Material mat = Material(r,g,b,kd,ks,shine,T,ior);
	contextBuffer[currContext]->setMaterial(mat);
}

void sglPointLight(const float x,
	const float y,
	const float z,
	const float r,
	const float g,
	const float b) {
	if (!transactionSceneEnabled || transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->addLight(Vector4f(x,y,z,1),Color(r,g,b));
	}

void sglRayTraceScene() {
	if (transactionSceneEnabled || transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->renderRayTrace();

}

//NOT USED
void sglRasterizeScene() {}

/// Environment map specification.
/**
  Sets the HDR environment map defining the "background" using a rectangular
  texture. If defined it replaces the background color (set with sglClearColor())
  for both primary and secondary rays.

  @param width [in] texture width
  @param height [in] texture height
  @param texels [in] texture elements (width*height RGB float triplets)

  ERRORS:
   - SGL_INVALID_OPERATION
    No context has been allocated yet or sglEnvironmentMap() is called within a
    sglBegin() / sglEnd() sequence.
*/
void sglEnvironmentMap(const int width,
	const int height,
	float *texels)
{
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	contextBuffer[currContext]->setEnviromentMap(width,height,texels);
}

void sglEmissiveMaterial(
	const float r,
	const float g,
	const float b,
	const float c0,
	const float c1,
	const float c2
	){
	if (transactionFragmentEnabled || contextBuffer[currContext] == NULL){
		setErrCode(SGL_INVALID_OPERATION);
		return;
	}
	EmissiveMaterial mat = EmissiveMaterial(r,g,b,c0,c1,c2);
	contextBuffer[currContext]->setEmissiveMaterial(mat);
}




















