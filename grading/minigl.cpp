/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

MGLpoly_mode curr_type;
MGLmatrix_mode curr_mode;

vec3 curr_color;
mat4 curr_proj;
mat4 curr_model;
mat4 *curr_matrix = &curr_proj;


vector<mat4> vec_ModelView;
vector<mat4> vec_ProjMatrix;
vector<mat4> *curr_vec;

vector<vector<MGLfloat>> vec_minZ;

struct Vertex {
	//MGLfloat w, x, y, z;
	vec4 vertices;
	vec3 color;
		
	Vertex() {

		vertices = vec4(0,0,0,1);
		color = vec3(0,0,0);
	}
	Vertex(const vec4 &vec, const vec3 &m_color) {
		vertices = vec;
		color = m_color;
	}
}; 

// Global List of VERTEX
vector<Vertex> vec_vertex;

struct Triangle {
	Vertex a, b, c;

	Triangle() {
		//a = b = c = 0;
	}
	Triangle(Vertex m_a, Vertex m_b, Vertex m_c) : a(m_a), b(m_b), c(m_c){}
};

// Global List of TRIANGLE
vector<Triangle> vec_triangle;

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

// Helper Functions ------------------------------------------------------------------------

MGLfloat getArea(vec2 a, vec2 b, vec2 c) {
	//return (a[0]*b[1] - a[1]*b[0] + b[0]*c[1] - b[1]*c[0] + c[0]*a[1] - c[1]*a[0]) * 0.5;
	return (a[0]*(b[1] - c[1]) + a[1]*(c[0] - b[0]) + (b[0]*c[1] - b[1]*c[0]));
}

void mult(Vertex &vec, mat4 mat) {
	int hold, count = 0;
	MGLfloat sum = 0;
	Vertex temp;

	for(unsigned int i = 0; i < 16; i+= 4) {
		for(unsigned int j = 0; j < 4; j++) {
			sum += vec.vertices[j] * mat.values[i + hold];
			hold++;
		}
		temp.vertices[count] = sum;
		count++;

		hold = 0;
		sum = 0;
	}

	for(unsigned int i = 0; i < 4; i++) {
		vec.vertices[i] = temp.vertices[i];
	}
}

mat4 top_of_active_matrix_stack() {
	if(curr_mode == MGL_PROJECTION) {
		return vec_ProjMatrix.back();
	}
	else if(curr_mode == MGL_MODELVIEW) {
		return vec_ModelView.back();
	}
	else {
		MGL_ERROR("Invalid matrix type.");
		exit(1);
	}
}

Triangle getNDC(const Triangle &tri) {
	Triangle temp = tri;

	// First point
	temp.a.vertices[0] = tri.a.vertices[0] / tri.a.vertices[3];
	temp.a.vertices[1] = tri.a.vertices[1] / tri.a.vertices[3];
	temp.a.vertices[2] = tri.a.vertices[2] / tri.a.vertices[3];

	// Second point
	temp.b.vertices[0] = tri.b.vertices[0] / tri.b.vertices[3];
	temp.b.vertices[1] = tri.b.vertices[1] / tri.b.vertices[3];
	temp.b.vertices[2] = tri.b.vertices[2] / tri.b.vertices[3];

	// Third poing
	temp.c.vertices[0] = tri.c.vertices[0] / tri.c.vertices[3];
	temp.c.vertices[1] = tri.c.vertices[1] / tri.c.vertices[3];
	temp.c.vertices[2] = tri.c.vertices[2] / tri.c.vertices[3];

	// temp.a.vertices = tri.a.vertices * (1/3);
	// temp.b.vertices = tri.b.vertices * (1/3);
	// temp.c.vertices = tri.c.vertices * (1/3);

	return temp;
}

vector<MGLfloat> getColor(MGLfloat alpha, MGLfloat beta, MGLfloat gamma, Triangle tri) {
	MGLfloat divisor = (alpha / tri.a.vertices[3]) + (beta / tri.b.vertices[3]) + (gamma / tri.c.vertices[3]);
	
	vector<MGLfloat> color;
	color.push_back((alpha / tri.a.vertices[3]) / divisor);
	color.push_back((beta / tri.b.vertices[3]) / divisor);
	color.push_back((gamma / tri.c.vertices[3]) / divisor);

	return color;
}

// End of Helpers --------------------------------------------------------------------------

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {

	// Pixel A
	MGLfloat fi_a = ((tri.a.vertices[0] / tri.a.vertices[3] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_a = ((tri.a.vertices[1] / tri.a.vertices[3] + 1.0) * height * 0.5)- 0.5;
	vec2 Pixel_A = vec2(fi_a, fj_a);

	// Pixel B
	MGLfloat fi_b = ((tri.b.vertices[0] / tri.b.vertices[3] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_b = ((tri.b.vertices[1] / tri.b.vertices[3] + 1.0) * height * 0.5) - 0.5;
	vec2 Pixel_B = vec2(fi_b, fj_b);

	// Pixel C
	MGLfloat fi_c = ((tri.c.vertices[0] / tri.c.vertices[3] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_c = ((tri.c.vertices[1] / tri.c.vertices[3] + 1.0) * height * 0.5) - 0.5;
	vec2 Pixel_C = vec2(fi_c, fj_c);
	
	MGLfloat area = getArea(Pixel_A, Pixel_B, Pixel_C);

	// Triangle temp = getNDC(tri);

	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < height; j++) {
			MGLfloat alpha = getArea(vec2(i, j), Pixel_B, Pixel_C) / area;
			MGLfloat beta = getArea(Pixel_A, vec2(i, j), Pixel_C) / area;
			MGLfloat gamma = getArea(Pixel_A, Pixel_B, vec2(i, j)) / area;

			if((alpha >= 0) && (beta >= 0) && (gamma >= 0)) {
				MGLfloat z_depth = alpha * (tri.a.vertices[2]/tri.a.vertices[3]) + beta * (tri.b.vertices[2]/tri.b.vertices[3]) + gamma * (tri.c.vertices[2]/tri.c.vertices[3]);

				if(z_depth >= -1 && z_depth <= 1 && (z_depth < vec_minZ[i][j])) {
					vector<MGLfloat> color = getColor(alpha, beta, gamma, tri);

					MGLfloat red = tri.a.color[0] * color[0] + tri.b.color[0] * color[1] + tri.c.color[0] * color[2];
					MGLfloat green = tri.a.color[1] * color[0] + tri.b.color[1] * color[1] + tri.c.color[1] * color[2];
					MGLfloat blue = tri.a.color[2] * color[0] + tri.b.color[2] * color[1] + tri.c.color[2] * color[2];

					data[i + j * width] = Make_Pixel(255 * red, 255 * green, 255 * blue);
					vec_minZ[i][j] = z_depth;
				}
			}
		}
	}
}

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
					MGLsize height,
                   MGLpixel *data)
{
	Make_Pixel(0, 0, 0);

	// // RESIZING
	vec_minZ.resize(width);

	// Initialize it to a large value (2 is large enough for this case)
	for(unsigned int j = 0; j < width; j++) {
		vec_minZ[j].resize(height);
		for(unsigned int k = 0; k < height; k++) {
			vec_minZ[j][k] = 2;
		}
	}

	for(unsigned int i = 0; i < vec_triangle.size(); i++) {
		Rasterize_Triangle(vec_triangle[i], width, height, data);
	}
	vec_triangle.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	if((mode == MGL_TRIANGLES) || (mode == MGL_QUADS)) {
		curr_type = mode;
	}
	else {
		MGL_ERROR("Invalid Polygon Mode.");
		exit(1);
	}
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(curr_type == MGL_TRIANGLES) {
		if((vec_vertex.size() % 3) == 0) {
			for(unsigned int i = 0; i < vec_vertex.size(); i++) {
				if(((i + 1) % 3) == 0) {
					//cout << i << endl;
					Triangle triang = Triangle(vec_vertex[(i + 1) - 3], vec_vertex[(i+1) - 2], vec_vertex[i]);
					vec_triangle.push_back(triang);
				}
				else {
					// Do Nothing
				}
			}
		}
	}

	else if(curr_type == MGL_QUADS) {
		if((vec_vertex.size() % 4) == 0) {
			for(unsigned int i = 0; i < vec_vertex.size(); i++) {
				if(((i + 1) % 4) == 0) {
					Triangle triang = Triangle(vec_vertex[(i+1)-4], vec_vertex[(i+1) - 2], vec_vertex[(i+1) - 3]);
					vec_triangle.push_back(triang);
					triang = Triangle(vec_vertex[(i+1) - 4], vec_vertex[(i+1) - 2], vec_vertex[i]);
					vec_triangle.push_back(triang);
				}
				else {
					// Do Nothing
				}
			}
		}
	}
	
	// Clean vector
	vec_vertex.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x, y, 0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{

	Vertex vec(curr_proj * curr_model * vec4(x,y,z,1), curr_color);
	// // // Multiply
	// mult(vec, curr_proj);
	// mult(vec, curr_model);

	vec_vertex.push_back(vec);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	if(mode == MGL_PROJECTION) {
		curr_matrix = &curr_proj;
		curr_vec = &vec_ProjMatrix;
	}
	else if (mode == MGL_MODELVIEW) {
		curr_matrix = &curr_model;
		curr_vec = &vec_ModelView;
	}
	else {
		MGL_ERROR("Invalid matrix mode.");
		exit(1);
	}
	curr_mode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if(curr_mode == MGL_PROJECTION || curr_mode == MGL_MODELVIEW) {
		curr_vec->push_back(*curr_matrix);
	}
	else {
		MGL_ERROR("Invalid matrix mode.");
		exit(1);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(!curr_vec->empty()) {
		*curr_matrix = curr_vec->back();
		curr_vec->pop_back();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	*curr_matrix = {{1, 0, 0, 0, 
					0, 1, 0 , 0,
					0, 0 , 1, 0,
					0, 0, 0, 1}};
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	mat4 temp;
	temp.make_zero();

	for(unsigned int i = 0; i < 4; i++) {
		for(unsigned int j = 0; j < 4; j++) {
			temp(i, j) = *(matrix + (i + j * 4));
		}
	}

	*curr_matrix = temp;
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	mat4 temp;
	temp.make_zero();

	for(unsigned int i = 0; i < 4; i++) {
		for(unsigned int j = 0; j < 4; j++) {
			temp(i, j) = *(matrix + (i + j * 4));
		}
	}

	*curr_matrix = *curr_matrix * temp;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
					MGLfloat y,
					MGLfloat z)
{
	mat4 translate = {{1, 0, 0, 0,
						0, 1, 0, 0,
						0, 0, 1, 0,
						x, y, z, 1}};

	*curr_matrix = *curr_matrix * translate;
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
				MGLfloat x,
				MGLfloat y,
				MGLfloat z)
{
	MGLfloat c = cos(angle * M_PI /180);
	MGLfloat s = sin(angle * M_PI/180);

	vec3 temp = vec3(x,y,z).normalized();

	mat4 rotate = {{temp[0]*temp[0]*(1 - c) + c, temp[1]*temp[0]*(1 - c) + temp[2]*s, temp[0]*temp[2]*(1 - c) - temp[1]*s, 0,
					temp[0]*temp[1]*(1 - c) - temp[2]*s, temp[1]*temp[1]*(1 - c) + c, temp[1]*temp[2]*(1 - c) + temp[0]*s, 0,
					temp[0]*temp[2]*(1 - c) + temp[1]*s, temp[1]*temp[2]*(1 - c) - temp[0]*s, temp[2]*temp[2]*(1 - c) + c, 0,
					0, 0, 0, 1}};

	*curr_matrix = *curr_matrix * rotate;
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
				MGLfloat y,
				MGLfloat z)
{
	mat4 scalar = {{x, 0, 0, 0,
					0, y, 0, 0,
					0, 0, z, 0,
					0, 0, 0, 1}};

	mglMultMatrix(scalar.values);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	MGLfloat A = (right + left) / (right - left);
	MGLfloat B = (top + bottom) / (top - bottom);
	MGLfloat C = -(far + near) / (far - near);
	MGLfloat D = - (2 * far * near) / (far - near);

	// corrected
	mat4 frustum = {{(2 * near) / (right - left), 0, 0, 0,
					0, (2 * near)/ (top - bottom), 0, 0,
					A, B, C, -1,
					0, 0, D, 0}};

	*curr_matrix = *curr_matrix * frustum;
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
				MGLfloat right,
				MGLfloat bottom,
				MGLfloat top,
				MGLfloat near,
				MGLfloat far)
{
	if(curr_mode == MGL_PROJECTION) {
		// corrected
		mat4 ortho = {{(2/(right - left)), 0, 0, 0,
						0, (2/(top - bottom)), 0, 0,
						0, 0, (-2/(far - near)), 0,
						-(right + left)/(right - left),  -(top + bottom)/(top - bottom), -(far + near)/(far - near), 1}};

		*curr_matrix = *curr_matrix * ortho;
	}
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
				MGLfloat green,
				MGLfloat blue)
{
	curr_color[0] = red;
	curr_color[1] = green;
	curr_color[2] = blue;
}
