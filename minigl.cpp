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

vec3 curr_color;
mat4 curr_proj;

struct Vertex {
	//MGLfloat w, x, y, z;
	vec4 vertices;
	// vec3 color;
		
	Vertex() {
		vertices[0] = 0;
		vertices[1] = 0;
		vertices[2] = 0;
		vertices[3] = 0;
	}
	Vertex(MGLfloat x, MGLfloat y, MGLfloat z, MGLfloat w){
		vertices[0] = x;
		vertices[1] = y;
		vertices[2] = z;
		vertices[3] = w;
		// color = m_color;
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

// End of Helpers --------------------------------------------------------------------------

void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {
	// Pixel A
	MGLfloat fi_a = ((tri.a.vertices[0] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_a = ((tri.a.vertices[1] + 1.0) * height * 0.5)- 0.5;
	vec2 Pixel_A = vec2(fi_a, fj_a);

	// Pixel B
	MGLfloat fi_b = ((tri.b.vertices[0] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_b = ((tri.b.vertices[1] + 1.0) * height * 0.5) - 0.5;
	vec2 Pixel_B = vec2(fi_b, fj_b);

	// Pixel C
	MGLfloat fi_c = ((tri.c.vertices[0] + 1.0) * width * 0.5) - 0.5;
	MGLfloat fj_c = ((tri.c.vertices[1] + 1.0) * height * 0.5) - 0.5;
	vec2 Pixel_C = vec2(fi_c, fj_c);
	
	MGLfloat area = getArea(Pixel_A, Pixel_B, Pixel_C);

	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < height; j++) {
			MGLfloat alpha = getArea(vec2(i, j), Pixel_B, Pixel_C) / area;
			MGLfloat beta = getArea(Pixel_A, vec2(i, j), Pixel_C) / area;
			MGLfloat gamma = getArea(Pixel_A, Pixel_B, vec2(i, j)) / area;

			if((alpha >= 0) && (beta >= 0) && (gamma >= 0)) {
				data[i + (j* width)] = Make_Pixel(255, 255, 255);
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
	// Vertex v1 = Vertex(1.0, x, y, 0.0, curr_color);
	// vec_vertex.push_back(v1);
	Vertex vec(x, y, 0.0, 1.0);

	// Multiply
	mult(vec, curr_proj);

	vec_vertex.push_back(vec);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	// Vertex v1 = Vertex(1.0, x, y, z, curr_color);
	// vec_vertex.push_back(v1);
	Vertex vec(x, y, z, 1.0);

	// Multiply
	mult(vec, curr_proj);

	vec_vertex.push_back(vec);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{

}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	
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
	
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	
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
	
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	
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

	curr_proj = {{(2/(right - left)), 0, 0, -(right + left)/(right - left), 
		0, (2/(top - bottom)), 0, -(top + bottom)/(top - bottom),
		0, 0, (-2/(far - near)), -(far + near)/(far - near),
		0, 0, 0, 1}};
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
