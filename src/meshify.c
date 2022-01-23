#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#ifdef HAVE_ZLIB
	#include <zlib.h>
	#ifdef HAVE_JSON
		#include "cJSON.h"
	#endif
#endif
#include "meshify.h"
#include "base64.h" //required for GIfTI
#include "nifti1.h"
#include "bwlabel.h"
#include "radixsort.h"

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

typedef struct {
	vec3d p000;
	double val[8];
} GRIDCELL;

#define ABS(x) (x < 0 ? -(x) : (x))

double sqr(double x) {
	return x * x;
}

double dx(vec3d p0, vec3d p1) {
	return sqrt( sqr(p0.x - p1.x) + sqr(p0.y - p1.y) + sqr(p0.z - p1.z));
}


vec3d setVec3d(double x, double y,double z) {
	return (vec3d){.x = x, .y = y, .z = z};
}

/*-------------------------------------------------------------------------
	Return the point between two points in the same ratio as
	isolevel is between valp1 and valp2
*/
vec3d VertexInterp(double isolevel, vec3d p1,vec3d p2,double valp1,double valp2)
{
	double mu;
	vec3d p;
	#ifdef MARCHING_CUBE_CONDITIONALS
	//pipelined computation faster than branched conditionals?
	if (ABS(isolevel-valp1) < 0.00001)
		return(p1);
	if (ABS(isolevel-valp2) < 0.00001)
		return(p2);
	if (ABS(valp1-valp2) < 0.00001)
		return(p1); //<- in theory, avoids divide by zero, but isolevel crossing detected
	#endif
	mu = (isolevel - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);
	return(p);
}

/*-------------------------------------------------------------------------
	Given a grid cell and an isolevel, calculate the triangular
	facets requied to represent the isosurface through the cell.
	Return the number of triangular facets, the array "triangles"
	will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
	of totally below the isolevel.
*/
int PolygoniseCube(GRIDCELL g,double iso,vec3d *tri)
{
	int cubeindex;
	vec3d vertlist[12];
/*
	int edgeTable[256].  It corresponds to the 2^8 possible combinations of
	of the eight (n) vertices either existing inside or outside (2^n) of the
	surface.  A vertex is inside of a surface if the value at that vertex is
	less than that of the surface you are scanning for.  The table index is
	constructed bitwise with bit 0 corresponding to vertex 0, bit 1 to vert
	1.. bit 7 to vert 7.  The value in the table tells you which edges of
	the table are intersected by the surface.  Once again bit 0 corresponds
	to edge 0 and so on, up to edge 12.
	Constructing the table simply consisted of having a program run thru
	the 256 cases and setting the edge bit if the vertices at either end of
	the edge had different values (one is inside while the other is out).
	The purpose of the table is to speed up the scanning process.  Only the
	edges whose bit's are set contain vertices of the surface.
	Vertex 0 is on the bottom face, back edge, left side.
	The progression of vertices is clockwise around the bottom face
	and then clockwise around the top face of the cube.  Edge 0 goes from
	vertex 0 to vertex 1, Edge 1 is from 2->3 and so on around clockwise to
	vertex 0 again. Then Edge 4 to 7 make up the top face, 4->5, 5->6, 6->7
	and 7->4.  Edge 8 thru 11 are the vertical edges from vert 0->4, 1->5,
	2->6, and 3->7.
		 4--------5     *---4----*
		/|       /|    /|       /|
	  / |      / |   7 |      5 |
	 /  |     /  |  /  8     /  9
	7--------6   | *----6---*   |
	|   |    |   | |   |    |   |
	|   0----|---1 |   *---0|---*
	|  /     |  /  11 /     10 /
	| /      | /   | 3      | 1
	|/       |/    |/       |/
	3--------2     *---2----*
*/
int edgeTable[256]={
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

/*
	int triTable[256][16] also corresponds to the 256 possible combinations
	of vertices.
	The [16] dimension of the table is again the list of edges of the cube
	which are intersected by the surface.  This time however, the edges are
	enumerated in the order of the vertices making up the triangle mesh of
	the surface.  Each edge contains one vertex that is on the surface.
	Each triple of edges listed in the table contains the vertices of one
	triangle on the mesh.  The are 16 entries because it has been shown that
	there are at most 5 triangles in a cube and each "edge triple" list is
	terminated with the value -1.
	For example triTable[3] contains
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
	This corresponds to the case of a cube whose vertex 0 and 1 are inside
	of the surface and the rest of the verts are outside (00000001 bitwise
	OR'ed with 00000010 makes 00000011 == 3).  Therefore, this cube is
	intersected by the surface roughly in the form of a plane which cuts
	edges 8,9,1 and 3.  This quadrilateral can be constructed from two
	triangles: one which is made of the intersection vertices found on edges
	1,8, and 3; the other is formed from the vertices on edges 9,8, and 1.
	Remember, each intersected edge contains only one surface vertex.  The
	vertex triples are listed in counter clockwise order for proper facing.
*/
int triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

	/*
		Determine the index into the edge table which
		tells us which vertices are inside of the surface
	*/
	cubeindex = 0;
	if (g.val[0] < iso) cubeindex |= 1;
	if (g.val[1] < iso) cubeindex |= 2;
	if (g.val[2] < iso) cubeindex |= 4;
	if (g.val[3] < iso) cubeindex |= 8;
	if (g.val[4] < iso) cubeindex |= 16;
	if (g.val[5] < iso) cubeindex |= 32;
	if (g.val[6] < iso) cubeindex |= 64;
	if (g.val[7] < iso) cubeindex |= 128;
	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		return(0);
	vec3d p[8];
	p[0] = setVec3d(g.p000.x  ,g.p000.y  ,g.p000.z  );
	p[1] = setVec3d(g.p000.x+1,g.p000.y  ,g.p000.z  );
	p[2] = setVec3d(g.p000.x+1,g.p000.y+1,g.p000.z  );
	p[3] = setVec3d(g.p000.x  ,g.p000.y+1,g.p000.z  );
	p[4] = setVec3d(g.p000.x  ,g.p000.y  ,g.p000.z+1);
	p[5] = setVec3d(g.p000.x+1,g.p000.y  ,g.p000.z+1);
	p[6] = setVec3d(g.p000.x+1,g.p000.y+1,g.p000.z+1);
	p[7] = setVec3d(g.p000.x  ,g.p000.y+1,g.p000.z+1);
	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
		vertlist[0] = VertexInterp(iso,p[0],p[1],g.val[0],g.val[1]);
	if (edgeTable[cubeindex] & 2)
		vertlist[1] = VertexInterp(iso,p[1],p[2],g.val[1],g.val[2]);
	if (edgeTable[cubeindex] & 4)
		vertlist[2] = VertexInterp(iso,p[2],p[3],g.val[2],g.val[3]);
	if (edgeTable[cubeindex] & 8)
		vertlist[3] = VertexInterp(iso,p[3],p[0],g.val[3],g.val[0]);
	if (edgeTable[cubeindex] & 16)
		vertlist[4] = VertexInterp(iso,p[4],p[5],g.val[4],g.val[5]);
	if (edgeTable[cubeindex] & 32)
		vertlist[5] = VertexInterp(iso,p[5],p[6],g.val[5],g.val[6]);
	if (edgeTable[cubeindex] & 64)
		vertlist[6] = VertexInterp(iso,p[6],p[7],g.val[6],g.val[7]);
	if (edgeTable[cubeindex] & 128)
		vertlist[7] = VertexInterp(iso,p[7],p[4],g.val[7],g.val[4]);
	if (edgeTable[cubeindex] & 256)
		vertlist[8] = VertexInterp(iso,p[0],p[4],g.val[0],g.val[4]);
	if (edgeTable[cubeindex] & 512)
		vertlist[9] = VertexInterp(iso,p[1],p[5],g.val[1],g.val[5]);
	if (edgeTable[cubeindex] & 1024)
		vertlist[10] = VertexInterp(iso,p[2],p[6],g.val[2],g.val[6]);
	if (edgeTable[cubeindex] & 2048)
		vertlist[11] = VertexInterp(iso,p[3],p[7],g.val[3],g.val[7]);
	int nvert = 0;
	/* Create the triangles */
	for (int i=0;triTable[cubeindex][i]!=-1;i+=3) {
		tri[nvert] = vertlist[triTable[cubeindex][i  ]];
		nvert++;
		tri[nvert] = vertlist[triTable[cubeindex][i+1]];
		nvert++;
		tri[nvert] = vertlist[triTable[cubeindex][i+2]];
		nvert++;
	}
	return(nvert);
}

int unify_vertices(vec3d **inpt, vec3i *tris, int ntri, bool verbose) {
	//"vertex welding": reduces the number of vertices, number of faces unchanged
	double startTime = clockMsec();
	vec3d *pts = *inpt;
	int npt = ntri * 3;
	int *old2new = (int *)malloc(npt * sizeof(int));
	vec3d pt0 = pts[0];
	float* dx_in = (float*)malloc(npt*sizeof(float));
	float* dx_out = (float*)malloc(npt*sizeof(float));
	uint32_t* idx_in = (uint32_t*)malloc(npt*sizeof(uint32_t));
	uint32_t* idx_out = (uint32_t*)malloc(npt*sizeof(uint32_t));
	for (int i = 0; i < npt; i++) {
		dx_in[i] = dx(pt0, pts[i]);
		idx_in[i] = i;
		old2new[i] = -1;//not yet set
	}
	radix11sort_f32(dx_in, dx_out, idx_in, idx_out, npt);
	free(dx_in);
	free(idx_in);
	int nnew = 0; //number of unique vertices
	float tol = 0.00001; //tolerance: accept two vertices as identical if they are nearer
	for (int i=0;i<npt;i++) {
		if (old2new[idx_out[i]] >= 0)
			continue; //already assigned
		float dx0 = dx_out[i];
		vec3d pt0 = pts[idx_out[i]];
		int j = i;
		while ((j < npt) && ((dx_out[j] - dx0) < tol)) {
			if (dx(pt0, pts[idx_out[j]]) < tol) {
				old2new[idx_out[j]] = nnew;
			}
			j++;
		}
		nnew++;
	}
	free(dx_out);
	free(idx_out);
	if (npt == nnew) {
		if (verbose)
			printf("Unify vertices found no shared vertices\n");
		free(old2new);
		return npt;
	}
	for (int i=0;i<ntri;i++) { //remap face indices
		tris[i].x = old2new[tris[i].x];
		tris[i].y = old2new[tris[i].y];
		tris[i].z = old2new[tris[i].z];
	}
	vec3d *oldpts = (vec3d *) malloc(npt * sizeof(vec3d));
	for (int i=0;i<npt;i++)
		oldpts[i] = pts[i];
	free(*inpt);
	*inpt = (vec3d *) malloc(nnew * sizeof(vec3d));
	pts = *inpt;
	for (int i=0;i<npt;i++)
		pts[old2new[i]] = oldpts[i];
	free(oldpts);
	free(old2new);
	if (verbose)
		printf("vertex welding %d -> %d: %ld ms\n", npt, nnew, timediff(startTime, clockMsec()));
	return nnew;
}

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209290e-07F // float
//#define DBL_EPSILON 2.2204460492503131e-16 // double
#endif

int remove_degenerate_triangles(vec3d *pts, vec3i **intris, int ntri, bool verbose) {
	//reduces the number of triangles, number of vertices unchanged
	double startTime = clockMsec();
	vec3i *tris = *intris;
	int *isdegenerate = (int *) malloc(ntri * sizeof(int));
	int ndegenerate = 0;
	for (int i=0;i<ntri;i++) {
		//sorted lengths a ≥ b ≥ c
		isdegenerate[i] = 0;
		double l = dx(pts[tris[i].x], pts[tris[i].y]);
		double m = dx(pts[tris[i].x], pts[tris[i].z]);
		double n = dx(pts[tris[i].y], pts[tris[i].z]);
		double c = fmin(fmin(l, m), n);
		double a = fmax(fmax(l, m), n);
		double b = l+m+n-a-c;
		if ((c-(a-b)) <= 0.0) {
			isdegenerate[i] = 1;
			ndegenerate ++;
			continue;
		}
		#define REQUIRE_SIGNIFICANT_AREA
		#ifdef REQUIRE_SIGNIFICANT_AREA
		//use Heron’s Formula to eliminate triangles of tiny area
		// see Kahan: Miscalculating Area and Angles of a Needle-like Triangle
		// https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
		double area4 = 0.25 * sqrt( (a+(b+c)) * (c-(a-b)) * (c+(a-b)) * (a+(b-c)) );
		if (area4 < FLT_EPSILON) {
			isdegenerate[i] = 1;
			ndegenerate ++;
		}
		#endif //REQUIRE_SIGNIFICANT_AREA
	}
	if (ndegenerate == 0) {
		free(isdegenerate);
		return ntri;
	}
	int newtri = ntri - ndegenerate;
	vec3i *oldtris = (vec3i *) malloc(ntri * sizeof(vec3i));
	for (int i=0;i<ntri;i++)
		oldtris[i] = tris[i];
	free(*intris);
	*intris = (vec3i *) malloc(newtri * sizeof(vec3i));
	tris = *intris;
	int j = 0;
	for (int i=0;i<ntri;i++) {
		if (isdegenerate[i])
			continue;
		tris[j] = oldtris[i];
		j++;
	}
	free(oldtris);
	free(isdegenerate);
	if (verbose)
		printf("remove degenerates %d -> %d: %ld ms\n", ntri, newtri, timediff(startTime, clockMsec()));
	return newtri;
}

int quick_smooth(float * img, int nx, int ny, int nz) {
	if ((nx < 5) || (ny < 5) || (nz < 5))
		return EXIT_FAILURE;
	int nvox = nx * ny * nz;
	float *tmp = (float *) malloc(nvox * sizeof(float));
	#define kwid 2
	#define k0 0.45
	#define k1 0.225
	#define k2 0.05
	int nxy = nx * ny;
	int nxy2 = nxy * 2;
	int nx2 = nx * 2;
	//smooth column direction
	memcpy(tmp, img, nvox * sizeof(float)); //dst,src,n
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			int zy = (y * nx) + (z * nxy);
			for (int x = kwid; x < (nx-kwid); x++) {
				int v = zy + x;
				tmp[v] = (img[v-2]*k2)+(img[v-1]*k1)+(img[v]*k0)+(img[v+1]*k1)+(img[v+2]*k2);
			}
		}
	}
	//smooth row direction:
	memcpy(img, tmp, nvox * sizeof(float)); //dst,src,n
	for (int z = 0; z < nz; z++) {
		for (int x = 0; x < nx; x++) {
			int xz = x + (z * nxy);
			for (int y = kwid; y < (ny - kwid); y++) {
				int v = xz + (y * nx);
				tmp[v] = (img[v-nx2]*k2)+(img[v-nx]*k1)+(img[v]*k0)+(img[v+nx]*k1)+(img[v+nx2]*k2);
			}
		}
	}
	memcpy(img, tmp, nvox * sizeof(float)); //dst,src,n
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			int yx = (y * nx) + x;
			for (int z = kwid; z < (nz-kwid); z++) {
				int v = yx + (z * nxy);
				img[v] = (tmp[v-nxy2]*k2)+(tmp[v-nxy]*k1)+(tmp[v]*k0)+(tmp[v+nxy]*k1)+(tmp[v+nxy2]*k2);
			}
		}
	}
	free(tmp);
	return EXIT_SUCCESS;
}

void dilate(float * img, size_t dim[3], bool is26) {
	int nx = dim[0];
	int ny = dim[1];
	int nz = dim[2];
	int nxy = nx * ny;
	int nvox = nx * ny * nz;
	uint8_t* mask = (uint8_t *) malloc(nvox * sizeof(uint8_t));
	memset(mask, 0, nvox * sizeof(uint8_t));
	int numk = 6;
	if (is26)
		numk = 26;
	int32_t *k = (int32_t *)malloc(numk * sizeof(int32_t)); //queue with untested seed
	if (is26) {
		int j = 0;
		for (int z = -1; z <= 1; z++)
			for (int y = -1; y <= 1; y++)
				for (int x = -1; x <= 1; x++) {
					if ((x == 0) && (y == 0) && (z == 0)) continue;
					k[j] = x + (y * nx) + (z * nx * ny);
					j++;
				} //for x
	} else { //if 26 neighbors else 6..
		k[0] = nx * ny; //up
		k[1] = -k[0]; //down
		k[2] = nx; //anterior
		k[3] = -k[2]; //posterior
		k[4] = 1; //left
		k[5] = -1;
	}
	for (int z = 1; z < (nz - 1); z++) {
		for (int y = 1; y < (ny - 1); y++) {
			size_t iyz = + (z * nxy) + (y * nx);
			for (int x = 1; x < (nx - 1); x++) {
				size_t vx = iyz + x;
				for (int n = 1; n < numk; n++) {
					if (img[vx + k[n]] > 0)
						mask[vx] = 1;
				} //check all neighbors
			} //x
		} //y
	} //z
	for (int v = 1; v < nvox; v++)
		if (mask[v] > 0)
			img[v] = 1;
	free(mask);
	free(k);
}

double clockMsec() { //return milliseconds since midnight
	struct timespec _t;
	clock_gettime(CLOCK_MONOTONIC, &_t);
	return _t.tv_sec*1000.0 + (_t.tv_nsec/1.0e6);
}

long timediff(double startTimeMsec, double endTimeMsec) {
	return round(endTimeMsec - startTimeMsec);
}

int meshify(float * img, size_t dim[3], float isolevel, vec3i **t, vec3d **p, int *nt, int *np, int preSmooth, bool onlyLargest, bool fillBubbles, bool verbose) {
// img: input volume
// hdr: nifti header
// isolevel: air/surface threshold
// t: triangle indices e.g. [0,1,3] indicates triangle composed of vertices 0,1,3
// p: 3D points, aka vertices
// nt: number of triangles, aka faces
// np: number of points
// preSmooth: Gaussian blur to soften image
	int NX = dim[0];
	int NY = dim[1];
	int NZ = dim[2];
	if (preSmooth) {
		double startTime = clockMsec();
		quick_smooth(img, NX, NY, NZ);
		if (verbose)
			printf("pre-smooth: %ld ms\n", timediff(startTime, clockMsec()));
	}
	int nvox = NX*NY*NZ;
	float mx = img[0];
	float mn = mx;
	for (int i=0; i< nvox; i++) {
		mx = fmaxf(mx, img[i]);
		mn = fminf(mn, img[i]);
	}
	if (mn == mx) {
		printf("Error: No variability in image intensity.\n");
		return EXIT_FAILURE;
	}
	if ((isolevel <= mn) || (isolevel > mx)) {
		isolevel = 0.5 * (mn + mx);
		printf("Suggested isolevel out of range. Intensity range %g..%g, setting isolevel to %g\n", mn, mx, isolevel);
	}
	if (verbose)
		printf("intensity range %g..%g, isolevel %g\n", mn, mx, isolevel);
	int NXY = NX * NY;
	//fill
	if ((onlyLargest) || (fillBubbles)) {
		double startTime = clockMsec();
		float* mask = (float *) malloc(nvox * sizeof(float));
		size_t dim[3] = {(size_t)NX, (size_t)NY, (size_t)NZ};
		memset(mask, 0, nvox * sizeof(float));
		for (int i=0;i<nvox;i++)
			if (img[i] >= isolevel)
				mask[i] = 1;
		bwlabel(mask, 18, dim, onlyLargest, fillBubbles);
		if (fillBubbles) {
			for (int i=0;i<nvox;i++)
				if (mask[i] != 0)
					img[i] = fmax(img[i], isolevel);
		}
		if (onlyLargest) {
			dilate(mask, dim, true); //expand by one voxel to preserve subvoxel edges
			for (int i=0;i<nvox;i++)
				if (mask[i] == 0)
					img[i] = mn;
		}
		free(mask);
		if (verbose)
			printf("voxel clustering (largest cluster, bubbles): %ld ms\n", timediff(startTime, clockMsec()));
	}
	//edge darken
	float edgeMax = 0.75 * (mn + isolevel);
	int vx = 0;
	int lo[3] = {NX,NY,NZ};
	int hi[3] = {0,0,0};
	for (int z=0;z<NZ;z++) //darken edges
		for (int y=0;y<NY;y++)
			for (int x=0;x<NX;x++) {
				if (img[vx] >= isolevel) {
					lo[0] = MIN(x, lo[0]);
					lo[1] = MIN(y, lo[1]);
					lo[2] = MIN(z, lo[2]);
					hi[0] = MAX(x, hi[0]);
					hi[1] = MAX(y, hi[1]);
					hi[2] = MAX(z, hi[2]);
				}
				if ((x == 0) || (y == 0) || (z == 0) || (x == (NX-1)) || (y == (NY-1)) || (z == (NZ-1)) )
					img[vx] = fminf(edgeMax, img[vx]);
				vx++;
			}
	//printf("Bounding box for bright voxels: %d..%d %d..%d %d..%d\n", lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
	for (int i=0;i<3;i++) {
		lo[i] = MAX(lo[i] - 1, 0);
		hi[i] = MIN(hi[i] + 2, dim[i]);
	}
	// Polygonise the grid
	double startTimeMC = clockMsec();
	GRIDCELL grid;
	int ptsCapacity = 65536 * 16;
	vec3d *pts = (vec3d *) malloc(ptsCapacity * sizeof(vec3d));
	int npt = 0;
	for (int z=lo[2];z<hi[2]-1;z++) {
		for (int y=lo[1];y<hi[1]-1;y++) {
			int zy = (z * NXY) + (y * NX);
			for (int x=lo[0];x<hi[0]-1;x++) {
				grid.val[0] = img[zy+x];
				grid.val[1] = img[zy+x+1];
				grid.val[2] = img[zy+x+1+NX];
				grid.val[3] = img[zy+x+NX];
				grid.val[4] = img[zy+x+NXY];
				grid.val[5] = img[zy+x+1+NXY];
				grid.val[6] = img[zy+x+1+NX+NXY];
				grid.val[7] = img[zy+x+NX+NXY];
				grid.p000 = setVec3d(x  ,y  ,z  );
				if ((npt + 30) > ptsCapacity) {
					ptsCapacity = round (ptsCapacity * 1.5);
					pts = (vec3d *)realloc(pts,ptsCapacity*sizeof(vec3d));
				}
				vec3d *ptsn = &pts[npt];
				npt += PolygoniseCube(grid,isolevel,ptsn);
			}
		}
	}
	if (npt < 3) {
		free(pts);
		return EXIT_FAILURE;
	}
	pts = (vec3d *)realloc(pts,npt*sizeof(vec3d));
	int ntri = npt / 3;
	vec3i *tris = (vec3i *) malloc(ntri * sizeof(vec3i));
	int j = 0;
	for (int i=0;i<ntri;i++) {
		tris[i].x = j;
		j++;
		tris[i].y = j;
		j++;
		tris[i].z = j;
		j++;
	}
	if (verbose)
		printf("marching cubes: %ld ms\n", timediff(startTimeMC, clockMsec()));
	npt = unify_vertices(&pts, tris, ntri, verbose);
	if (npt < 3) return EXIT_FAILURE;
	ntri = remove_degenerate_triangles(pts, &tris, ntri, verbose);
	*t = tris;
	*p = pts;
	*nt = ntri;
	*np = npt;
	return EXIT_SUCCESS;
}

bool littleEndianPlatform () {
	uint32_t value = 1;
	return (*((char *) &value) == 1);
}

void swap_4bytes( size_t n , void *ar ) { // 4 bytes at a time
	size_t ii ;
	unsigned char * cp0 = (unsigned char *)ar, * cp1, * cp2 ;
	unsigned char tval ;
	for( ii=0 ; ii < n ; ii++ ){
		cp1 = cp0; cp2 = cp0+3;
		tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
		cp1++;  cp2--;
		tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
		cp0 += 4;
	}
	return ;
}

typedef struct {
	float x,y,z;
} vec3s; //single precision (float32)

vec3s vec3d2vec4s(vec3d v) {
	return (vec3s){.x = (float)v.x, .y = (float)v.y, .z = (float)v.z};
} // convert float64 to float32

int save_freesurfer(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
//FreeSurfer Triangle Surface Binary Format http://www.grahamwideman.com/gw/brain/fs/surfacefileformats.htm
	uint8_t magic[3] = {0xFF, 0xFF, 0xFE};
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fwrite(magic, 3, 1, fp);
	time_t t = time(NULL);
	char s[128] = "";
	struct tm *tm = localtime(&t);
	strftime(s, sizeof(s), "created by niimath on %c\n\n", tm);
	fwrite(s, strlen(s), 1, fp);
	int32_t VertexCount = npt;
	int32_t FaceCount = ntri;
	if ( &littleEndianPlatform) {
		swap_4bytes(1, &VertexCount);
		swap_4bytes(1, &FaceCount);
	}
	fwrite(&VertexCount, sizeof(int32_t), 1, fp);
	fwrite(&FaceCount, sizeof(int32_t), 1, fp);
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++) //double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (&littleEndianPlatform)
		swap_4bytes(npt * 3, pts32);
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	if (&littleEndianPlatform) {
		vec3i *trisSwap = (vec3i *) malloc(ntri * sizeof(vec3i));
		for (int i = 0; i < ntri; i++)
			trisSwap[i] = tris[i];
		swap_4bytes(ntri * 3, trisSwap);
		fwrite(trisSwap, ntri * sizeof(vec3i), 1, fp);
		free(trisSwap);
	} else
		fwrite(tris, ntri * sizeof(vec3i), 1, fp);
	fclose(fp);
	return EXIT_SUCCESS;
}

#ifdef HAVE_ZLIB
#ifdef HAVE_JSON

enum TZipMethod {zmZlib, zmGzip, zmBase64, zmLzip, zmLzma, zmLz4, zmLz4hc};

int zmat_run(const size_t inputsize, unsigned char *inputstr, size_t *outputsize, unsigned char **outputbuf, const int zipid, int *ret, const int iscompress){
	z_stream zs;
	size_t buflen[2]={0};
	*outputbuf=NULL;
	zs.zalloc = Z_NULL;
	zs.zfree = Z_NULL;
	zs.opaque = Z_NULL;
	if(inputsize==0)
		return -1;
	if(iscompress){
		/** perform compression or encoding   */
		if(zipid==zmBase64){
			/** base64 encoding  */
			*outputbuf=base64_encode((const unsigned char*)inputstr, inputsize, outputsize);
		}else if(zipid==zmZlib){
			/** zlib (.zip) or gzip (.gz) compression  */
			if(deflateInit(&zs,  (iscompress>0) ? Z_DEFAULT_COMPRESSION : (-iscompress)) != Z_OK)
				return -2;
			buflen[0] =deflateBound(&zs,inputsize);
			*outputbuf=(unsigned char *)malloc(buflen[0]);
			zs.avail_in = inputsize; /* size of input, string + terminator*/
			zs.next_in = (Bytef *)inputstr; /* input char array*/
			zs.avail_out = buflen[0]; /* size of output*/
			zs.next_out =  (Bytef *)(*outputbuf); /*(Bytef *)(); // output char array*/
			*ret=deflate(&zs, Z_FINISH);
			*outputsize=zs.total_out;
			if(*ret!=Z_STREAM_END && *ret!=Z_OK)
				return -3;
			deflateEnd(&zs);
		}else{
			return -7;
		}
	}else{
		/** perform decompression or decoding */
		if(zipid==zmBase64){
			/** base64 decoding  */
			*outputbuf=base64_decode((const unsigned char*)inputstr, inputsize, outputsize);
		}else if(zipid==zmZlib){
			/** zlib (.zip) or gzip (.gz) decompression */
			int count=1;
			if(zipid==zmZlib)
				if(inflateInit(&zs) != Z_OK)
					return -2;
			buflen[0] =inputsize*20;
			*outputbuf=(unsigned char *)malloc(buflen[0]);
			zs.avail_in = inputsize; /* size of input, string + terminator*/
			zs.next_in =inputstr; /* input char array*/
			zs.avail_out = buflen[0]; /* size of output*/
			zs.next_out =  (Bytef *)(*outputbuf); /*(Bytef *)(); // output char array*/
			while((*ret=inflate(&zs, Z_SYNC_FLUSH))!=Z_STREAM_END && count<=10){
				*outputbuf=(unsigned char *)realloc(*outputbuf, (buflen[0]<<count));
				zs.next_out =  (Bytef *)(*outputbuf+(buflen[0]<<(count-1)));
				zs.avail_out = (buflen[0]<<(count-1)); /* size of output*/
				count++;
			}
			*outputsize=zs.total_out;

			if(*ret!=Z_STREAM_END && *ret!=Z_OK)
				return -3;
			inflateEnd(&zs);
		}else{
			return -7;
		}
	}
	return 0;
}

int save_jmsh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	FILE *fp;
	cJSON *root=NULL, *hdr=NULL, *node=NULL, *face=NULL;
	char *jsonstr=NULL;
	int dim[2]={0,3}, len[2]={1,0};
	size_t compressedbytes, totalbytes;
	unsigned char *compressed=NULL, *buf=NULL;
	int ret=0, status=0;
	root=cJSON_CreateObject();
	cJSON_AddItemToObject(root,  "_DataInfo_", hdr = cJSON_CreateObject());
	cJSON_AddStringToObject(hdr, "JMeshVersion", "0.5");
	cJSON_AddStringToObject(hdr, "Comment", "Created by nii2mesh");
	cJSON_AddItemToObject(root,  "MeshVertex3", node = cJSON_CreateObject());
	cJSON_AddStringToObject(node,"_ArrayType_","double");
	dim[0]=npt;
	cJSON_AddItemToObject(node,  "_ArraySize_",cJSON_CreateIntArray(dim,2));
	cJSON_AddStringToObject(node,"_ArrayZipType_","zlib");
	len[1]=dim[0]*dim[1];
	cJSON_AddItemToObject(node,  "_ArrayZipSize_",cJSON_CreateIntArray(len,2));
	totalbytes=dim[0]*dim[1]*sizeof(pts[0].x);
	ret=zmat_run(totalbytes, (unsigned char *)&(pts[0].x), &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
	if(!ret){
		 ret=zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status,1);
		 cJSON_AddStringToObject(node,  "_ArrayZipData_",(char *)buf);
	}
	if(compressed){
		free(compressed);
		compressed=NULL;
	}
	if(buf){
		free(buf);
		buf=NULL;
	}
	cJSON_AddItemToObject(root,  "MeshTri3", face = cJSON_CreateObject());
	cJSON_AddStringToObject(face,"_ArrayType_","uint32");
	dim[0]=ntri;
	cJSON_AddItemToObject(face,  "_ArraySize_",cJSON_CreateIntArray(dim,2));
	cJSON_AddStringToObject(face,"_ArrayZipType_","zlib");
	len[1]=dim[0]*dim[1];
	cJSON_AddItemToObject(face,  "_ArrayZipSize_",cJSON_CreateIntArray(len,2));
	totalbytes=dim[0]*dim[1]*sizeof(tris[0].x);
	ret=zmat_run(totalbytes, (unsigned char *)&(tris[0].x), &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
	if(!ret){
		ret=zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status,1);
		cJSON_AddStringToObject(face,  "_ArrayZipData_",(char *)buf);
	}
	if(compressed)
		free(compressed);
	if(buf)
		free(buf);
	jsonstr=cJSON_Print(root);
	if(jsonstr==NULL)
		return EXIT_FAILURE;
	fp=fopen(fnm,"wt");
	if(fp==NULL)
		return EXIT_FAILURE;
	fprintf(fp,"%s\n",jsonstr);
	fclose(fp);
	if(jsonstr)
		free(jsonstr);
	if(root)
		cJSON_Delete(root);
	return EXIT_SUCCESS;
}
#endif //HAVE_JSON
#endif //HAVE_ZLIB

int save_json(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	FILE *fp = fopen(fnm,"w");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp,"{\n");
	fprintf(fp,"\t\"_DataInfo_\":{\n\t\t\"JMeshVersion\":\"0.5\",\n\t\t\"Comment\":\"Created by nii2mesh\"\n\t},\n");
	fprintf(fp,"\t\"MeshVertex3\":[\n");
	for (int i=0;i<npt;i++)
		fprintf(fp, "[%g,\t%g,\t%g],\n", pts[i].x, pts[i].y,pts[i].z);
	fprintf(fp,"\t],\n\t\"MeshTri3\":[\n");
	for (int i=0;i<ntri;i++)
		fprintf(fp, "[%d,\t%d,\t%d],\n", tris[i].x+1, tris[i].y+1, tris[i].z+1);
	fprintf(fp,"\t]\n}\n");
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_mz3(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz) {
//https://github.com/neurolabusc/surf-ice/tree/master/mz3
	struct __attribute__((__packed__)) mz3hdr {
		uint16_t SIGNATURE, ATTR;
		uint32_t NFACE, NVERT, NSKIP;
	};
	struct mz3hdr h;
	h.SIGNATURE = 0x5A4D;
	h.ATTR = 3;//isFACE +1 isVERT +2
	h.NFACE = ntri;
	h.NVERT = npt;
	h.NSKIP = 0;
	if (! &littleEndianPlatform)
		swap_4bytes(3, &h.NFACE);
	FILE *fp;
	#ifdef HAVE_ZLIB
	gzFile fgz;
	if (isGz) {
		fgz = gzopen(fnm, "w");
		if (! fgz)
			return EXIT_FAILURE;
		gzwrite(fgz, &h, sizeof(struct mz3hdr));
	} else
	#endif
	{
		fp = fopen(fnm,"wb");
		if (fp == NULL)
			return EXIT_FAILURE;
		fwrite(&h, sizeof(struct mz3hdr), 1, fp);
	}
	if (! &littleEndianPlatform) {
		vec3i *trisSwap = (vec3i *) malloc(ntri * sizeof(vec3i));
		for (int i = 0; i < ntri; i++)
			trisSwap[i] = tris[i];
		swap_4bytes(ntri * 3, trisSwap);
		#ifdef HAVE_ZLIB
		if (isGz)
			gzwrite(fgz, trisSwap, ntri * sizeof(vec3i));
		else
		#else
			fwrite(trisSwap, ntri * sizeof(vec3i), 1, fp);
		#endif
		free(trisSwap);
	} else {
		#ifdef HAVE_ZLIB
		if (isGz)
			gzwrite(fgz, tris, ntri * sizeof(vec3i));
		else
		#endif
			fwrite(tris, ntri * sizeof(vec3i), 1, fp);
	}
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++) //double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (! &littleEndianPlatform)
		swap_4bytes(npt * 3, pts32);
	#ifdef HAVE_ZLIB
	if (isGz) {
		gzwrite(fgz, pts32, npt * sizeof(vec3s));
		gzclose(fgz);
	} else
	#endif
	{
		fwrite(pts32, npt * sizeof(vec3s), 1, fp);
		fclose(fp);
	}
	free(pts32);
	return EXIT_SUCCESS;
}

int save_off(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	FILE *fp = fopen(fnm,"w");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp,"OFF\n%d\t%d\t0\n",npt,ntri);
	for (int i=0;i<npt;i++)
		fprintf(fp, "%g %g %g\n", pts[i].x, pts[i].y,pts[i].z);
	for (int i=0;i<ntri;i++)
		fprintf(fp, "%d %d %d\n", tris[i].x+1, tris[i].y+1, tris[i].z+1);
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_obj(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	FILE *fp = fopen(fnm,"w");
	if (fp == NULL)
		return EXIT_FAILURE;
	for (int i=0;i<npt;i++)
		fprintf(fp, "v %g %g %g\n", pts[i].x, pts[i].y,pts[i].z);
	for (int i=0;i<ntri;i++)
		fprintf(fp, "f %d %d %d\n", tris[i].x+1, tris[i].y+1, tris[i].z+1);
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_stl(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	//binary STL http://paulbourke.net/dataformats/stl/
	//n.b. like other tools, ignores formal restriction that all adjacent facets must share two common vertices.
	//n.b. does not write normal
	typedef struct  __attribute__((__packed__)) {
		vec3s norm, pts[3];
		uint16_t spacer;
	} tfacet;
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	uint8_t hdr[80] = { 0 };
	fwrite(hdr, 80, 1, fp);
	int32_t nf = ntri;
	fwrite(&nf, sizeof(int32_t), 1, fp);
	tfacet *facets = (tfacet *) malloc(ntri * sizeof(tfacet));
	vec3s n0 = (vec3s){.x = 0.0, .y = 0.0, .z = 0.0};
	for (int i = 0; i < ntri; i++) { //double->single precision
		facets[i].norm = n0;
		facets[i].pts[0] = vec3d2vec4s(pts[tris[i].x]);
		facets[i].pts[1] = vec3d2vec4s(pts[tris[i].y]);
		facets[i].pts[2] = vec3d2vec4s(pts[tris[i].z]);
		facets[i].spacer = 0;
	}
	fwrite(facets, ntri * sizeof(tfacet), 1, fp);
	free(facets);
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_ply(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	typedef struct  __attribute__((__packed__)) {
		uint8_t n;
		int32_t x,y,z;
	} vec1b3i;
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("ply\n",fp);
	if (&littleEndianPlatform)
		fputs("format binary_little_endian 1.0\n",fp);
	else
		fputs("format binary_big_endian 1.0\n",fp);
	fputs("comment niimath\n",fp);
	char vpts[80];
	sprintf(vpts, "element vertex %d\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	fputs("property float x\n",fp);
	fputs("property float y\n",fp);
	fputs("property float z\n",fp);
	char vfc[80];
	sprintf(vfc, "element face %d\n", ntri);
	fwrite(vfc, strlen(vfc), 1, fp);
	fputs("property list uchar int vertex_indices\n",fp);
	fputs("end_header\n",fp);
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++) { //double->single precision
		pts32[i].x = pts[i].x;
		pts32[i].y = pts[i].y;
		pts32[i].z = pts[i].z;
	}
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	vec1b3i *tris4 = (vec1b3i *) malloc(ntri * sizeof(vec1b3i));
	for (int i = 0; i < ntri; i++) { //double->single precision
		tris4[i].n = 3;
		tris4[i].x = tris[i].x;
		tris4[i].y = tris[i].y;
		tris4[i].z = tris[i].z;
	}
	fwrite(tris4, ntri * sizeof(vec1b3i), 1, fp);
	free(tris4);
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_gii(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz){
	//https://www.nitrc.org/projects/gifti/
	//https://stackoverflow.com/questions/342409/how-do-i-base64-encode-decode-in-c
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",fp);
	fputs("<!DOCTYPE GIFTI SYSTEM \"http://www.nitrc.org/frs/download.php/115/gifti.dtd\">\n",fp);
	fputs("<GIFTI Version=\"1.0\"  NumberOfDataArrays=\"2\">\n",fp);
	fputs("   <MetaData>\n",fp);
	fputs("      <MD>\n",fp);
	fputs("          <Name><![CDATA[nii2mesh-version]]></Name>\n",fp);
	fputs("          <Value><![CDATA[nii2mesh, 1 Jan 2022]]></Value>\n",fp);
	fputs("      </MD>\n",fp);
	fputs("   </MetaData>\n",fp);
	fputs("   <LabelTable/>\n",fp);
	fputs("   <DataArray  ArrayIndexingOrder=\"RowMajorOrder\"\n",fp);
	fputs("               DataType=\"NIFTI_TYPE_INT32\"\n",fp);
	char fc[80];
	sprintf(fc, "               Dim0=\"%d\"\n", ntri);
	fwrite(fc, strlen(fc), 1, fp);
	fputs("               Dim1=\"3\"\n",fp);
	fputs("               Dimensionality=\"2\"\n",fp);
	#ifdef HAVE_ZLIB
	if (isGz)
		fputs("               Encoding=\"GZipBase64Binary\"\n" ,fp);
	else
	#endif
		fputs("               Encoding=\"Base64Binary\"\n" ,fp);
	if (&littleEndianPlatform)
		fputs("               Endian=\"LittleEndian\"\n",fp);
	else
		fputs("               Endian=\"BigEndian\"\n",fp);
	fputs("               ExternalFileName=\"\"\n",fp);
	fputs("               ExternalFileOffset=\"\"\n",fp);
	fputs("               Intent=\"NIFTI_INTENT_TRIANGLE\">\n" ,fp);
	fputs("      <MetaData>\n",fp);
	fputs("      </MetaData>\n",fp);
	fputs("      <Data>",fp);
	size_t out_len;
	unsigned char * fcs;
	#ifdef HAVE_ZLIB
	if (isGz) {
		unsigned long srcLen = ntri * sizeof(vec3i);
		unsigned long destLen = compressBound(srcLen);
		unsigned char* ostream = (unsigned char*) malloc(destLen);
		int res = compress(ostream, &destLen,(const unsigned char *) tris, srcLen);
		if (res != Z_OK)
			printf("Compression error\n");
		fcs = base64_encode(ostream, destLen, &out_len);
		free(ostream);
	} else
	#endif
		fcs = base64_encode((const unsigned char *) tris, ntri * sizeof(vec3i), &out_len);
	fwrite(fcs, out_len, 1, fp);
	free(fcs);
	fputs("</Data>\n",fp);
	fputs("   </DataArray>\n",fp);
	fputs("   <DataArray  ArrayIndexingOrder=\"RowMajorOrder\"\n",fp);
	fputs("               DataType=\"NIFTI_TYPE_FLOAT32\"\n",fp);
	char vpts[80];
	sprintf(vpts, "               Dim0=\"%d\"\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	fputs("               Dim1=\"3\"\n",fp);
	fputs("               Dimensionality=\"2\"\n",fp);
	#ifdef HAVE_ZLIB
	if (isGz)
		fputs("               Encoding=\"GZipBase64Binary\"\n" ,fp);
	else
	#endif
		fputs("               Encoding=\"Base64Binary\"\n" ,fp);
	if (&littleEndianPlatform)
		fputs("               Endian=\"LittleEndian\"\n",fp);
	else
		fputs("               Endian=\"BigEndian\"\n",fp);
	fputs("               ExternalFileName=\"\"\n",fp);
	fputs("               ExternalFileOffset=\"\"\n",fp);
	fputs("               Intent=\"NIFTI_INTENT_POINTSET\">\n",fp);
	fputs("      <MetaData>\n",fp);
	fputs("      </MetaData>\n",fp);
	fputs("      <CoordinateSystemTransformMatrix>\n",fp);
	fputs("         <DataSpace><![CDATA[NIFTI_XFORM_UNKNOWN]]></DataSpace>\n",fp);
	fputs("         <TransformedSpace><![CDATA[NIFTI_XFORM_UNKNOWN]]></TransformedSpace>\n",fp);
	fputs("         <MatrixData>1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 </MatrixData>\n",fp);
	fputs("      </CoordinateSystemTransformMatrix>\n",fp);
	fputs("      <Data>",fp);
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++)//double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	unsigned char * vts;
	#ifdef HAVE_ZLIB
	if (isGz) {
		unsigned long srcLen = npt * sizeof(vec3s);
		unsigned long destLen = compressBound(srcLen);
		unsigned char* ostream = (unsigned char*) malloc(destLen);
		int res = compress(ostream, &destLen,(const unsigned char *) pts32, srcLen);
		if (res != Z_OK)
			printf("Compression error\n");
		vts = base64_encode(ostream, destLen, &out_len);
		free(ostream);
	} else
	#endif
		vts = base64_encode((const unsigned char *) pts32, npt * sizeof(vec3s), &out_len);
	free(pts32);
	fwrite(vts, out_len, 1, fp);
	free(vts);
	fputs("</Data>\n",fp);
	fputs("   </DataArray>\n",fp);
	fputs("</GIFTI>\n",fp);
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_vtk(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt){
	typedef struct {
		uint32_t n, x,y,z;
	} vec4i;
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("# vtk DataFile Version 3.0\n",fp);
	fputs("this file was written using niimath\n",fp);
	fputs("BINARY\n",fp);
	fputs("DATASET POLYDATA\n",fp);
	char vpts[80];
	sprintf(vpts, "POINTS %d float\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++)//double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (&littleEndianPlatform)
		swap_4bytes(3*npt, pts32);
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	char vfac[80];
	sprintf(vfac, "POLYGONS %d %d\n", ntri, ntri * 4);
	fwrite(vfac, strlen(vfac), 1, fp);
	vec4i *tris4 = (vec4i *) malloc(ntri * sizeof(vec4i));
	for (int i = 0; i < ntri; i++) { //double->single precision
		tris4[i].n = 3;
		tris4[i].x = tris[i].x;
		tris4[i].y = tris[i].y;
		tris4[i].z = tris[i].z;
	}
	if (&littleEndianPlatform)
		swap_4bytes(4*ntri, tris4);
	fwrite(tris4, ntri * sizeof(vec4i), 1, fp);
	free(tris4);
	fclose(fp);
	return EXIT_SUCCESS;
}

void strip_ext(char *fname){
	char *end = fname + strlen(fname);
	while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
		--end;
	}
	if ((end > fname && *end == '.') &&
		(*(end - 1) != '\\' && *(end - 1) != '/')) {
		*end = '\0';
	}
}

int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz){
	char basenm[768], ext[768] = "";
	strcpy(basenm, fnm);
	strip_ext(basenm); // ~/file.nii -> ~/file
	if (strlen(fnm) > strlen(basenm))
		strcpy(ext, fnm + strlen(basenm));
	if (strstr(ext, ".gii"))
		return save_gii(fnm, tris, pts, ntri, npt, isGz);
	else if ((strstr(ext, ".inflated")) || (strstr(ext, ".pial")))
		return save_freesurfer(fnm, tris, pts, ntri, npt);
#ifdef HAVE_ZLIB
#ifdef HAVE_JSON
	else if (strstr(ext, ".jmsh"))
		return save_jmsh(fnm, tris, pts, ntri, npt);
#endif //HAVE_JSON
#endif //HAVE_ZLIB
	else if (strstr(ext, ".json"))
		return save_json(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".mz3"))
		return save_mz3(fnm, tris, pts, ntri, npt, isGz);
	else if (strstr(ext, ".off"))
		return save_off(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".obj"))
		return save_obj(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".ply"))
		return save_ply(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".stl"))
		return save_stl(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".vtk"))
		return save_vtk(fnm, tris, pts, ntri, npt);
	strcpy(basenm, fnm);
	strcat(basenm, ".obj");
	return save_obj(basenm, tris, pts, ntri, npt);
}

double sform(vec3d p, float srow[4]) {
	return (p.x * srow[0])+(p.y * srow[1])+(p.z * srow[2])+srow[3];
}

void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]){
	for (int i = 0; i < np; i++) {
		vec3d v = p[i];
		p[i].x = sform(v, srow_x);
		p[i].y = sform(v, srow_y);
		p[i].z = sform(v, srow_z);
	}
	//detect determinant
	vec3d p0;
	p0.x = srow_x[0] + srow_x[1] + srow_x[2];
	p0.y = srow_y[0] + srow_y[1] + srow_y[2];
	p0.z = srow_z[0] + srow_z[1] + srow_z[2];
	float det = p0.x * p0.y * p0.z;
	if (det >= 0.0) return; //positive volume
	//negative volume: we need to reverse the triangle winding
	for (int i = 0; i < nt; i++) {
		vec3i f = t[i];
		t[i].x = f.y;
		t[i].y = f.x;
	}
}
