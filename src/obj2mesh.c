// gcc -O3 -DNII2MESH -DHAVE_ZLIB MarchingCubes.c obj2mesh.c isolevel.c meshify.c quadric.c base64.c bwlabel.c radixsort.c -o obj2mesh -lz -lm
// obj2mesh dragon.obj reduced.obj
//
//For debugging ()
//gcc  -O1 -g -fsanitize=address -fno-omit-frame-pointer   -DNII2MESH -DHAVE_ZLIB MarchingCubes.c obj2mesh.c isolevel.c meshify.c quadric.c base64.c bwlabel.c radixsort.c -o obj2mesh -lz -lm
// obj2mesh -t 3749 -v 2 sphere.obj smoother.obj

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "meshify.h"
#include "quadric.h"
#if defined(_OPENMP)
	#include <omp.h>
#endif
#ifdef HAVE_ZLIB
	#include <zlib.h>
#endif

void read_obj(const char* filename, vec3d **verts, vec3i **tris, int* nvert, int* ntri){
	//printf ( "Loading Objects %s ... \n",filename);
	*ntri = 0;
	*nvert = 0;
	FILE* fn;
	if(filename==NULL) return ;
	if((char)filename[0]==0) return ;
	if ((fn = fopen(filename, "rb")) == NULL) {
		printf ( "File %s not found!\n" ,filename );
		return;
	}
	int tCapacity = 65536 * 16;
	*tris = (vec3i *) malloc(tCapacity * sizeof(vec3i));
	vec3i *ts = *tris;
	int vCapacity = 65536 * 16;
	*verts = (vec3d *) malloc(vCapacity * sizeof(vec3d));
	vec3d *vs = *verts;
	char line[1000];
	memset ( line,0,1000 );
	int vertex_cnt = 0;
	while(fgets( line, 1000, fn ) != NULL) {
		vec3d v;
		if (( line[0] == 'v' ) && ( line[1] == ' ' )) {
			if(sscanf(line,"v %lf %lf %lf", &v.x, &v.y, &v.z)==3) {
				if (*nvert >= vCapacity) {
					vCapacity = round (vCapacity * 1.5);
					*verts = (vec3d *)realloc(*verts,vCapacity*sizeof(vec3d));
					vs = *verts;
				}
				vs[*nvert] = v;
				*nvert = *nvert + 1;
			}
		}
		vec3i t;
		if ( line[0] == 'f' ) {
			if(sscanf(line,"f %d %d %d", &t.x, &t.y, &t.z)==3) {
				if (*ntri >= tCapacity) {
					tCapacity = round (tCapacity * 1.5);
					*tris = (vec3i *)realloc(*tris,tCapacity*sizeof(vec3i));
					ts = *tris;
				}
				t = (vec3i){ .x = t.x-1, .y = t.y-1, .z = t.z-1 };
				ts[*ntri] = t;
				*ntri = *ntri + 1;
			}
		}
	}
	fclose(fn);
	*verts = (vec3d *)realloc(*verts,*nvert*sizeof(vec3d));
	*tris = (vec3i *)realloc(*tris,*ntri*sizeof(vec3i));
} // read_obj()

int main(int argc,char **argv) {
	double agressiveness = 7.0;
	double reduceFraction = 0.5;
	int smoothIter = 0;
	int quality = 1;
	int target_count = 0;
	bool verbose = false;
	if (argc < 3) {
		printf("Converts a NIfTI voxelwise image to triangulated mesh.\n");
		printf("Usage: %s [options] niftiname meshname\n",argv[0]);
		printf("Options\n");
		printf("    -q v    quality (0=lossy, 1=lossy then losslesss, default %d)\n", quality);
		printf("    -r v    reduction factor (default %g)\n", reduceFraction);
		printf("    -s v    smoothing iterations (default %d)\n", smoothIter);
		printf("    -t v    target count (unused by default overrides reduction factor)\n");
		printf("    -v v    verbose (0=silent, 1=verbose, default %d)\n", verbose);
		printf("mesh extension sets format (.gii, .mz3, .obj, .ply, .pial, .stl, .vtk)\n");
		printf("Example: '%s dragon.obj small.obj'\n",argv[0]);
		printf("Example: '%s -v 1 -r 0.1 dragon.obj smaller.gii'\n",argv[0]);
		printf("Example: '%s -s 30 dragon.obj smoother.gii'\n",argv[0]);
		exit(-1);
	}
	for (int i=1;i<argc;i++) {
		if (strcmp(argv[i],"-q") == 0)
			quality = atoi(argv[i+1]);
		if (strcmp(argv[i],"-r") == 0)
			reduceFraction = atof(argv[i+1]);
		if (strcmp(argv[i],"-s") == 0)
			smoothIter = atoi(argv[i+1]);
		if (strcmp(argv[i],"-t") == 0)
			target_count = atoi(argv[i+1]);
		if (strcmp(argv[i],"-v") == 0)
			verbose = atoi(argv[i+1]);
	}
	vec3d *verts = NULL;
	vec3i *tris = NULL;
	int nvert;
	int ntri;
	read_obj(argv[argc-2], &verts, &tris, &nvert, &ntri);
	if ((nvert < 4) || (ntri < 2))
		exit(1);
	int startTri = ntri;
	int startVert = nvert;
	if (target_count <= 0)
		target_count = round((float)ntri * reduceFraction);
	double startTime = clockMsec();
	if (smoothIter > 0) {
		laplacian_smoothHC(verts, tris, nvert, ntri, 0.1, 0.5, smoothIter, true);
		if (verbose)
			printf("%d smoothing iterations: %ld ms\n", smoothIter, timediff(startTime, clockMsec()));
		startTime = clockMsec();
	}
	quadric_simplify_mesh(&verts, &tris, &nvert, &ntri, target_count, agressiveness, verbose, (quality > 0));
	if (verbose)
		printf("simplify: %ld ms\n", timediff(startTime, clockMsec()));
	if (verbose)
		printf("simplify vertices %d->%d triangles: %d->%d (r = %g)\n", startVert, nvert, startTri, ntri, (float)ntri / (float) startTri);
	save_mesh(argv[argc-1], tris, verts, ntri, nvert, true);
	free(tris);
	free(verts);
}
