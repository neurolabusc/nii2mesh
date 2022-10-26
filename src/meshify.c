#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _MSC_VER

#else
 #include <unistd.h>
#endif
#include <time.h>
#ifdef HAVE_ZLIB
	#include <zlib.h>
	#ifdef HAVE_JSON
		#include "cJSON.h"
		#define FIND_JSON_KEY(id,parent,fallback,val) \
			((tmp=cJSON_GetObjectItem(root,id))==0 ? fallback : tmp->val)
	#endif
#endif
#include "meshify.h"
#include "base64.h" //required for GIfTI
#include "bwlabel.h"
#include "radixsort.h"
#include "meshtypes.h"
#ifdef USE_CLASSIC_CUBES
	#include "oldcubes.h"
#else
	#include "MarchingCubes.h"
#endif


#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

double sqr(double x) {
	return x * x;
}

double dx(vec3d p0, vec3d p1) {
	return sqrt( sqr(p0.x - p1.x) + sqr(p0.y - p1.y) + sqr(p0.z - p1.z));
}

int unify_vertices(vec3d **inpt, vec3i *tris, int npt, int ntri, bool verbose) {
	//"vertex welding": reduces the number of vertices, number of faces unchanged
	double startTime = clockMsec();
	vec3d *pts = *inpt;
	int *old2new = (int *)malloc(npt * sizeof(int));
	vec3d pt0 = pts[0];
	float* dx_in = (float*)malloc(npt*sizeof(float));
	float* dx_out = (float*)malloc(npt*sizeof(float));
	uint32_t* idx_in = (uint32_t*)malloc(npt*sizeof(uint32_t));
	uint32_t* idx_out = (uint32_t*)malloc(npt*sizeof(uint32_t));
	for (int i = 0; i < npt; i++) {
		//printf("%d %d\n", i, npt);
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
		printf("remove degenerate triangles %d -> %d: %ld ms\n", ntri, newtri, timediff(startTime, clockMsec()));
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
#ifdef _MSC_VER
	clock_t t = clock();
	return (double)((double)t) / (CLOCKS_PER_SEC / 1000.0);
#else
	#ifdef __MINGW32__
		time_t seconds_since_midnight = time(NULL) % 86400;
		return seconds_since_midnight;
	#else
		struct timespec _t;
		clock_gettime(CLOCK_MONOTONIC, &_t);
		return _t.tv_sec*1000.0 + (_t.tv_nsec/1.0e6);
	#endif
#endif
}

long timediff(double startTimeMsec, double endTimeMsec) {
	return round(endTimeMsec - startTimeMsec);
}

int meshify(float * img, size_t dim[3], int originalMC, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, int preSmooth, bool onlyLargest, bool fillBubbles, bool verbose) {
// img: input volume
// hdr: nifti header
// isolevel: air/surface threshold
// t: triangle indices e.g. [0,1,3] indicates triangle composed of vertices 0,1,3
// p: 3D points, aka vertices
// nt: number of triangles, aka faces
// np: number of points
	int NX = dim[0];
	int NY = dim[1];
	int NZ = dim[2];
	int nvox = NX*NY*NZ;
	// preSmooth: Gaussian blur to soften image
	if (preSmooth) {
		double startTime = clockMsec();
		quick_smooth(img, NX, NY, NZ);
		if (verbose)
			printf("pre-smooth: %ld ms\n", timediff(startTime, clockMsec()));
	}
	//determine image intensity range - ensure isolvel will detect edge
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
	//(optional) fill bubbles and only extract largest contiguous object
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
		hi[i] = MIN(hi[i] + 2, dim[i]-1);
	}
	double startTimeMC = clockMsec();
	vec3d *pts = NULL;
	vec3i *tris = NULL;
	int ntri;
	int npt;
	if (marchingCubes(img, dim, lo, hi, originalMC, isolevel, &pts, &tris, &npt, &ntri) != EXIT_SUCCESS)
		return EXIT_FAILURE;
	if (verbose)
		printf("marching cubes (%dx%dx%d): %ld ms\n", NX, NY, NZ, timediff(startTimeMC, clockMsec()));
	npt = unify_vertices(&pts, tris, npt, ntri, verbose);
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

void swap_8bytes( size_t n , void *ar )    // 8 bytes at a time
{
    size_t ii ;
    unsigned char * cp0 = (unsigned char *)ar, * cp1, * cp2 ;
    unsigned char tval ;
    for( ii=0 ; ii < n ; ii++ ){
        cp1 = cp0; cp2 = cp0+7;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp1++;  cp2--;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp1++;  cp2--;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp1++;  cp2--;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp0 += 8;
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
	if (littleEndianPlatform()) {
		swap_4bytes(1, &VertexCount);
		swap_4bytes(1, &FaceCount);
	}
	fwrite(&VertexCount, sizeof(int32_t), 1, fp);
	fwrite(&FaceCount, sizeof(int32_t), 1, fp);
	vec3s *pts32 = (vec3s *) malloc(npt * sizeof(vec3s));
	for (int i = 0; i < npt; i++) //double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (littleEndianPlatform())
		swap_4bytes(npt * 3, pts32);
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	if (littleEndianPlatform()) {
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
		}else if(zipid==zmZlib || zipid==zmGzip){
			/** zlib (.zip) or gzip (.gz) compression  */
			if(zipid==zmZlib){
				if(deflateInit(&zs,  (iscompress>0) ? Z_DEFAULT_COMPRESSION : (-iscompress)) != Z_OK)
					return -2;
			}else{
				if(deflateInit2(&zs, (iscompress>0) ? Z_DEFAULT_COMPRESSION : (-iscompress), Z_DEFLATED, 15|16, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK)
					return -2;
			}
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
		}else if(zipid==zmZlib || zipid==zmGzip){
			/** zlib (.zip) or gzip (.gz) decompression */
			int count=1;
			if(zipid==zmZlib){
				if(inflateInit(&zs) != Z_OK)
					return -2;
			}else{
				if(inflateInit2(&zs, 15|32) != Z_OK)
					return -2;
			}
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

#ifdef HAVE_JSON

const char *zipformat[]={"zlib","gzip","base64","lzip","lzma","lz4","lz4hc",""};

int key_lookup(char *origkey, const char *table[]){
	int i=0;
	while(table[i] && table[i][0]!='\0'){
		if(strcmp(origkey,table[i])==0)
			return i;
		i++;
	}
	return -1;
}

// decoding JData ND array construct {"_ArraySize_":,"_ArrayType_":"_ArrayZipType_":"_ArrayZipSize_":"_ArrayZipData_":}

int  jdata_decode(void **vol, unsigned short *ndim, unsigned short *dims, int maxdim, char **type, cJSON *obj){
	int ret=0;
	cJSON * ztype=NULL;
	cJSON * vsize=cJSON_GetObjectItem(obj,"_ArraySize_");
	cJSON * vtype=cJSON_GetObjectItem(obj,"_ArrayType_");
	cJSON * vdata=cJSON_GetObjectItem(obj,"_ArrayData_");
	if(!vdata){
		ztype=cJSON_GetObjectItem(obj,"_ArrayZipType_");
		vdata=cJSON_GetObjectItem(obj,"_ArrayZipData_");
	}
	if(vtype)
		*type=vtype->valuestring;
	if(vdata){
		if(vsize){
			cJSON *tmp=vsize->child;
			*ndim=cJSON_GetArraySize(vsize);
			for(int i=0;i<MIN(maxdim,*ndim);i++){
				dims[i]=tmp->valueint;
				tmp=tmp->next;
			}
		}
		if(ztype){
			size_t len, newlen;
			int status=0;
			char *buf=NULL;
			int zipid=key_lookup((char *)(ztype->valuestring),zipformat);
			if(zipid<0 ||zipid>zmBase64)
				 return -1;
			if(zipid==zmBase64)
				return zmat_run(strlen(vdata->valuestring), (unsigned char *)vdata->valuestring, &len, (unsigned char **)vol, zmBase64, &status, 0);
			else
				ret=zmat_run(strlen(vdata->valuestring), (unsigned char *)vdata->valuestring, &len, (unsigned char **)&buf, zmBase64, &status, 0);
			if(!ret && vsize){
				if(*vol)
					free(*vol);
				ret=zmat_run(len, (unsigned char *)buf, &newlen, (unsigned char **)(vol), zipid, &status, 0);
			}
			if(buf)
				free(buf);
		}else
			return -1;
        }else
		 return -1;
	return ret;
}

void read_vec4float(float *vec4, cJSON *vec){
	for(int i=0;i<4;i++){
		vec4[i]=vec->valuedouble;
		vec=vec->next;
	}
}

void array3d_row2col(float **vol, unsigned short *dim){
	unsigned int x,y,z;
	unsigned int dimxy,dimyz;
	float *newvol=NULL;

	if(*vol==NULL || dim[0]==0 || dim[1]==0 || dim[2]==0){
		return;
	}
	newvol=(float *)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
	dimxy=dim[0]*dim[1];
	dimyz=dim[1]*dim[2];
	for(x=0;x<dim[0];x++)
		for(y=0;y<dim[1];y++)
			for(z=0;z<dim[2];z++){
				newvol[z*dimxy+y*dim[0]+x]=(*vol)[x*dimyz+y*dim[2]+z];
			}
	free(*vol);
	*vol=newvol;
}

// parsing a JSON/JNIfTI-encoded jnii volume file

float * load_jnii(const char *fnm, nifti_1_header * hdr) {
	char *jbuf;
	int len;
	float * img32=NULL;
	cJSON *root, *jniihead, *jniidata, *tmp;

	// reading JNIfTI/JSON file to a string
	FILE *fp=fopen(fnm,"rb");
	if(!fp){
		printf("Unable to open the specified file %s\n", fnm);
		return NULL;
	}
	fseek (fp, 0, SEEK_END);
	len=ftell(fp)+1;
	jbuf=(char *)malloc(len);
	rewind(fp);
	if(fread(jbuf,len-1,1,fp)!=1)
		return NULL;
	jbuf[len-1]='\0';
	fclose(fp);

	// parse JNIfTI/JSON file
	root = cJSON_Parse(jbuf);

	if(!root){ // if not a valid JSON file, print error and return
		char *ptrold, *ptr=(char*)cJSON_GetErrorPtr();
		if(ptr) ptrold=strstr(jbuf,ptr);
		if(fp!=NULL) fclose(fp);
		if(ptr && ptrold){
			char *offs=(ptrold-jbuf>=50) ? ptrold-50 : jbuf;
			while(offs<ptrold){
				fprintf(stderr,"%c",*offs);
				offs++;
			}
			fprintf(stderr,"<error>%.50s\n",ptrold);
		}
		if(fp!=NULL)
			free(jbuf);
		return NULL;
	}
	free(jbuf);

	jniihead = cJSON_GetObjectItem(root,"NIFTIHeader");
	jniidata = cJSON_GetObjectItem(root,"NIFTIData");
	memset(hdr, 0, sizeof(nifti_1_header));
	if(jniihead){
		hdr->sizeof_hdr=FIND_JSON_KEY("NIIHeaderSize",jniihead,348,valueint);
		hdr->scl_slope=FIND_JSON_KEY("ScaleSlope",jniihead,1.0,valuedouble);
		hdr->scl_inter=FIND_JSON_KEY("ScaleOffset",jniihead,0.0,valuedouble);
		tmp=cJSON_GetObjectItem(jniihead,"Affine");
		if(tmp && cJSON_IsArray(tmp) && cJSON_GetArraySize(tmp)==3){
			cJSON *row=tmp->child;
			cJSON *elem=row->child;
			read_vec4float(hdr->srow_x, elem);
			elem=row->next->child;
			read_vec4float(hdr->srow_y, elem);
			elem=row->next->next->child;
			read_vec4float(hdr->srow_z, elem);
		}else{
			hdr->srow_x[0]=1.f;
			hdr->srow_y[1]=1.f;
			hdr->srow_z[2]=1.f;
		}
	}else
		return NULL;

	if(jniidata){
		char *type=NULL;
		void *imgRaw=NULL;

		tmp=cJSON_GetObjectItem(jniidata, "_ArrayType_");
		if(!tmp)
			return NULL;

		if(jdata_decode((void **)&imgRaw, (unsigned short *)hdr->dim, (unsigned short *)hdr->dim+1, 3, &type, jniidata)!=0){
			if(imgRaw)
				free(imgRaw);
			return NULL;
		}
		int nvox = hdr->dim[1]*hdr->dim[2]*hdr->dim[3];
		img32 = (float *) malloc(nvox*sizeof(float));
		if(strcmp(type,"uint8")==0){
			for(int i=0; i< nvox; i++)
				img32[i]=((unsigned char *)imgRaw)[i];
		}else if(strcmp(type,"int8")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((char *)imgRaw)[i];
		}else if(strcmp(type,"uint16")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((unsigned short *)imgRaw)[i];
		}else if(strcmp(type,"int16")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((short *)imgRaw)[i];
		}else if(strcmp(type,"uint32")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((unsigned int *)imgRaw)[i];
		}else if(strcmp(type,"int32")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((int *)imgRaw)[i];
		}else if(strcmp(type,"uint64")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((long long *)imgRaw)[i];
		}else if(strcmp(type,"int64")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((unsigned long long *)imgRaw)[i];
		}else if(strcmp(type,"double")==0){
			for(int i=0; i<nvox; i++)
				img32[i]=((double *)imgRaw)[i];
		}else if(strcmp(type,"single")==0){
			memcpy(img32, imgRaw, nvox*sizeof(float));
		}else{
			if(imgRaw)
				free(imgRaw);
			return NULL;
		}
		if(imgRaw)
			free(imgRaw);

		tmp=cJSON_GetObjectItem(jniidata, "_ArrayOrder_");
		if(!tmp || (cJSON_IsString(tmp) && ((tmp->valuestring)[0]=='r' || (tmp->valuestring)[0]=='R')))
			array3d_row2col(&img32, (unsigned short *)hdr->dim+1);
	}else
		return NULL;
	cJSON_Delete(root);
	return img32;
}
#endif // HAVE_JSON

int save_jmsh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isdouble){
	FILE *fp;
	const char *pyparsers[]={"https://pypi.org/project/jdata","https://pypi.org/project/bjdata"};
	const char *jsparsers[]={"https://www.npmjs.com/package/jda","https://www.npmjs.com/package/bjd"};
	const char *cparsers[]={"https://github.com/DaveGamble/cJSON","https://github.com/NeuroJSON/ubj"};
	size_t compressedbytes, totalbytes;
	unsigned char *compressed=NULL, *buf=NULL;
	int ret=0, status=0;

	float *floatpts=NULL;
	if(!isdouble){
		floatpts=(float *)malloc(npt*3*sizeof(float));
		for(int i=0; i<npt; i++){
			floatpts[i*3]=pts[i].x;
			floatpts[i*3+1]=pts[i].y;
			floatpts[i*3+2]=pts[i].z;
		}
	}

	fp=fopen(fnm,"wt");
	if(fp==NULL)
		return EXIT_FAILURE;

	fprintf(fp, "{\n\t\"_DataInfo_\":{\n");
	fprintf(fp, "\t\t\"JMeshVersion\":\"0.5\",\n");
	fprintf(fp, "\t\t\"Comment\":\"Created by nii2mesh with NeuroJSON JMesh format (https://neurojson.org)\",\n");
	fprintf(fp, "\t\t\"AnnotationFormat\":\"https://neurojson.org/jmesh/draft1\",\n");
	fprintf(fp, "\t\t\"SerialFormat\":\"https://json.org\",\n");
	fprintf(fp, "\t\t\"Parser\":{\n");
	fprintf(fp, "\t\t\t\"Python\":[\"%s\",\t\"%s\"],\n",pyparsers[0],pyparsers[1]);
	fprintf(fp, "\t\t\t\"MATLAB\":\"https://github.com/NeuroJSON/jsonlab\",\n");
	fprintf(fp, "\t\t\t\"JavaScript\":[\"%s\",\t\"%s\"],\n",jsparsers[0],jsparsers[1]);
	fprintf(fp, "\t\t\t\"CPP\":\"https://github.com/NeuroJSON/json\",\n");
	fprintf(fp, "\t\t\t\"C\":[\"%s\",\t\"%s\"]\n",cparsers[0],cparsers[1]);
	fprintf(fp, "\t\t}\n\t},\n");
	fprintf(fp, "\t\"MeshVertex3\":{\n");
	fprintf(fp, "\t\t\"_ArrayType_\":\"%s\",\n", (isdouble ? "double" : "single"));
	fprintf(fp, "\t\t\"_ArraySize_\":[%d, 3],\n", npt);
	fprintf(fp, "\t\t\"_ArrayZipType_\":\"zlib\",\n");
	fprintf(fp, "\t\t\"_ArrayZipSize_\":[1, %d],\n", npt*3);

	totalbytes=npt*3*(isdouble ? sizeof(pts[0].x) : sizeof(float));
	ret=zmat_run(totalbytes, (isdouble ? (unsigned char *)&(pts[0].x) : (unsigned char *)floatpts), &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
	if(!ret){
		ret=zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status,1);
		fprintf(fp, "\t\t\"_ArrayZipData_\":\"");
		fwrite(buf, 1, totalbytes, fp);
		fprintf(fp, "\"\n");
	}
	if(compressed){
		free(compressed);
		compressed=NULL;
	}
	if(buf){
		free(buf);
		buf=NULL;
	}
	fprintf(fp, "\t},\n");
	fprintf(fp, "\t\"MeshTri3\":{\n");
	fprintf(fp, "\t\t\"_ArrayType_\":\"uint32\",\n");
	fprintf(fp, "\t\t\"_ArraySize_\":[%d, 3],\n", ntri);
	fprintf(fp, "\t\t\"_ArrayZipType_\":\"zlib\",\n");
	fprintf(fp, "\t\t\"_ArrayZipSize_\":[1, %d],\n", ntri*3);

	totalbytes=ntri*3*sizeof(tris[0].x);
	unsigned int *val=(unsigned int *)malloc(totalbytes);
	memcpy(val,&(tris[0].x),totalbytes);
	for(int i=0;i<ntri*3;i++)
		val[i]++;
	ret=zmat_run(totalbytes, (unsigned char *)val, &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
	free(val);
	if(!ret){
		ret=zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status,1);
		fprintf(fp, "\t\t\"_ArrayZipData_\":\"");
		fwrite(buf, 1, totalbytes, fp);
		fprintf(fp, "\"\n");
	}
	if(compressed)
		free(compressed);
	if(buf)
		free(buf);
	if(floatpts)
		free(floatpts);

	fprintf(fp, "\t}\n}\n");
	fclose(fp);
	return EXIT_SUCCESS;
}

void write_ubjsonint(int len, int *dat, FILE *fp){
	if (!littleEndianPlatform())
		swap_4bytes(len, dat);
	fwrite(dat,len,4,fp);
}

void write_ubjsonfloat(int len, float *dat, FILE *fp){
	if (!littleEndianPlatform())
		swap_4bytes(len, dat);
	fwrite(dat,len,4,fp);
}

void write_ubjsondouble(int len, double *dat, FILE *fp){
	if (!littleEndianPlatform())
		swap_8bytes(len, dat);
	fwrite(dat,len,8,fp);
}

int save_bmsh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz, bool isdouble){
	int markerlen=0;
	const char *output[]={
	"{",
		"N","","_DataInfo_","{",
			"N","","JMeshVersion","S","","0.5",
			"N","","Comment","S","","Created by nii2mesh with NeuroJSON Binay JMesh format (https://neurojson.org)",
			"N","","AnnotationFormat","S","","https://neurojson.org/jmesh/draft1",
			"N","","SerialFormat","S","","https://neurojson.org/bjdata/draft2",
			"N","","Parser","{",
				"N","","Python","[",
					"S","","https://pypi.org/project/jdata",
					"S","","https://pypi.org/project/bjdata",
				"]",
				"N","","MATLAB","S","","https://github.com/NeuroJSON/jsonlab",
				"N","","JavaScript","[",
					"S","","https://www.npmjs.com/package/jda",
					"S","","https://www.npmjs.com/package/bjd",
				"]",
				"N","","CPP","S","","https://github.com/NeuroJSON/json",
				"N","","C","[",
					"S","","https://github.com/DaveGamble/cJSON",
					"S","","https://github.com/NeuroJSON/ubj",
				"]",
			"}",
		"}",
		"N","","MeshVertex3","?1",
			"N","","_ArrayType_","S","",(isdouble)?"double":"single",
			"N","","_ArraySize_","[$l#U\x2","?2",
			"N","","_ArrayZipType_","S","","zlib",
			"N","","_ArrayZipSize_","[$l#U\x1","?3",
			"N","","_ArrayZipData_","S","","?4",
		"}",
		"N","","MeshTri3","{",
			"N","","_ArrayType_","S","","uint32",
			"N","","_ArraySize_","[$l#U\x2","?5",
			"N","","_ArrayZipType_","S","","zlib",
			"N","","_ArrayZipSize_","[$l#U\x1","?6",
			"N","","_ArrayZipData_","S","","?7",
		"}",
	"}"
	};

	float *floatpts=NULL;
	if(!isdouble){
		floatpts=(float *)malloc(npt*3*sizeof(float));
		for(int i=0; i<npt; i++){
			floatpts[i*3]=pts[i].x;
			floatpts[i*3+1]=pts[i].y;
			floatpts[i*3+2]=pts[i].z;
		}
	}
	unsigned int *newtris=(unsigned int *)malloc(ntri*3*sizeof(unsigned int));
	memcpy(newtris,&(tris[0].x),ntri*3*sizeof(unsigned int));
	for(int i=0;i<ntri*3;i++)
		newtris[i]++;

	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;

	markerlen=sizeof(output)/sizeof(char*);

	for(int i=0;i<markerlen;i++){
		int slen=strlen(output[i]);
		if(slen>0 && output[i][0]!='?'){
			if(!(slen==1 && output[i][0]=='N'))
				fwrite(output[i],1,slen,fp);
			if(slen==1 && (output[i][0]=='N' || output[i][0]=='S') && i+2<markerlen && output[i+1][0]=='\0' && output[i+2][0]!='?'){
				unsigned int keylen=strlen(output[i+2]);
				if(keylen<256){
					unsigned char keylenbyte=keylen;
					fputc('U',fp);
					fwrite(&keylenbyte,1,sizeof(keylenbyte),fp);
				}else{
					fputc('l',fp);
					write_ubjsonint(1,(int *)(&keylen),fp);
				}
			}
		}else{
			if(slen>0){
				int slotid=0;
				if(sscanf(output[i],"\?%d",&slotid)==1 && slotid>0){
					unsigned char *compressed=NULL;
					size_t compressedbytes, totalbytes;
					int dim[2]={0,3};
					int ret=0, status=0;
					switch(slotid){
						case 1: {
							if(isGz){
								fputc('{',fp);
							}else{ // write node and face data in BJData strongly-typed ND array construct
								int val[2]={0};
								val[0]=npt;
								val[1]=3;

								fwrite("[$", 1, 2, fp);
								fputc(isdouble ? 'D' : 'd', fp);
								fwrite("#[$l#U\x02", 1, 7, fp);
								write_ubjsonint(2,val,fp);
								if(isdouble)
									write_ubjsondouble(npt*3, &(pts[0].x), fp);
								else
									write_ubjsonfloat(npt*3, floatpts, fp);

								val[0]=ntri;
								fwrite("U\x08MeshTri3[$l#[$l#U\x02", 1, 20, fp);
								write_ubjsonint(2, val, fp);
								write_ubjsonint(ntri*3, (int *)newtris, fp);
								fputc('}', fp);
							}
							break;
						}
						case 2: {int val[2]; val[0]=npt; val[1]=3; write_ubjsonint(2,val,fp);break;}
						case 3: {int val=npt*3;	write_ubjsonint(1,&val,fp);break;}
						case 4:
							dim[0]=npt;

							totalbytes=dim[0]*dim[1]*(isdouble? sizeof(pts[0].x) : sizeof(float));
							ret=zmat_run(totalbytes, (isdouble ? (unsigned char *)&(pts[0].x) : (unsigned char *)floatpts) , &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
							if(!ret){
								int clen=compressedbytes;
								fputc('l',fp);
								write_ubjsonint(1,&clen,fp);
								fwrite(compressed,1,compressedbytes,fp);
							}
							if(compressed)
								free(compressed);
							break;
						case 5: {int val[2]; val[0]=ntri;val[1]=3; write_ubjsonint(2,val,fp);break;}
						case 6: {int val=ntri*3; write_ubjsonint(1,&val,fp);break;}
						case 7:
							dim[0]=ntri;
							totalbytes=dim[0]*dim[1]*sizeof(tris[0].x);
							ret=zmat_run(totalbytes, (unsigned char *)newtris, &compressedbytes, (unsigned char **)&compressed, zmZlib, &status,1);
							if(!ret){
								int clen=compressedbytes;
								fputc('l',fp);
								write_ubjsonint(1,&clen,fp);
								fwrite(compressed,1,compressedbytes,fp);
							}

							if(compressed)
								free(compressed);
							break;
					}
					if(slotid==1 && !isGz)
						break;
				}
			}
		}
	}
	fclose(fp);

	free(newtris);
	if(floatpts)
		free(floatpts);

	return EXIT_SUCCESS;
}
#endif //HAVE_ZLIB

int save_json(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isdouble){
	FILE *fp = fopen(fnm,"w");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp,"{\n");
	fprintf(fp,"\t\"_DataInfo_\":{\n\t\t\"JMeshVersion\":\"0.5\",\n\t\t\"Comment\":\"Created by nii2mesh with NeuroJSON JMesh format (https://neurojson.org)\"\n\t},\n");
	fprintf(fp,"\t\"MeshVertex3\":[\n");
        if(isdouble){
                for (int i=0;i<npt;i++)
                        fprintf(fp, "[%.16g,\t%.16g,\t%.16g],\n", pts[i].x, pts[i].y,pts[i].z);
        }else{
                for (int i=0;i<npt;i++)
                        fprintf(fp, "[%.7g,\t%.7g,\t%.7g],\n", pts[i].x, pts[i].y,pts[i].z);
        }
	fprintf(fp,"\t],\n\t\"MeshTri3\":[\n");
	for (int i=0;i<ntri;i++)
		fprintf(fp, "[%d,\t%d,\t%d],\n", tris[i].x+1, tris[i].y+1, tris[i].z+1);
	fprintf(fp,"\t]\n}\n");
	fclose(fp);
	return EXIT_SUCCESS;
}

int save_mz3(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz) {
//https://github.com/neurolabusc/surf-ice/tree/master/mz3
	#ifdef _MSC_VER
#pragma pack(2)
	struct mz3hdr {
		uint16_t SIGNATURE, ATTR;
		uint32_t NFACE, NVERT, NSKIP;
	};
#pragma pack()	
	#else
	struct __attribute__((__packed__)) mz3hdr {
		uint16_t SIGNATURE, ATTR;
		uint32_t NFACE, NVERT, NSKIP;
	};
	#endif
	struct mz3hdr h;
	h.SIGNATURE = 0x5A4D;
	h.ATTR = 3;//isFACE +1 isVERT +2
	h.NFACE = ntri;
	h.NVERT = npt;
	h.NSKIP = 0;
	if (!littleEndianPlatform())
		swap_4bytes(3, &h.NFACE);
	FILE *fp=NULL;
	#ifdef HAVE_ZLIB
	gzFile fgz=NULL;
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
	if (!littleEndianPlatform()) {
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
	if (!littleEndianPlatform())
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
	#ifdef _MSC_VER
#pragma pack(2)
	typedef struct  {
		vec3s norm, pts[3];
		uint16_t spacer;
	} tfacet;
#pragma pack()	
	#else
	typedef struct  __attribute__((__packed__)) {
		vec3s norm, pts[3];
		uint16_t spacer;
	} tfacet;
	#endif
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
	#ifdef _MSC_VER
#pragma pack(1)
	typedef struct  {
		uint8_t n;
		int32_t x,y,z;
	} vec1b3i;
#pragma pack()	
	#else
	typedef struct  __attribute__((__packed__)) {
		uint8_t n;
		int32_t x,y,z;
	} vec1b3i;
	#endif
	FILE *fp = fopen(fnm,"wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("ply\n",fp);
	if (littleEndianPlatform())
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
	if (littleEndianPlatform())
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
	if (littleEndianPlatform())
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
	if (littleEndianPlatform())
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
	if (littleEndianPlatform())
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

int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz, bool isdouble){
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
	else if (strstr(ext, ".jmsh"))
		return save_jmsh(fnm, tris, pts, ntri, npt, isdouble);
	else if (strstr(ext, ".bmsh"))
		return save_bmsh(fnm, tris, pts, ntri, npt, isGz, isdouble);
#endif //HAVE_ZLIB
	else if (strstr(ext, ".json"))
		return save_json(fnm, tris, pts, ntri, npt, isdouble);
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
