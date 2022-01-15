/*
 * Base64 encoding/decoding (RFC1341)
 * Copyright (c) 2005, Jouni Malinen <j@w1.fi>
 *
 * This software may be distributed under the terms of the BSD license.
 * See README for more details.
 */

#ifndef MESHIFY_H
#define MESHIFY_H

#include <stdbool.h>
#include "nifti1.h"

typedef struct {
	double x,y,z;
} vec3d;

typedef struct {
	int x,y,z;
} vec3i;


void strip_ext(char *fname);
int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz);
int meshify(float * img, nifti_1_header * hdr, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, int preSmooth, bool onlyLargest, bool fillBubbles, bool verbose);
void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]);
double clockMsec();
long timediff(double startTimeMsec, double endTimeMsec);

#endif /* MESHIFY_H */