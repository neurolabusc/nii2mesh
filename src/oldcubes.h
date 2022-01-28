/*
Use Otsu's method to detect a reasonable threshold
*/

#ifndef OLDCUBES_H
#define OLDCUBES_H

#include "meshtypes.h"
int marchingCubes(float * img, size_t dim[3], int lo[3], int hi[3], int originalMC, float isolevel, vec3d **vs, vec3i **ts, int *nv, int *nt);

#endif /* OLDCUBES_H */