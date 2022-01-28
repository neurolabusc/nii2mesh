#ifndef MESHIFY_H
#define MESHIFY_H

#include <stdbool.h>
#include "meshtypes.h"

void strip_ext(char *fname);
int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz);
int meshify(float * img, size_t dim[3], int originalMC, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, int preSmooth, bool onlyLargest, bool fillBubbles, bool verbose);
void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]);
double clockMsec();
long timediff(double startTimeMsec, double endTimeMsec);

#endif /* MESHIFY_H */