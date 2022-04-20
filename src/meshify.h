#ifndef MESHIFY_H
#define MESHIFY_H

#include <stdbool.h>
#include "meshtypes.h"

#ifdef HAVE_JSON
        #include "nifti1.h"
#endif

void strip_ext(char *fname);
int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz, bool isdouble);
int meshify(float * img, size_t dim[3], int originalMC, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, int preSmooth, bool onlyLargest, bool fillBubbles, bool verbose);
void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]);
double clockMsec();
long timediff(double startTimeMsec, double endTimeMsec);
#ifdef HAVE_JSON
float * load_jnii(const char *fnm, nifti_1_header * hdr);
#endif

#endif /* MESHIFY_H */