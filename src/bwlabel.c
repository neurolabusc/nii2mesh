#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _MSC_VER

#else 
 #include <unistd.h>
#endif
#include <stdint.h>
#include "bwlabel.h"

#define printfx(...) fprintf(stderr, __VA_ARGS__)

//Jesper Andersson has acknowledged that this port of spm_bwlabel.c may be released using the BSD 2-Clause license

//Copyright 2021 Jesper Andersson
//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//for usage, see https://en.wikibooks.org/wiki/SPM/How-to#How_to_remove_clusters_under_a_certain_size_in_a_binary_mask.3F

/****************************************************************
 **
 ** Set of routines implementing a 2D or 3D connected component
 ** labelling algorithm. Its interface is modelled on bwlabel
 ** (which is a routine in the image processing toolbox) and
 ** takes as input a binary image and (optionally) a connectednes
 ** criterion (6, 18 or 26, in 2D 6 will correspond to 4 and 18
 ** and 26 to 8). It will output an image/volume where each
 ** connected component will have a unique label.
 **
 ** The implementation is not recursive (i.e. will no crash for
 ** large connected components) and is loosely based on
 ** Thurfjell et al. 1992, A new three-dimensional connected
 ** components labeling algorithm with simultaneous object
 ** feature extraction capability. CVGIP: Graphical Models
 ** and Image Processing 54(4):357-364.
 **
 ***************************************************************/

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

static void fill_tratab(uint32_t  *tt, uint32_t  *nabo, uint32_t  nr_set)  {
/*
*tt   Translation table
*nabo   Set of neighbours
nr_set   Number of neighbours in nabo
*/
   int           i = 0, j = 0, cntr = 0;
   uint32_t  tn[9];
   uint32_t  ltn = UINT_MAX;

   /*
   Find smallest terminal number in neighbourhood
   */

   for (i=0; i<nr_set; i++)
   {
      j = nabo[i];
      cntr=0;
      while (tt[j-1] != j)
      {
         j = tt[j-1];
         cntr++;
         if (cntr>100) {printfx("\nOoh no!!"); break;}
      }
      tn[i] = j;
      ltn = MIN(ltn,j);
   }
   /*
   Replace all terminal numbers in neighbourhood by the smallest one
   */
   for (i=0; i<nr_set; i++)
   {
      tt[tn[i]-1] = ltn;
   }

   return;
}

#define idx(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

static uint32_t check_previous_slice(uint32_t  *il,     /* Initial labelling map */
                                  uint32_t  r,       /* row */
                                  uint32_t  c,       /* column */
                                  uint32_t  sl,      /* slice */
                                  size_t        dim[3],  /* dimensions of il */
                                  uint32_t  conn,    /* Connectivity criterion */
                                  uint32_t  *tt)     /* Translation table */
//                                  uint32_t  ttn)     /* Size of translation table */
{
   uint32_t l=0;
   uint32_t nabo[9];
   uint32_t nr_set = 0;

   if (!sl) return(0);
   if (conn >= 6)
   {
      if ((l = il[idx(r,c,sl-1,dim)])) {nabo[nr_set++] = l;}
   }
   if (conn >= 18)
   {
      if (r) {if ((l = il[idx(r-1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c) {if ((l = il[idx(r,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r < dim[0]-1) {if ((l = il[idx(r+1,c,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (c < dim[1]-1) {if ((l = il[idx(r,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }
   if (conn == 26)
   {
      if (r && c) {if ((l = il[idx(r-1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && c) {if ((l = il[idx(r+1,c-1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if (r && (c < dim[1]-1)) {if ((l = il[idx(r-1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
      if ((r < dim[0]-1) && (c < dim[1]-1)) {if ((l = il[idx(r+1,c+1,sl-1,dim)])) {nabo[nr_set++] = l;}}
   }

   if (nr_set)
   {
      fill_tratab(tt,/*ttn,*/nabo,nr_set);
      return(nabo[0]);
   }
   else {return(0);}
}

static void * mxRealloc(void *oldArray, size_t oldBytes, size_t newBytes) {
   // https://octave.org/doxygen/3.8/df/d4e/mex_8cc_source.html
   //reallocate memory, preserve previous bytes
   if (newBytes <= 0) {
      free(oldArray);
      return NULL;
   }
   void *newArray = (void *)malloc(newBytes);
   memset(newArray, 0, newBytes);
   if (oldBytes > 0) {
     //void * memcpy ( void * destination, const void * source, size_t num );
     oldBytes = MIN(oldBytes, newBytes);
     //void * memcpy ( void * destination, const void * source, size_t num );
     memcpy(newArray, oldArray, oldBytes);
     free(oldArray);
   }
   return newArray;
}

/* do_initial_labelling */

static uint32_t do_initial_labelling(uint8_t        *bw,   /* Binary map */
                                  size_t        *dim,  /* Dimensions of bw */
                                  uint32_t  conn,  /* Connectivity criterion */
                                  uint32_t  *il,   /* Initially labelled map */
                                  uint32_t  **tt)  /* Translation table */
{
   uint32_t  i = 0, j = 0;
   uint32_t  nabo[8];
   uint32_t  label = 1;
   uint32_t  nr_set = 0;
   uint32_t  l = 0;
   int32_t       sl, r, c;
   uint32_t  ttn = 1000;
   *tt = (uint32_t *)malloc(ttn * sizeof(uint32_t));
   memset(*tt, 0, ttn * sizeof(uint32_t));
   for (sl=0; sl<dim[2]; sl++)
   {
      for (c=0; c<dim[1]; c++)
      {
         for (r=0; r<dim[0]; r++)
         {
            nr_set = 0;
            if (bw[idx(r,c,sl,dim)])
            {
               nabo[0] = check_previous_slice(il,r,c,sl,dim,conn,*tt /*,ttn*/);
               if (nabo[0]) {nr_set += 1;}
               /*
                  For six(surface)-connectivity
               */
               if (conn >= 6)
               {
                  if (r)
                  {
                     if ((l = il[idx(r-1,c,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c)
                  {
                     if ((l = il[idx(r,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               /*
                  For 18(edge)-connectivity
                  N.B. In current slice no difference to 26.
               */
               if (conn >= 18)
               {
                  if (c && r)
                  {
                     if ((l = il[idx(r-1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
                  if (c && (r < dim[0]-1))
                  {
                     if ((l = il[idx(r+1,c-1,sl,dim)])) {nabo[nr_set++] = l;}
                  }
               }
               if (nr_set)
               {
                  il[idx(r,c,sl,dim)] = nabo[0];
                  fill_tratab(*tt,/*ttn,*/nabo,nr_set);
               }
               else
               {
                  il[idx(r,c,sl,dim)] = label;
                  if (label >= ttn) {ttn += 1000; *tt = (uint32_t*)mxRealloc(*tt, (ttn - 1000)*sizeof(uint32_t), ttn*sizeof(uint32_t));}
                  (*tt)[label-1] = label;
                  label++;
               }
            }
         }
      }
   }

   /*
      Finalise translation table
   */

   for (i=0; i<(label-1); i++)
   {
      j = i;
      while ((*tt)[j] != j+1)
      {
         j = (*tt)[j]-1;
      }
      (*tt)[i] = j+1;
   }
   return(label-1);
}

/* translate_labels */

static int translate_labels(uint32_t  *il,     /* Map of initial labels. */
                        size_t        dim[3],  /* Dimensions of il. */
                        uint32_t  *tt,     /* Translation table. */
                        uint32_t  ttn,     /* Size of translation table. */
                        uint32_t        *l)      /* Final map of labels. */
{
   int            n=0;
   int            i=0;
   int   ml=0;
   int         cl = 0;
   n = dim[0]*dim[1]*dim[2];
   for (i=0; i<ttn; i++) {ml = MAX(ml,tt[i]);}
   uint32_t *fl = (uint32_t *)malloc(ml * sizeof(uint32_t));
   memset(fl, 0, ml * sizeof(uint32_t));
   for (i=0; i<n; i++)
   {
      if (il[i])
      {
         if (!fl[tt[il[i]-1]-1])
         {
            cl += 1;
            fl[tt[il[i]-1]-1] = cl;
         }
         l[i] = fl[tt[il[i]-1]-1];
      }
   }
   free(fl);
   return(cl);
}

static void fillh(uint32_t* imgBin, size_t dim[3], int is26, int nLabels) {
  //aka nifti_fillh
  //all given binary image, interior 0 voxels set to 1
  int nx = dim[0];
  int ny = dim[1];
  int nz = dim[2];
  if ((nx < 3) || (ny < 3) || (nz < 3) || (nLabels < 1))
    return;
  int nvox = dim[0] * dim[1] * dim[2];
  //uint8_t *vxv = (uint8_t *)malloc(nvox * sizeof(uint8_t));
  //memset(vxv, 0, nvox * sizeof(uint8_t));
  //set up kernel to search for neighbors. Since we already included sides, we do not worry about A<->P and L<->R wrap
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
  //https://en.wikipedia.org/wiki/Flood_fill
  int32_t *q = (int32_t *)malloc(nvox * sizeof(int32_t)); //queue with untested seed
  uint8_t *vxs = (uint8_t *)malloc(nvox * sizeof(uint8_t));
  for (int label = 1; label <= nLabels; label++) {
    for (size_t i = 0; i < nvox; i++)
      vxs[i] = (imgBin[i] == label);
    int qlo = 0;
    int qhi = -1; //ints always signed in C!
    //load edges
    size_t i = 0;
    for (int z = 0; z < nz; z++) {
      int zedge = 0;
      if ((z == 0) || (z == (nz - 1)))
        zedge = 1;
      for (int y = 0; y < ny; y++) {
        int yedge = 0;
        if ((y == 0) || (y == (ny - 1)))
          yedge = 1;
        for (int x = 0; x < nx; x++) {
          if ((vxs[i] == 0) && (zedge || yedge || (x == 0) || (x == (nx - 1)))) { //found new seed
            vxs[i] = 1; //do not find again
            qhi++;
            q[qhi] = i;
          } // new seed
          i++;
        } //for x
      } //y
    } //z
    //printf("seeds %d kernel %d\n", qhi+1, numk);
    //run a 'first in, first out' queue
    while (qhi >= qlo) {
      //retire one seed, add 0..6 new ones (fillh) or 0..26 new ones (fillh26)
      for (int j = 0; j < numk; j++) {
        int jj = q[qlo] + k[j];
        if ((jj < 0) || (jj >= nvox))
          continue;
        if (vxs[jj] != 0)
          continue;
        //add new seed;
        vxs[jj] = 1;
        qhi++;
        q[qhi] = jj;
      }
      qlo++;
    } //while qhi >= qlo: continue until all seeds tested
    for (size_t i = 0; i < nvox; i++) {
      if (vxs[i] == 0)
        imgBin[i] = label; //hidden internal voxel not found from the fill
    }
  }
  //for (size_t i = 0; i < nvox; i++)
  //  imgBin[i] = (vxs[i] == 0); //hidden internal voxel not found from the fill
  free(vxs);
  free(q);
  free(k);
}

int bwlabel(float *img, int conn, size_t dim[3], bool onlyLargest, bool fillBubbles) {
  if ((conn!=6) && (conn!=18) && (conn!=26)) {
     printfx("bwlabel: conn must be 6, 18 or 26.\n");
     return 0;
  }
  if ((dim[0] < 2) || (dim[1] < 2) || (dim[2] < 1)) {
    printfx("bwlabel: img must be 2 or 3-dimensional\n");
    return 0;
  }
  size_t nvox = dim[0] * dim[1] * dim[2];
  uint32_t *l = (uint32_t *)malloc(nvox * sizeof(uint32_t)); //output image
  memset(l, 0, nvox * sizeof(uint32_t));
  uint32_t *il = (uint32_t *)malloc(nvox * sizeof(uint32_t));
  memset(il, 0, nvox * sizeof(uint32_t));
  uint8_t *bw = (uint8_t *)malloc(nvox * sizeof(uint8_t));
  memset(bw, 0, nvox * sizeof(uint8_t));
  for (size_t i = 0; i < nvox; i++)
    if (img[i] != 0.0) bw[i] = 1;
  uint32_t  *tt = NULL;
  uint32_t ttn = do_initial_labelling(bw,dim,conn,il,&tt);
  free(bw);
  int nl = translate_labels(il,dim,tt,ttn,l);
  free(il);
  free(tt);
  if ((nl > 0) && (onlyLargest)){
    int mxL = 0;
    int mxN = 0;
    for (int j = 1; j <= nl; j++) {
      int n = 0;
      for (int i = 0; i < nvox; i++)
        if (l[i] == j)
          n++;
      if (n > mxN) {
        mxN = n;
        mxL = j;
      }
    } //for j: each label
    for (int i = 0; i < nvox; i++)
      l[i] = (l[i] == mxL);
    nl = 1;
  } //if labels found
  if (fillBubbles)
      fillh(l, dim, 1, nl);
  for (size_t i = 0; i < nvox; i++)
    img[i] = l[i];
  free(l);
  return nl;
}
