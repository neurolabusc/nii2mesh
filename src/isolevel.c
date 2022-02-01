/*
Use Otsu's method to provided reasonable isolevel thresholds
 */
#include <stdint.h>
#include "isolevel.h"
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef MAX //from Christian Gaser's TFCE example
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

static bool isnanx(float f) { //isnan disabled by gcc -Ofast and -ffinite-math-only
	return isnan(f);
}
/*
static bool isnanx(float f) { //isnan disabled by gcc -Ofast and -ffinite-math-only
//4byte IEEE: msb[31] = signbit, bits[23-30] exponent, bits[0..22] mantissa
//exponent of all 1s =   Infinity, NAN or Indeterminate
	uint32_t i = *(long *)&f;
	#define MY_LITTLE (((union { unsigned x; unsigned char c; }){1}).c)
	#ifdef MY_LITTLE
		return ((i&0x7f800000)==0x7f800000)&&(i&0x7fffff);
	#else
		return ((i&0x0000807f)==0x0000807f)&&(i&0xffff7f00);
	#endif
}*/

static int nifti_robust_range(float* img, int nvox, float *pct2, float *pct98, int ignoreZeroVoxels) {
	//https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;31f309c1.1307
	// robust range is essentially the 2nd and 98th percentiles
	// "but ensuring that the majority of the intensity range is captured, even for binary images."
	// fsl uses 1000 bins, also limits for volumes less than 100 voxels taylor.hanayik@ndcn.ox.ac.uk 20190107
	//fslstats trick -r
	// 0.000000 1129.141968
	//niimath >fslstats trick -R
	// 0.000000 2734.000000
	*pct2 = 0.0;
	*pct98 = 1.0;
	if (nvox < 1)
		return 1;
	float mn = INFINITY;
	float mx = -INFINITY;
	size_t nZero = 0;
	size_t nNan = 0;
	for (size_t i = 0; i < nvox; i++) {
		if (isnanx(img[i])) {
			nNan++;
			continue;
		}
		if (img[i] == 0.0) {
			nZero++;
			if (ignoreZeroVoxels)
				continue;
		}
		mn = fmin(img[i], mn);
		mx = fmax(img[i], mx);
	}
	if ((nZero > 0) && (mn > 0.0) && (!ignoreZeroVoxels))
		mn = 0.0;
	if (mn > mx)
		return 0; //all NaN
	if (mn == mx) {
		*pct2 = mn;
		*pct98 = mx;
		return 0;
	}
	if (!ignoreZeroVoxels)
		nZero = 0;
	nZero += nNan;
	size_t n2pct = round((nvox - nZero) * 0.02);
	if ((n2pct < 1) || (mn == mx) || ((nvox - nZero) < 100)) { //T Hanayik mentioned issue with very small volumes
		*pct2 = mn;
		*pct98 = mx;
		return 0;
	}
#define nBins 1001
	float scl = (nBins - 1) / (mx - mn);
	int hist[nBins];
	for (int i = 0; i < nBins; i++)
		hist[i] = 0;
	if (ignoreZeroVoxels) {
		for (int i = 0; i < nvox; i++) {
			if (isnanx(img[i]))
				continue;
			if (img[i] == 0.0)
				continue;
			hist[(int)round((img[i] - mn) * scl)]++;
		}
	} else {
		for (int i = 0; i < nvox; i++) {
			if (isnanx(img[i]))
				continue;
			hist[(int)round((img[i] - mn) * scl)]++;
		}
	}
	size_t n = 0;
	size_t lo = 0;
	while (n < n2pct) {
		n += hist[lo];
		//if (lo < 10)
		//	printf("%zu %zu %zu %d\n",lo, n, n2pct, ignoreZeroVoxels);
		lo++;
	}
	lo--; //remove final increment
	n = 0;
	int hi = nBins;
	while (n < n2pct) {
		hi--;
		n += hist[hi];
	}
	if (lo == hi) { //MAJORITY are not black or white
		int ok = -1;
		while (ok != 0) {
			if (lo > 0) {
				lo--;
				if (hist[lo] > 0)
					ok = 0;
			}
			if ((ok != 0) && (hi < (nBins - 1))) {
				hi++;
				if (hist[hi] > 0)
					ok = 0;
			}
			if ((lo == 0) && (hi == (nBins - 1)))
				ok = 0;
		} //while not ok
	}//if lo == hi
	*pct2 = (lo) / scl + mn;
	*pct98 = (hi) / scl + mn;
	//printf("full range %g..%g (voxels 0 or NaN =%zu) robust range %g..%g\n", mn, mx, nZero, *pct2, *pct98);
	return 0;
}

static int nii_otsu(int* H, int nBin, int mode, int *dark, int *mid, int *bright) {
//H: Histogram H[0..nBin-1] with each bin storing nuumber of pixels of this brightness
//nBin: number of bins in histogram, e.g. 256 for H[0..255]
//mode: segment and levels 1: 3/4, 2: 2/3 3: 1/2, 4: 1/3, 5: 1/4
//dark/bright/high: set to threshold
	int thresh = 0;
	*dark = 0;
	*mid = 0;
	*bright = 0;
	double Sum = 0.0;
	for (int v = 0; v < nBin; v++)
		Sum = Sum + H[v];
	if (Sum <= 0)
		return 0;
	double *P = (double*) malloc(nBin * nBin * sizeof(double));
	double *S = (double*) malloc(nBin * nBin * sizeof(double));
	P[0] = H[0];
	S[0] = H[0];
	for (int v = 1; v < nBin; v++) {
		double Prob = H[v]/Sum;
		P[v] = P[v-1]+Prob;
		S[v] = S[v-1]+(v+1)*Prob;
	}
	for (int u = 1; u < nBin; u++) {
		for (int v = u; v < nBin; v++) {
			P[(u*nBin)+v] = P[v]-P[u-1];
			S[(u*nBin)+v] = S[v]-S[u-1];
		}
	}
	//result is eq 29 from Liao
	for (int u = 0; u < nBin; u++) {
		for (int v = u; v < nBin; v++) {
			if (P[(u*nBin)+v] != 0) //avoid divide by zero errors...
				P[(u*nBin)+v] = (S[(u*nBin)+v]*S[(u*nBin)+v]) / P[(u*nBin)+v];
		}
	}
	if ((mode == 1) || (mode == 5)) {
		int lo = (int)(0.25*nBin);
		int mi = (int)(0.50*nBin);
		int hi = (int)(0.75*nBin);
		//double max = P[0][lo] + P[lo+1][hi] + P[hi+1][nBin-1];
		double max = P[lo] + P[((lo+1)*nBin)+mi] + P[((mi+1)*nBin)+hi] + P[((hi+1)*nBin)+255];
		for (int l = 0; l < (nBin-3); l++) {
			for (int m = l + 1; m < (nBin-2); m++) {
				for (int h = m + 1; h < (nBin-1); h++) {
					//double v = P[0][l]+P[l+1][h]+P[h+1][nBin-1];
					double v = P[l] + P[((l+1)*nBin)+m] + P[((m+1)*nBin)+h] + P[((h+1)*nBin)+255];
					if (v > max) {
						lo = l;
						mi = m;
						hi = h;
						max = v;
					} //new max
				}//for h -> hi
			} //for m -> mi
		} //for l -> low
		if (mode == 1)
			thresh = hi;
		else
			thresh = lo;
		*dark = lo;
		*mid = mi;
		*bright = hi;
	} else if ((mode == 2) || (mode == 4)) {
		int lo = (int)(0.33*nBin);
		int hi = (int)(0.67*nBin);
		double max = P[lo] + P[((lo+1)*nBin)+hi] + P[((hi+1)*nBin)+nBin-1];
		for (int l = 0; l < (nBin-2); l++) {
			for (int h = l + 1; h < (nBin-1); h++) {
				double v = P[l]+P[((l+1)*nBin)+h]+P[((h+1)*nBin)+nBin-1];
				if (v > max) {
					lo = l;
					hi = h;
					max = v;
				} //new max
			}//for h -> hi
		} //for l -> low
		if (mode == 1)
			thresh = hi;
		else
			thresh = lo;
		*dark = lo;
		*mid = thresh;
		*bright = hi;
	} else { //two levels:
		thresh = (int)(0.25*nBin); //nBin / 2;
		double max = P[thresh]+P[((thresh+1)*nBin)+nBin-1];
		//exhaustively search
		for (int i = 0; i < (nBin-1); i++) {
			double v = P[i]+P[((i+1)*nBin)+nBin-1];
			if (v > max) {
				thresh = i;
				max = v;
			}//new max
		}
		*dark = thresh;
		*mid = thresh;
		*bright = thresh;
	}
	free(P);
	free(S);
	return thresh;
}

float setThreshold(float*img, int nvox, int darkMediumBright123) {
	float mn, mx;
	if (nifti_robust_range(img, nvox, &mn, &mx, 0) != 0)
		return 1; //no variability
	#define kOtsuBins 256
	float scl = (kOtsuBins - 1) / (mx - mn);
	//create histogram
	int *hist = (int*) malloc(kOtsuBins * sizeof(int));
	for (int i = 0; i < kOtsuBins; i++)
		hist[i] = 0;
	for (int i = 0; i < nvox; i++) {
		if (isnanx(img[i]))
			continue;
		int idx = (int)round((img[i] - mn) * scl);
		idx = MIN(idx, kOtsuBins - 1);
		idx = MAX(idx, 0);
		hist[idx]++;
	}
	int dark, mid, bright;
	if ((darkMediumBright123 == 1) || (darkMediumBright123 == 3)) {
		//mode 5: Otsu multi-level segment to 4 intensities with 3 boundaries
		nii_otsu(hist, kOtsuBins, 5, &dark, &mid, &bright);
		free(hist);
		if (darkMediumBright123 == 1)
			return (dark / scl) + mn;
		if (darkMediumBright123 == 3)
			return (bright / scl) + mn;
	}
	//mode 3: classic Otsu binary threshold segment to 2 intensities with 1 boundary
	nii_otsu(hist, kOtsuBins, 3, &dark, &mid, &bright);
	free(hist);
	return (mid / scl) + mn;
}
