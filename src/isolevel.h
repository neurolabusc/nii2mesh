/*
Use Otsu's method to detect a reasonable threshold
*/

#ifndef ISOLEVEL_H
#define ISOLEVEL_H

#include <stdio.h>

float setThreshold(float*img, int nvox, int darkMediumBright123);

#endif /* ISOLEVEL_H */