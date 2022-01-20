#ifndef EX_BWLABEL_H
#define EX_BWLABEL_H

#ifdef  __cplusplus
//extern "C" {
#endif

#include <stdbool.h>

int bwlabel(float *img, int conn, size_t dim[3], bool onlyLargest, bool fillBubbles);
//if 
// onlyLargest = 0: return image where voxel intensity is cluster number
// onlyLargest = 1: return image where voxels of largest cluster is 1, all other voxels 0
// onlyLargest = 2: return image where voxels of and inside largest cluster is 1, all other voxels 0


#ifdef  __cplusplus
//}
#endif

#endif // EX_BWLABEL_H
