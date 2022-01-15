/*
 * Copyright (c) 2014 Cameron Hart
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */
#ifndef BITS_RADIXSORT_H
#define BITS_RADIXSORT_H

#include <stdint.h>
#include <stdlib.h>

#if _MSC_VER || __cplusplus
#define restrict __restrict
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern uint32_t radix8sort_u32(uint32_t* restrict keys_in_out, uint32_t* restrict keys_temp,
    uint32_t* restrict values_in_out, uint32_t* values_temp, uint32_t size);

extern uint32_t radix8sort_u64(uint64_t* restrict keys_in_out, uint64_t* restrict keys_temp,
    uint32_t* restrict values_in_out, uint32_t* values_temp, uint32_t size);

extern uint32_t radix8sort_f32(float* restrict keys_in_out, float* restrict keys_temp,
    uint32_t* restrict values_in_out, uint32_t* restrict values_temp, uint32_t size);

extern uint32_t radix11sort_u32(uint32_t* restrict keys_in, uint32_t* restrict keys_out,
    uint32_t* restrict values_in, uint32_t* restrict values_out, uint32_t size);

extern uint32_t radix11sort_u64(uint64_t* restrict keys_in_out, uint64_t* restrict keys_temp,
    uint32_t* restrict values_in_out, uint32_t* values_temp, uint32_t size);

extern uint32_t radix11sort_f32(float* restrict keys_in, float* restrict keys_out,
    uint32_t* restrict values_in, uint32_t* restrict values_out, uint32_t size);

#ifdef __cplusplus
}
#endif

#endif // BITS_RADIXSORT_H
