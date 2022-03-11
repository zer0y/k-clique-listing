#ifndef _SET_INTERSECTION_H
#define _SET_INTERSECTION_H
#include <x86intrin.h>
#include <unistd.h>
#include <cstdint>
#include <sys/time.h>
#include "immintrin.h"
int intersect(const int *set_a, int size_a, const int *set_b, int size_b, int *set_c);
int intersect_count(const int *set_a, int size_a, const int *set_b, int size_b);
int intersect_simd4x(const int *set_a, int size_a, const int *set_b, int size_b, int *set_c);
int intersect_simd4x_count(const int* set_a, int size_a, const int* set_b, int size_b);
#endif

