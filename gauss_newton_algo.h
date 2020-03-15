//
// Created by efti-nile on 07.03.2020.
//

#ifndef SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H
#define SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H

#include <stdint-gcc.h>
#include <glob.h>
#include <math.h>

#include "LinearAlgebra/declareFunctions.h"

#define NUM_ITERATIONS 20
#define GAUSS_NEWTON_EPSILON 1e-6

typedef struct {
    double ax;
    double ay;
} angles_t;

typedef struct {
    double B[3];
    double C[3];
} temp_coefs_t;

typedef struct {
    double K[3*3];
    double b[3];
} axes_coefs_t;

// Direction of g-vector
typedef enum {
    X_POSITIVE=0, Y_POSITIVE,   Z_POSITIVE,
    X_NEGATIVE,   Y_NEGATIVE,   Z_NEGATIVE
} g_direction_t;

int gauss_newton_calc(double *, double temperatures[], unsigned  int num_measurements, angles_t *pangles,
        temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs, g_direction_t dir);

// Local functions
static void reorder_axes(double *vect3, unsigned int order[3], bool flips[3]);
static void convert_temp(double temperatures[], unsigned int num_measurements);
static void compensate_temp(double *accelerations, double temperatures[],
                            unsigned int num_measurements, temp_coefs_t *ptemp_coefs);
static g_direction_t normalize(double *accelerations, unsigned int num_measurements);
static double l2_norm(double vect[], unsigned int len);
static void compensate_axes(double *accelerations, unsigned int num_measurements, axes_coefs_t *paxes_coefs);

#endif //SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H
