//
// Created by efti-nile on 07.03.2020.
//

#ifndef SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H
#define SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H

#include <stdint-gcc.h>

typedef struct {
    int16_t acc[3];  // raw values
    float temp;  // in Celsius degrees
} measurement_t;

typedef struct {
    float ax;
    float ay;
} angles_t;

typedef struct {
    float B[3];
    float C[3];
} temp_coefs_t;

typedef struct {
    float K[3*3];
    float b[3];
} axes_coefs_t;

int gauss_newton_calc(measurement_t *pmeas, int num_measurements, angles_t *angles,
        temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs);

#endif //SHMS_ORIENTATION_PROCESSING_CMODULE_GAUSS_NEWTON_ALGO_H
