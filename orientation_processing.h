//
// Created by efti-nile on 07.03.2020.
//

#ifndef SHMS_ORIENTATION_PROCESSING_CMODULE_ORIENTATION_PROCESSING_H
#define SHMS_ORIENTATION_PROCESSING_CMODULE_ORIENTATION_PROCESSING_H

#include <stdio.h>
#include <string.h>

#include "uthash.h"
#include "gauss_newton_algo.h"

#define CALIB_COEFS_PATH "./calib_coefs"
#define ACCEL_NAME "adxl355"
#define GAUSS_NEWTON_WINDOW_SIZE 25

#define PROCESS_RESULT_ERROR (-1)
#define PROCESS_RESULT_OK 0
#define PROCESS_RESULT_NEW_ORIENTATION 1

typedef struct {
    angles_t angles;
    long int timestamp;
} orientation_t;

typedef struct {
    int id;
    measurement_t data_buf[GAUSS_NEWTON_WINDOW_SIZE];
    long int timestamps[GAUSS_NEWTON_WINDOW_SIZE];
    int data_size;
    temp_coefs_t temp_coefs;
    axes_coefs_t axes_coefs;

    UT_hash_handle hh;
} sensor_t;

int process(measurement_t *pmeas, int id, orientation_t *porient);
sensor_t *find_sensor(int id);
int load_calib_coefs(int id, temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs);
int load_temp_coefs(const char *file_name, temp_coefs_t *ptemp_coefs);
int load_K_matrix(const char *file_name, float *matrix);
int load_b_vector(const char *file_name, float *vector);
void add_sensor(sensor_t *psensor);

#endif //SHMS_ORIENTATION_PROCESSING_CMODULE_ORIENTATION_PROCESSING_H
