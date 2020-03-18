//
// Created by efti-nile on 07.03.2020.
//

#include "orientation_processing.h"

sensor_t *sensors = NULL;

int process(measurement_t *pmeas, int id, long int timestamp, orientation_t *porient) {
    sensor_t *psensor = find_sensor(id);  // try to find...

    if (psensor == NULL) {  // if not found create a new sensor instance
        psensor = (sensor_t *) malloc(sizeof(sensor_t));
        psensor->data_size = 0;
        psensor->id = id;
        if (load_calib_coefs(id, &psensor->temp_coefs,  &psensor->axes_coefs)) {
            printf("Can't load calibration coefficients!\n");
            return PROCESS_RESULT_ERROR;
        }
        add_sensor(psensor);
    }

    // Add new measurement
    for (int i = 0; i < 3; ++i) {
        psensor->accelerations[3 * psensor->data_size + i] = (double) pmeas->acc[i];
    }
    psensor->temperatures[psensor->data_size] = (double) pmeas->temp;
    psensor->timestamps[psensor->data_size] = timestamp;
    psensor->data_size++;

    if (psensor->data_size >= GAUSS_NEWTON_WINDOW_SIZE) {
        psensor->data_size = 0;
        int retval = gauss_newton_calc(psensor->accelerations, psensor->temperatures, GAUSS_NEWTON_WINDOW_SIZE,
                &porient->angles, &psensor->temp_coefs, &psensor->axes_coefs);
        if (retval == 0) {
            porient->timestamp = psensor->timestamps[GAUSS_NEWTON_WINDOW_SIZE-1];
            return PROCESS_RESULT_NEW_ORIENTATION;
        } else {
            if (retval == ERROR_DOESNT_CONVERGE) {
                printf("Gauss-Newton doesn't converge!\n");
            } else if (retval == ERROR_ZERO_DET) {
                printf("Gauss-Newton cant't inverse matrix G!\n");
            } else {
                printf("Gauss-Newton's circuits dead!\n");
            }

            return PROCESS_RESULT_ERROR;
        }
    }

    return PROCESS_RESULT_OK;
}

int load_calib_coefs(int id, temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs) {
    char file_pref[256] = "";
    sprintf(file_pref, CALIB_COEFS_PATH "/%d " ACCEL_NAME " ", id);

    // Read temperature coefficients
    char file_name[256] = "";
    strcpy(file_name, file_pref);
    strcat(file_name, "temp coefs.csv");
    if (load_temp_coefs(file_name, ptemp_coefs)) {
        return -3;
    };

    // Read K-matrix of axes calibration
    strcpy(file_name, file_pref);
    strcat(file_name, "K-matrix.csv");
    if (load_K_matrix(file_name, paxes_coefs->K)) {
        return -2;
    };

    // Read b-vector of axes calibration
    strcpy(file_name, file_pref);
    strcat(file_name, "b-vector.csv");
    if (load_b_vector(file_name, paxes_coefs->b)) {
        return -1;
    }

    return 0;
}

/*
 An example of temperature coefficients file:

,B,C
ADXL_ACC_X,0.2969321110239852,3.769493334012686e-05
ADXL_ACC_Y,-1.5055177940329225,4.2737141097391814e-05
ADXL_ACC_Z,1.166537726128811,3.0121018928494783e-06

*/
int load_temp_coefs(const char *file_name, temp_coefs_t *ptemp_coefs) {
    FILE *pfile = fopen(file_name, "r");
    if (pfile != NULL) {
        for (int line_num = 0; line_num < 4; ++line_num) {
            char line[256];
            if (fgets(line, sizeof(line), pfile) != NULL) {
                if (line_num > 0) {  // skip the first line
                    char *rest = line, *token;
                    for (int col_num = 0; col_num < 3; ++col_num) {
                        if ((token = strtok_r(rest, ",", &rest)) != NULL) {  // get CSV-token
                            if (col_num == 1) {
                                ptemp_coefs->B[line_num - 1] = atof(token);
                            } else if (col_num == 2) {
                                ptemp_coefs->C[line_num - 1] = atof(token);
                            }
                        } else {
                            printf("%s:%d: not enough values!", file_name, line_num);
                            return -3;
                        }
                    }
                }
            } else {
                printf("%s: not enough lines!", file_name);
                return -2;
            }
        }
        fclose(pfile);
        return 0;
    } else {
        printf("%s: can't open!", file_name);
        return -1;
    }
}

/*
 An example of K-matrix file:

1.000017558451608002e+00,4.368807634697606697e-04,1.099593122021705376e-04
4.369595974757343082e-04,9.999815531531026647e-01,-1.149758119888395292e-03
1.099419090058375035e-04,-1.149933568443346115e-03,9.999328195863038671e-01

*/
int load_K_matrix(const char *file_name, double *matrix) {
    FILE *pfile = fopen(file_name, "r");
    if (pfile != NULL) {
        for (int line_num = 0; line_num < 3; ++line_num) {
            char line[256];
            if (fgets(line, sizeof(line), pfile) != NULL) {
                char *rest = line, *token;
                for (int col_num = 0; col_num < 3; ++col_num) {
                    if ((token = strtok_r(rest, ",", &rest)) != NULL) {  // get CSV-token
                        matrix[line_num * 3 + col_num] = atof(token);
                    } else {
                        printf("%s:%d: not enough values!", file_name, line_num);
                        return -3;
                    }
                }
            } else {
                printf("%s: not enough lines!", file_name);
                return -2;
            }
        }
        fclose(pfile);
        return 0;
    } else {
        printf("%s: can't open!", file_name);
        return -1;
    }
}

/*
 An example of b-vector file:

-1.501509374113828288e-05
-3.548651774164593431e-05
1.436852599883027548e-05

*/
int load_b_vector(const char *file_name, double *vector) {
    FILE *pfile = fopen(file_name, "r");
    if (pfile != NULL) {
        for (int line_num = 0; line_num < 3; ++line_num) {
            char line[256];
            if (fgets(line, sizeof(line), pfile) != NULL) {
                vector[line_num] = atof(line);
            } else {
                printf("%s: not enough lines!", file_name);
                return -2;
            }
        }
        fclose(pfile);
        return 0;
    } else {
        printf("%s: can't open!", file_name);
        return -1;
    }
}

sensor_t *find_sensor(int id) {
    sensor_t *psensor;
    HASH_FIND_INT(sensors, &id, psensor);
    return psensor;
}

void add_sensor(sensor_t *psensor) {
    HASH_ADD_INT(sensors, id, psensor);
}
