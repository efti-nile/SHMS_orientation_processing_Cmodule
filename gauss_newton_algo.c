//
// Created by efti-nile on 07.03.2020.
//

#include "gauss_newton_algo.h"

int gauss_newton_calc(measurement_t *measurements, int num_measurements, angles_t *angles,
                      temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs) {
    angles->ax = -1.;
    angles->ay = 1.;
    return 0;
}