//
// Created by efti-nile on 07.03.2020.
//

#include "gauss_newton_algo.h"

static inline void get_rotmat(double *rotmat, double cos_ax, double sin_ax, double cos_ay, double sin_ay) {
    rotmat[0] = cos_ay;  rotmat[1] = sin_ax * sin_ay; rotmat[2] = cos_ax * sin_ay;
    rotmat[3] = 0;       rotmat[4] = cos_ax;          rotmat[5] = -sin_ax;
    rotmat[6] = -sin_ay; rotmat[7] = cos_ay * sin_ax; rotmat[8] = cos_ax * cos_ay;
}

static inline void get_dgx(double *dgx, double cos_ax, double sin_ax, double cos_ay, double sin_ay) {
    dgx[0] = -sin_ax * sin_ay;
    dgx[1] = cos_ax * cos_ay;
}

static inline void get_dgy(double *dgy, double cos_ax) {
    dgy[0] = -cos_ax;
    dgy[1] = 0;
}
static inline void get_dgz(double *dgz, double cos_ax, double sin_ax, double cos_ay, double sin_ay) {
    dgz[0] = -cos_ay * sin_ax;
    dgz[1] = -cos_ax * sin_ay;
}

int gauss_newton_calc(double *accelerations, double temperatures[], unsigned  int num_measurements, angles_t *pangles,
                      temp_coefs_t *ptemp_coefs, axes_coefs_t *paxes_coefs) {
    // Preprocess data
    convert_temp(temperatures, num_measurements);
    compensate_temp(accelerations, temperatures, num_measurements, ptemp_coefs);
    g_direction_t g_dir = normalize(accelerations, num_measurements);
    compensate_axes(accelerations, num_measurements, paxes_coefs);

    // Gauss-Newton

    double angles[2] = {0};

    for (unsigned int iter_num = 0; iter_num < NUM_ITERATIONS; iter_num++) {
        // Precalculate trigonometry
        double cos_ax = cos(angles[0]), sin_ax = sin(angles[0]);
        double cos_ay = cos(angles[1]), sin_ay = sin(angles[1]);

        double rotmat[3*3];  // rotation matrix
        get_rotmat(rotmat, cos_ax, sin_ax, cos_ay, sin_ay);

        double g_mod[3], g_ref[3] = {0., 0., 1.};
        mul(g_ref, rotmat, false, g_mod, 1, 3, 3);

        double delta[num_measurements*3], delta_sum[3] = {0};
        reorder_axes(g_mod, g_dir);
        for (unsigned int i = 0; i < num_measurements; ++i) {
            sub(accelerations + 3 * i, g_mod, delta + 3 * i, 3, 1, 1);
        }
        sumrows(delta, delta_sum, num_measurements, 3);
        bring_axes_back(delta_sum, g_dir);

        double dx_dy_dz[2*3];
        get_dgx(dx_dy_dz, cos_ax, sin_ax, cos_ay, sin_ay);
        get_dgy(dx_dy_dz + 2, cos_ax);
        get_dgz(dx_dy_dz + 4, cos_ax, sin_ax, cos_ay, sin_ay);


        double dJ[2];
        mul(delta_sum, dx_dy_dz, false, dJ, 1, 3, 2);
        scale(dJ, 2., 1, 2);

        double Gx[2*2], Gy[2*2], Gz[2*2];
        mul(dx_dy_dz, dx_dy_dz, false, Gx, 2, 1, 2);  // dx-column * dx-row
        mul(dx_dy_dz + 2, dx_dy_dz + 2, false, Gy, 2, 1, 2);  // dy-column * dy-row
        mul(dx_dy_dz + 4, dx_dy_dz + 4, false, Gz, 2, 1, 2);  // dz-column * dz-row


        double G[2*2] = {0};
        add(G, Gx, G, 2, 2, 2);
        add(G, Gy, G, 2, 2, 2);
        add(G, Gz, G, 2, 2, 2);
        scale(G, 2 * (double) num_measurements, 2, 2);

        int retval = inv_2d((double (*)[2]) G);  // 1d G[4] casted to 2d G[2][2]
        if (retval != 0) {
            return ERROR_ZERO_DET;
        }

        double d_angle[2];
        mul(G, dJ, false, d_angle, 1, 2, 2);
        sub(angles, d_angle, angles, 2, 1, 1);

        double d_angles_sqnorm = d_angle[0] * d_angle[0] + d_angle[1] * d_angle[1];
        if (d_angles_sqnorm < EPSILON) {
            pangles->ax = angles[0] * (180.0 / M_PI);
            pangles->ay = angles[1] * (180.0 / M_PI);
            return 0;
        }
    }

    return ERROR_DOESNT_CONVERGE;
}

static const unsigned int axe_reorder_table[6][3] = {
        {2, 1, 0},
        {2, 0, 1},
        {0, 1, 2},
        {2, 1, 0},
        {2, 0, 1},
        {0, 1, 2}
};

static const bool axe_flips[6][3] = {
        {false, true, false},
        {false, false, false},
        {false, false, false},
        {false, false, true},
        {false, true, true},
        {false, true, true}
};

/*
    Reorders a vector of length 3 according the direction of g vector. This enables the sensor to work at any position
    without changes in the algorithm math.

    Example of reorder definition via order and flips arrays:
        vect3 = {0., 1., 2.}, order = {2, 0, 1}, flips = {false, false, true} result in {1., -2., 0.}
*/
static void reorder_axes(double vect3[3], g_direction_t g_dir) {
    double tmp[3];
    for (int i = 0; i < 3; ++i) {
        tmp[axe_reorder_table[g_dir][i]] = vect3[i];
        if (axe_flips[g_dir][i]) {
            tmp[axe_reorder_table[g_dir][i]] *= -1;
        }
    }
    memcpy(vect3, tmp, sizeof(tmp));
}

// This function is an inverse one for reorder_axes function
static void bring_axes_back(double vect3[3], g_direction_t g_dir) {
    double tmp[3];
    for (int i = 0; i < 3; ++i) {
        tmp[i] = vect3[axe_reorder_table[g_dir][i]];
        if (axe_flips[g_dir][i]) {
            tmp[i] *= -1;
        }
    }
    memcpy(vect3, tmp, sizeof(tmp));
}

// ADXL355 conversion function
static void convert_temp(double temperatures[], unsigned int num_measurements) {
    const double intercept_deg = 25, intercept_lsb = 1852, lsb_per_deg = -9.05;
    for (unsigned int i = 0; i < num_measurements; ++i) {
        temperatures[i] = intercept_deg + (temperatures[i] - intercept_lsb) / lsb_per_deg;
    }
}

static void compensate_temp(double *accelerations, double temperatures[],
        unsigned int num_measurements, temp_coefs_t *ptemp_coefs) {
    for (unsigned int i = 0; i < num_measurements; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            double tc = ptemp_coefs->B[j] + ptemp_coefs->C[j] * accelerations[i * 3 + j];
            accelerations[i * 3 + j] -= tc * temperatures[i];
        }
    }
}

static g_direction_t normalize(double *accelerations, unsigned int num_measurements) {
    // Normalize
    for (unsigned int i = 0; i < num_measurements; ++i) {
        double norm = l2_norm(accelerations + 3 * i, 3);
        scale(accelerations + 3 * i, 1 / norm, 1, 3);
    }

    // Select axe with max abs value
    double sum[3] = {0};
    sumrows(accelerations, sum, num_measurements, 3);
    unsigned int argmax = 0;
    double max = sum[argmax];
    double max_abs = fabs(max);
    for (unsigned int i = 1; i < 3; ++i) {
        double tmp = fabs(sum[i]);
        if (tmp > max_abs) {
            argmax = i;
            max_abs = tmp;
            max = sum[i];
        }
    }

    // Get sensor position
    if (max > 0) {
        return (g_direction_t) argmax;  // positive directions
    } else {
        return (g_direction_t) (argmax + X_NEGATIVE);  // negative directions
    }
}

// EmbeddedLapack's norm function doesn't suit so here's mine
static double l2_norm(const double *vect, unsigned int len) {
    double sum = .0;
    for (unsigned int i = 0; i < len; ++i) {
        sum += vect[i] * vect[i];
    }
    return sqrt(sum);
}

// Inverts 2x2 matrix
static int inv_2d(double A[2][2]) {
    double T[2][2];
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (fabs(det) > FLOAT_ZERO) {
        T[0][0] = A[1][1]; T[0][1] = -A[0][1];
        T[1][0] = -A[1][0]; T[1][1] = A[0][0];
        scale((double *) T, 1 / det, 2, 2);  // 1d T[2][2] casted to 2d T[4]
        memcpy(A, T, sizeof(T));
        return 0;
    } else {
        return -1;  // The det is zero
    }
}

static void compensate_axes(double *accelerations, unsigned int num_measurements, axes_coefs_t *paxes_coefs) {
    for (unsigned int i = 0; i < num_measurements; ++i) {
        double K_mul_acc[3];  // a buffer
        mul(paxes_coefs->K, accelerations + 3 * i, false, K_mul_acc, 3, 3, 1);

        // write compensated values
        add(K_mul_acc, paxes_coefs->b, accelerations + 3 * i, 3, 1, 1);
    }
}
