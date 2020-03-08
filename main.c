#include <stdio.h>
#include <time.h>
//#include "LinearAlgebra/declareFunctions.h"

#include "orientation_processing.h"

int main() {
    int ids[3] = {5374047, 5374049, 5374077};
    orientation_t or;
    measurement_t m;
    m.acc[0] = 1.;
    m.acc[1] = 2.;
    m.acc[2] = 3.;
    m.temp = 25.7;
    for (int i = 0; i < 9; i++) {
        process(&m, ids[i % 3], &or);
    }
//    clock_t start, end;
//    float cpu_time_used;
//    start = clock();
//
//    // A matrix with size 6 x 4
//    double A[6*4] = {0.674878,   0.151285,   0.875139,   0.150518,
//                     0.828102,   0.150747,   0.934674,   0.474325,
//                     0.476510,   0.914686,   0.740681,   0.060455,
//                     0.792594,   0.471488,   0.529343,   0.743405,
//                     0.084739,   0.475160,   0.419307,   0.628999,
//                     0.674878,   0.151285,   0.875139,   0.150518};
//
//    double Q[6*6];
//    double R[6*4];
//
//    // Solve
//    qr(A, Q, R, 6, 4);
//
//    // Print
//    print(A, 6,4);
//    print(Q, 6,6);
//    print(R, 6,4);
//
//    end = clock();
//    cpu_time_used = ((float) (end - start)) / CLOCKS_PER_SEC;
//    printf("\nTotal speed  was %f,", cpu_time_used);
    return 0;
}