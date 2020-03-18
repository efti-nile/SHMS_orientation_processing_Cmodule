//
// Created by efti-nile on 17.03.2020.
//

// A little extract of EmbeddedLAPACK library

#ifndef SHMS_ORIENTATION_PROCESSING_CMODULE_LINEAR_ALGEBRA_H
#define SHMS_ORIENTATION_PROCESSING_CMODULE_LINEAR_ALGEBRA_H

#include <stdbool.h>

void scale(double* A, double scalar, int row, int column);
void mul(double* A, double* B, bool elementWise, double* C, int row_a, int column_a, int column_b);
void sub(double* A, double* B, double* C, int row_a, int column_a, int column_b);
void sumrows(double* A, double* B, int row, int column);
void add(double* A, double* B, double* C, int row_a, int column_a, int column_b);

#endif //SHMS_ORIENTATION_PROCESSING_CMODULE_LINEAR_ALGEBRA_H
