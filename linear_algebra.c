//
// Created by efti-nile on 17.03.2020.
//

// A little extract of EmbeddedLAPACK library

#include "linear_algebra.h"

/*
 * Add C = A+B with A size row_a x column_a and B have the size row_a x column_b and C have the size row_a x column_a
 */

void add(double* A, double* B, double* C, int row_a, int column_a, int column_b) {

    /*
     * This uses the formula C=A+B
     */

    // Check which one should be the column counter
    int columnMatrix = column_a; // Initial select
    if(column_b > column_a) // But if....
        columnMatrix = column_b;

    // Add all values
    for (int i = 0; i < row_a; i++) {
        for (int j = 0; j < columnMatrix; j++) {

            if (column_b == 1 && column_a > 1)
                *(C++) = *(A++) + *(B + i); // Matrix + Vector
            else if (column_b > 1 && column_a == 1)
                *(C++) = *(A + i) + *(B++); // Vector + Matrix
            else
                *(C++) = *(A++) + *(B++); // Matrix + Matrix or Vector + Vector

        }
    }
}

/*
 * Multiply A with size row_a x column_a with matrix B with size row_a x column_b and get matrix C with row_a x column_b
 */

void mul(double* A, double* B, bool elementWise, double* C, int row_a, int column_a, int column_b) {

    /*
     * C = A*B if elementWise = false
     * C = A.*B if elementWise = true
     */

    // Data matrix
    double* data_a = A;
    double* data_b = B;
    //double* C = C; // No need

    if (elementWise == true) {

        // Dot multiply all values
        if (column_b > 1) { // If matrix b is a matrix
            for (int i = 0; i < row_a; i++) {
                for (int j = 0; j < column_a; j++) {
                    // Do element wise mutiplication. In MATLAB it is A.*A
                    *(C++) = *(data_a++) * *(data_b++);
                }
            }
        } else {
            // If matrix b is a vector
            for (int i = 0; i < row_a; i++) {
                for (int j = 0; j < column_a; j++) {
                    // Do element wise mutiplication. In MATLAB it is A.*b
                    *(C++) = *(data_a++) * *(data_b + i);
                }
            }
        }

    } else {
        // Do regular mutiplication. In MATLAB it is A*A
        // Let's take our a matrix
        for (int i = 0; i < row_a; i++) {

            // Then we go through every column of b
            for (int j = 0; j < column_b; j++) {
                data_a = &A[i * column_a];
                data_b = &B[j];

                *C = 0; // Reset
                // And we multiply rows from a with columns of b
                for (int k = 0; k < column_a; k++) {
                    *C += *data_a * *data_b;
                    data_a++;
                    data_b += column_b;
                }
                C++;
            }
        }
    }
}

/*
 * Sub C = A-B with A size row_a x column_a and B have the size row_a x column_b and C have the size row_a x column_a
 */

void sub(double* A, double* B, double* C, int row_a, int column_a, int column_b) {

    /*
     * This uses the formula C=A-B
     */

    // Check which one should be the column counter
    int columnMatrix = column_a; // Initial select
    if(column_b > column_a) // But if....
        columnMatrix = column_b;

    // Add all values
    for (int i = 0; i < row_a; i++) {
        for (int j = 0; j < columnMatrix; j++) {

            if (column_b == 1 && column_a > 1)
                *(C++) = *(A++) - *(B + i); // Matrix - Vector
            else if (column_b > 1 && column_a == 1)
                *(C++) = *(A + i) - *(B++); // Vector - Matrix
            else
                *(C++) = *(A++) - *(B++); // Matrix - Matrix or Vector - Vector

        }
    }
}

/*
 * Scale a matrix A with a scalar value. Size of matrix A is row x column
 */
void scale(double* A, double scalar, int row, int column) {

    for (int i = 0; i < row*column; i++)
        *(A + i) = *(A + i) * scalar;

}

/*
 * Sum all rows from A size row x column, into a flat matrix B size 1 x column
 */

void sumrows(double* A, double* B, int row, int column) {

    for(int j = 0; j < column; j++){
        for(int i = 0; i < row; i++){
            *(B + j) += *((A + i*column) + j);
        }
    }

}
