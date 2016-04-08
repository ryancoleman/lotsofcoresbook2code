/* FILE LICENSE TAG: SAMPLE */

#include <stdlib.h> // for malloc, rand_r and RAND_MAX
#include <stddef.h> // for size_t
#include <mkl.h>

/**
 * Generate a double-precision real-valued hermitian matrix.
 *
 * @param   side_size Number of rows/columns
 * @return  A raw pointer to the contents of the matrix.
 *          The allocated memory is side_size*side_size*sizeof(double) bytes.
 *          Memory is allocated through malloc. Release it with free(), not delete[]
 *
 * The approach here is to generate a random symmetric matrix M and multiply it by its transpose.
 * The resulting matrix is Hermitian if the M matrix is nonsingular.
 * The probability of a random matrix being non-singular approaches 1 as the size of
 * the matrix grows.
 */
double *dpo_generate(size_t side_size)
{
    unsigned int seed = side_size;
#ifdef _WIN32
    srand(seed);
#endif
    // M is a (very) pseudo-random symmetric matrix
    double *M_matrix = new double[side_size * side_size];
    for (size_t row = 0; row < side_size; ++row) {
        for (size_t col = row; col < side_size; ++col) {
            M_matrix[col * side_size + row] = M_matrix[row * side_size + col]
                                              = (double)
#ifdef _WIN32
                                                rand()
#else
                                                rand_r(&seed)
#endif
                                                / RAND_MAX;
        }
    }
    double *ret_matrix = (double *) malloc(side_size * side_size * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                side_size, side_size, side_size,
                1.0,
                M_matrix, side_size,
                M_matrix, side_size,
                0.0,
                ret_matrix, side_size);

    //adjust diagonals (diag = sum (row entries) + 1.0)
    for (size_t row = 0; row < side_size; ++row) {
        double diag = 1.0; //start from 1.0
        for (size_t col = 0; col < side_size; ++col) {
            diag += ret_matrix[row * side_size + col];
        }
        //set the diag entry
        ret_matrix[row * side_size + row] = diag;
    }

    delete [] M_matrix;
    return ret_matrix;
}
