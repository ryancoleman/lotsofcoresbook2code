/* FILE LICENSE TAG: SAMPLE */

void mean_and_stdev(double time[], double &mean, double &stdev, int n)
{
    mean = 0.0;
    stdev = 0.0;
    int i;

    for (i = 1; i < n; ++i) {
        //ignoring first iteration
        mean += time[i];
    }

    mean = mean / (n - 1);

    for (i = 1; i < n; ++i) {
        stdev += pow((time[i] - mean), 2);
    }

    stdev = sqrt(stdev / (n - 2));
}

void eyeInit(double *mat, int size)
{
    int i, j;
    memset(mat, 0, size * size * sizeof(double));
    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j) {
            if (i == j) {
                mat[i * size + j] = 1.0;
            }
        }
}

void read_mat(double *mat, int size, char filename[])
{
    int i, j;
    FILE *fin;
    fin = fopen(filename, "r");
    if (fin == NULL) {
        perror("Error");
        exit(0);
    }
    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j) {
            fscanf(fin, "%lg\n", &mat[i * size + j]);
        }

    fclose(fin);
}

void print_result(double *mat, int size, int print_size)
{
    for (int i = 0; i < print_size; ++i) {
        for (int j = 0; j < print_size; ++j) {
            printf("%.4f  ", mat[i * size + j]);
        }
        printf("\n");
    }

}

bool verify_results(double *act, double *ref, int size)
{
    double diff;

    bool res = true;
    for (int i = 0; i < size; ++i) {
        diff = ref[i] - act[i];
        if (fabs(ref[i]) > 1e-5) {
            diff /= ref[i];
        }
        diff = fabs(diff);
        if (diff > 1.0e-5) {
            printf("\nError detected at i = %d: ref %g actual %g\n",
                   i, ref[i], act[i]);
            res = false;
            break;
        }
    }
    return res;
}

void fill_upper_triangle(double *mat, int size)
{
    unsigned int i, j;
    for (i = 0; i < size; ++i)
        for (j = i + 1; j < size; ++j) {
            mat[i * size + j] = mat[j * size + i];
        }
}

void split_into_blocks(double *mat, double *mat_split[], int num_tiles, int tile_size, int size, bool layRow)
{
    int itile, i, j, offset_tile;

    int tot_tiles = num_tiles * num_tiles;

    #pragma omp parallel for private(i, j, offset_tile) schedule (auto)
    for (itile = 0; itile < tot_tiles; ++itile) {
        if (layRow) {
            offset_tile = int(itile / num_tiles) * num_tiles * tile_size * tile_size
                          + int(itile % num_tiles) * tile_size;
        } else {
            offset_tile = int(itile % num_tiles) * num_tiles * tile_size * tile_size
                          + int(itile / num_tiles) * tile_size;
        }

        for (i = 0; i < tile_size; ++i)
#pragma simd
            for (j = 0; j < tile_size; ++j) {
                mat_split[itile][i * tile_size + j] = mat[offset_tile + i * size + j];
            }
    }
}

void assemble(double *mat_split[], double *mat, int num_tiles, int tile_size, int size, bool layRow)
{
    int i_tile, j_tile, tile, i, j, i_local, j_local;
    #pragma omp parallel for private(j, i_local, j_local, i_tile, j_tile, tile) schedule (auto)
    for (i = 0; i < size; ++i) {
        i_local = int(i % tile_size);
        i_tile = int(i / tile_size);
#pragma simd private(j_tile, tile, j_local)
        for (j = 0; j < size; ++j) {
            j_tile  = int(j / tile_size);
            if (layRow) {
                tile = i_tile * num_tiles + j_tile;
            } else {
                tile = j_tile * num_tiles + i_tile;
            }
            j_local = int(j % tile_size);
            mat[i * size + j] = mat_split[tile][i_local * tile_size + j_local];
        }
    }
}

void copy_mat(double *A, double *B, int size)
{
    int i, j;
    for (i = 0; i < size; ++i)
        for (j = 0; j < size; ++j) {
            B[i * size + j] = A[i * size + j];
        }
}
