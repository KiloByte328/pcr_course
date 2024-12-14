#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int rank, commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srand(rank);
    size_t n = 2000;
    size_t m = n + 1;
    size_t my_part = n / commsize;
    double* Matr = malloc(sizeof(double) * (m * my_part));
    double* Others = malloc(sizeof(double) * (m * commsize));
    for (size_t y = 0; y < my_part; y++)
    {
        for (size_t x = 0; x < m; x++)
        {
            Matr[y*m + x] = rand() % (RAND_MAX - 100);
        }
    }
    MPI_Datatype row;
    MPI_Type_contiguous(m, MPI_DOUBLE, row);
    MPI_Type_commit(&row);
    double now = MPI_Wtime();
    double Div;
    for (int grid_walker = 0; grid_walker < n; grid_walker++)
    {
        if (grid_walker <= rank == commsize - 1 ? n : ((rank + 1) * my_part) - 1)
        {
            if (Matr[grid_walker * m + grid_walker] != 1) {
            for (size_t x = 0; x < m; x++)
            {
                x == grid_walker ? x++ : x;
                Matr[grid_walker * m + x] = Matr[grid_walker * m + x] / Matr[grid_walker * m + grid_walker];
            }
            Matr[grid_walker * m + grid_walker] >= 1 ? (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] / Matr[grid_walker * m + grid_walker]) :
            (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] * (1 / Matr[grid_walker * m + grid_walker]));
        }
        MPI_Bcast(&Matr[grid_walker * m], 1, row, (grid_walker / my_part), MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast(&Others, 1, row, (grid_walker / my_part), MPI_COMM_WORLD);
        }
    }
/*     for (size_t grid_walker = 0; grid_walker < my_part; grid_walker++)
    {
         */
/*         if (rank == 0)
        {
            if (Matr[grid_walker * m + grid_walker] != 1) {
            for (size_t x = 0; x < m; x++)
            {
                x == grid_walker ? x++ : x;
                Matr[grid_walker * m + x] = Matr[grid_walker * m + x] / Matr[grid_walker * m + grid_walker];
            }
            Matr[grid_walker * m + grid_walker] >= 1 ? (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] / Matr[grid_walker * m + grid_walker]) :
            (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] * (1 / Matr[grid_walker * m + grid_walker]));
        }
        MPI_Bcast(&Matr[grid_walker * m], 1, row, 0, MPI_COMM_WORLD);
        for (size_t y = 0; y < my_part; y++)
        {
            y + 1 != n ? (y == grid_walker ? y++ : y) : y;
            Div = Matr[y * m + grid_walker] / Matr[grid_walker * m + grid_walker];
            for(size_t x = grid_walker; x < m; x++)
            {
                if (x == grid_walker)
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + grid_walker]);
                else
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + x]);
            }
        }
        }
        else
        {
        size_t true_grid_walker = grid_walker * rank;
        MPI_Bcast(&Others, 1, row, 0, MPI_COMM_WORLD);
        for (size_t y = 0; y < my_part; y++)
        {
            y + 1 != n ? (y == grid_walker ? y++ : y) : y;
            Div = Others[y * m + grid_walker] / Matr[grid_walker * m + grid_walker];
            for(size_t x = grid_walker; x < m; x++)
            {
                if (x == grid_walker)
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + grid_walker]);
                else
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + x]);
            }
        }
        } */
        /* if (Matr[grid_walker * m + grid_walker] != 1) {
            for (size_t x = 0; x < m; x++)
            {
                x == grid_walker ? x++ : x;
                Matr[grid_walker * m + x] = Matr[grid_walker * m + x] / Matr[grid_walker * m + grid_walker];
            }
            Matr[grid_walker * m + grid_walker] >= 1 ? (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] / Matr[grid_walker * m + grid_walker]) :
            (Matr[grid_walker * m + grid_walker] = Matr[grid_walker * m + grid_walker] * (1 / Matr[grid_walker * m + grid_walker]));
        }
        MPI_Bcast()
        for (size_t y = 0; y < my_part; y++)
        {
            y + 1 != n ? (y == grid_walker ? y++ : y) : y;
            Div = Matr[y * m + grid_walker] / Matr[grid_walker * m + grid_walker];
            for(size_t x = grid_walker; x < m; x++)
            {
                if (x == grid_walker)
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + grid_walker]);
                else
                    Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + x]);
            }
        } */
/*     }
    for (size_t grid_back = n - 1; grid_back >= 1; grid_back--)
    {
        for (size_t y = grid_back - 1; y > 0; y--) 
        {
            Div = Matr[(grid_back - 1) * m + grid_back] / Matr[grid_back * m + grid_back];
            for (size_t x = m - 1; x >= grid_back; x--)
            {
                Matr[y * m + x] = Matr[y * m + x] - Matr[grid_back * m + x] * Div;
            }
        }
    } */
    double end = MPI_Wtime() - now;
    printf("Time of work: %f\n", end);
    MPI_Finalize();
    return 0;
}