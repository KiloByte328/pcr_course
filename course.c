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
    double* Others = malloc(sizeof(double) * m);
    for (size_t y = 0; y < my_part; y++)
    {
        for (size_t x = 0; x < m; x++)
        {
            Matr[y * m + x] = rand() % 1000;
            while (Matr[y * m + x] == 0.0) Matr[y * m + x] = rand() % 1000;
        }
    }
    MPI_Datatype row;
    MPI_Type_contiguous(m, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);
    double now = MPI_Wtime();
    double Div;
    for (size_t grid_walker = 0; grid_walker < n; grid_walker++)
    {
        if (grid_walker / my_part == rank)
        {
            Div = Matr[(grid_walker % my_part) * m + grid_walker];
            if (Div == 0)
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            if (Div != 1)
            {
                for (size_t x = grid_walker; x < m; x++)
                {
                    Matr[(grid_walker % my_part) * m + x] = Matr[grid_walker % my_part * m + x] >= 1 ?
                    Matr[grid_walker % my_part * m + x] / Div 
                    : Matr[grid_walker % my_part * m + x] * (1 / Div);
                }
            }
            MPI_Bcast(&Matr[grid_walker % my_part], 1, row, grid_walker / my_part, MPI_COMM_WORLD);
            for (size_t y = grid_walker % my_part + 1; y < my_part; y++)
            {
                if (rank != 0)
                    Div = Matr[y * m + grid_walker] >= 1 ?
                    Matr[y * m + grid_walker] / Others[grid_walker] :
                    Matr[y * m + grid_walker] * (1 / Others[grid_walker]);
                else 
                    Div = Matr[y * m + grid_walker] >= 1 ?
                    Matr[y * m + grid_walker] / Matr[grid_walker % my_part * m + grid_walker] :
                    Matr[y * m + grid_walker] * (1 / Matr[grid_walker % my_part * m + grid_walker]);
                for (size_t x = grid_walker; x < m; x++)
                {
                    Matr[y * m + x] = Matr[y * m + x] - Div * Matr[grid_walker % my_part * m + x];
                }
            }
        }
        else
        {
            MPI_Bcast(&Others[0], 1, row, grid_walker / my_part, MPI_COMM_WORLD);
            if (grid_walker / my_part > rank)
            for (size_t y = grid_walker % my_part + 1; y < my_part; y++)
            {
                Div = Matr[y * m + grid_walker] >= 1 ?
                Matr[y * m + grid_walker] / Others[grid_walker] :
                Matr[y * m + grid_walker] * (1 / Others[grid_walker]);
                for (size_t x = grid_walker; x < m; x++)
                {
                    Matr[y * m + x] = Matr[y * m + x] - Div * Others[x];
                }
            }
        }
    }
    for (size_t grid_back = n - 1; grid_back >= 1; grid_back--)
    {
        if (grid_back / my_part == rank)
        {
            Div = Matr[grid_back % my_part * m + grid_back];
            for (size_t x = m; x >= grid_back; x--)
            {
                Matr[(grid_back - 1) % my_part * m + x] -= Matr[grid_back % my_part * m + x] * Div;
            }
            MPI_Bcast(&Matr[grid_back - 1], 1, row, grid_back / my_part, MPI_COMM_WORLD);
            for (size_t y = (grid_back - 2) % my_part; y > 0; y--)
            {
                Div = Matr[y * m + grid_back];
                for (size_t x = m - 1; x >= grid_back; x--)
                {
                    Matr[y * m + x] -= Matr[y * m + x] * Div;
                }
            }
        }
        else
        {
            MPI_Bcast(&Others[0], 1, row, grid_back / my_part, MPI_COMM_WORLD);
            for (size_t y = (grid_back - 1) % my_part; y > 0; y--)
            {
                Div = Matr[y  * m + grid_back];
                for (size_t x = m - 1; x >= grid_back; x--)
                {
                    Matr[y * m + x] -= Matr[y * m + x] * Div;
                }
            }
        }
    }
    if (rank == commsize - 1)
    {
        double end = MPI_Wtime() - now;
        printf("Time of work: %f\n", end);
    }
    MPI_Finalize();
    return 0;
}