#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int get_rnk(int commsize, int n, int who) {
  return who / (n / commsize);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank, commsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  srand(rank);
  int n = 1000;
  int m = n + 1;
  int my_begin = rank == 0 ? 0 : rank * (n / commsize);
  int my_end = rank == commsize - 1 ? n : (rank + 1) * (n / commsize);
  int my_part = my_end - my_begin;
  double* Matr = malloc(sizeof(double) * (m * my_part));
  double* Others = malloc(sizeof(double) * m);
  for (size_t y = 0; y < my_part; y++) {
    for (size_t x = 0; x < m; x++) {
      Matr[y * m + x] = rand() % 100;
      while (Matr[y * m + x] == 0.0)
        Matr[y * m + x] = rand() % 100;
    }
  }
  MPI_Datatype row;
  MPI_Type_contiguous(m, MPI_DOUBLE, &row);
  MPI_Type_commit(&row);
  double now = MPI_Wtime();
  double Div;
  for (int grid_walker = 0; grid_walker < n; grid_walker++) {
    if (grid_walker >= my_begin && grid_walker < my_end) {
      Div = Matr[(grid_walker % my_part) * m + grid_walker];

      for (size_t x = grid_walker; x < m; x++) {
        Matr[(grid_walker % my_part) * m + x] =
            Matr[grid_walker % my_part * m + x] >= 1
                ? Matr[grid_walker % my_part * m + x] / Div
                : Matr[grid_walker % my_part * m + x] * (1 / Div);
      }
      MPI_Bcast(&Matr[grid_walker % my_part], 1, row, rank, MPI_COMM_WORLD);
      for (size_t y = grid_walker % my_part + 1; y < my_part; y++) {
        Div = Matr[y * m + grid_walker] >= 1
                  ? Matr[y * m + grid_walker] /
                        Matr[grid_walker % my_part * m + grid_walker]
                  : Matr[y * m + grid_walker] *
                        (1 / Matr[grid_walker % my_part * m + grid_walker]);
        for (size_t x = grid_walker; x < m; x++) {
          Matr[y * m + x] =
              Matr[y * m + x] - (Div * Matr[grid_walker % my_part * m + x]);
        }
      }
    } else {
      MPI_Bcast(&Others[0], 1, row,
                (grid_walker >= n - (n % commsize)
                     ? commsize - 1
                     : get_rnk(commsize, n, grid_walker)),
                MPI_COMM_WORLD);
      if (my_end > grid_walker)
        for (size_t y = 0; y < my_part; y++) {
          Div = Matr[y * m + grid_walker] >= 1
                    ? Matr[y * m + grid_walker] / Others[grid_walker]
                    : Matr[y * m + grid_walker] * (1 / Others[grid_walker]);
          for (size_t x = grid_walker; x < m; x++) {
            Matr[y * m + x] = Matr[y * m + x] - (Div * Others[x]);
          }
        }
    }
  }
  for (int grid_back = n - 1; grid_back >= 1; grid_back--) {
    if (grid_back > my_begin && grid_back < my_end) {
      if (grid_back == my_begin) {
        int for_y = 1;
        Div = Matr[for_y * m + grid_back] * Matr[(for_y - 1) * m + grid_back];
        for (int x = my_part; x >= grid_back; x--) {
          Matr[(grid_back - 1) * m + x] -= Matr[for_y * m + x] * Div;
        }
        MPI_Bcast(&Matr[for_y], 1, row, rank, MPI_COMM_WORLD);
        for (int y = (for_y - 1); y >= 0; y--) {
          Div = Matr[y * m + grid_back] * Matr[for_y * m + grid_back];
          for (int x = m - 1; x >= grid_back; x--) {
            Matr[y * m + x] = Matr[y * m + x] - (Matr[y * m + x] * Div);
          }
        }
      } else {
        Div = Matr[grid_back % my_part * m + grid_back] *
              Matr[(grid_back % my_part - 1) * m + grid_back];
        for (int x = my_part; x >= grid_back; x--) {
          Matr[(grid_back % my_part - 1) * m + x] -=
              Matr[grid_back % my_part * m + x] * Div;
        }
        MPI_Bcast(&Matr[grid_back % my_part], 1, row, rank, MPI_COMM_WORLD);
        for (int y = (grid_back - 1) % my_part; y >= 0; y--) {
          Div = Matr[y * m + grid_back] *
                Matr[grid_back % my_part * m + grid_back];
          for (int x = m - 1; x >= grid_back; x--) {
            Matr[y * m + x] = Matr[y * m + x] - (Matr[y * m + x] * Div);
          }
        }
      }
    } else {
      MPI_Bcast(
          &Others[0], 1, row,
          (grid_back >= n - (n % commsize) ? commsize - 1
                                           : get_rnk(commsize, n, grid_back)),
          MPI_COMM_WORLD);
      if (my_begin < grid_back)
        for (int y = (grid_back - 1) % my_part; y >= 0; y--) {
          Div = Matr[y * m + grid_back] * Others[grid_back];
          for (int x = m - 1; x >= grid_back; x--) {
            Matr[y * m + x] = Matr[y * m + x] - (Others[x] * Div);
          }
        }
    }
  }
  MPI_Allreduce(&Others, MPI_IN_PLACE, 1, row, MPI_SUM, MPI_COMM_WORLD);
  rank == 0 ? printf("Time of work is: %f\n", MPI_Wtime() - now) : 0;
  free(Matr);
  free(Others);
  MPI_Finalize();
  return 0;
}
