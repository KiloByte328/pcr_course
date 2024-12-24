#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
  size_t n = 1000;
  size_t m = n + 1;
  double* Matr = malloc(sizeof(double) * (m * n));
  for (size_t y = 0; y < n; y++) {
    for (size_t x = 0; x < m; x++) {
      Matr[y * m + x] = rand() % 100;
    }
  }
  double now = MPI_Wtime();
  double Div;
  for (size_t grid_walker = 0; grid_walker < n; grid_walker++) {
    Div = Matr[grid_walker * m + grid_walker];
    if (Matr[grid_walker * m + grid_walker] == 0.0)
      return 1;
    if (Matr[grid_walker * m + grid_walker] != 1) {
      for (size_t x = grid_walker; x < m; x++) {
        Matr[grid_walker * m + x] =
            Matr[grid_walker * m + grid_walker] >= 1
                ? (Matr[grid_walker * m + x] / Div)
                : (Matr[grid_walker * m + x] * (1 / Div));
      }
    }
    for (size_t y = grid_walker + 1; y < n; y++) {
      Div =
          Matr[y * m + grid_walker] >= 1
              ? Matr[y * m + grid_walker] / Matr[grid_walker * m + grid_walker]
              : Matr[y * m + grid_walker] *
                    (1 / Matr[grid_walker * m + grid_walker]);
      for (size_t x = grid_walker; x < m; x++) {
        Matr[y * m + x] = Matr[y * m + x] - (Div * Matr[grid_walker * m + x]);
      }
    }
  }
  for (size_t grid_back = n - 1; grid_back >= 1; grid_back--) {
    for (size_t y = grid_back - 1; y > 0; y--) {
      Div = Matr[y * m + grid_back] * Matr[grid_back * m + grid_back];
      for (size_t x = m - 1; x >= grid_back; x--) {
        Matr[y * m + x] = Matr[y * m + x] - (Matr[grid_back * m + x] * Div);
      }
    }
  }
  double end = MPI_Wtime() - now;
  printf("Time of work is: %f \n", end);
  return 0;
}
