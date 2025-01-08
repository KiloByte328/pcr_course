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
  int m = n + 10;
  int my_begin = rank == 0 ? 0 : rank * (n / commsize);
  int my_end = rank == commsize - 1 ? n : (rank + 1) * (n / commsize);
  int my_part = my_end - my_begin;
  double* Matr = malloc(sizeof(double) * (m * my_part));
  double* Others = malloc(sizeof(double) * m);
  double* Res = malloc(sizeof(double) * n * (m - n));
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
      rank == commsize - 1 ? printf("i want to modyfi: %f, i get this to modyfi: %f\n", Matr[grid_walker % my_part + grid_walker], Others[grid_walker]) : 0;
      if (my_end > grid_walker)
        for (size_t y = 0; y < my_part; y++) {
          Div = Matr[y * m + grid_walker] >= 1 ? Matr[y * m + grid_walker] / Others[grid_walker] : Matr[y * m + grid_walker] * ( 1 / Others[grid_walker]);
          for (size_t x = grid_walker; x < m; x++) {
            Matr[y * m + x] = Matr[y * m + x] - (Div * Others[x]);
          }
        }
    }
    rank == commsize - 1 ? printf("i modyfied it: %f, grid_walker is: %d\n", Matr[grid_walker % my_part + grid_walker], grid_walker) : 0;
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
      rank == commsize - 1 ? printf("grid is: %d, i get: %f, i want modyfi: %f\n", grid_back, Others[grid_back], Matr[grid_back % my_part * m + grid_back]): 0;
      if (my_begin < grid_back)
        for (int y = (grid_back - 1) % my_part; y >= 0; y--) {
          Div = Matr[y * m + grid_back] * Others[grid_back];
          for (int x = m - 1; x >= grid_back; x--) {
            Matr[y * m + x] = Matr[y * m + x] - (Others[x] * Div);
          }
        } 
    }
  }
  //printf("my part is: %d, my start is: %d, my end is: %d \n", my_part, my_begin, my_end);
  for(int y = 0; y < n; y++)
  {
	  if (y >= my_begin && y < my_end)
		MPI_Bcast(&Matr[y % my_part * n + n - 1], m - n - 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
	  else
		MPI_Bcast(&Res[y * (m - n)],  m - n - 1, MPI_DOUBLE, (y >= n - n % commsize ? commsize - 1: get_rnk(commsize, n, y)), MPI_COMM_WORLD);
	  // поставь сюда дома Allreduce вместо бкастов
  }
  for (int y = 0; y < my_part; y++)
  {
	for(int x = 0; x < m - n; x++)
	{
		Res[my_begin + y + x] = Matr[y * m + x];
		rank == 1 ? printf("Copy form y = %d, x = %d, what copied: %f\n", y, x, Matr[y * m + x]) : 0;
	}
  }
  MPI_Barrier(MPI_COMM_WORLD);
  rank == 0 ? printf("Time of work is: %f\n", MPI_Wtime() - now) : 0;
  free(Matr);
  free(Res);
  free(Others);
  MPI_Finalize();
  return 0;
}
