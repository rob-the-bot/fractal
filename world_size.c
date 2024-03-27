#include <mpi.h>
#include <stdio.h>

int main(int argc, char const *argv[])
{
	MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    printf("%d\n", world_size);
    MPI_Finalize();
}
