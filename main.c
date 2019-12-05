#include </usr/include/openmpi/mpi.h>
#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "math.h"

#define MASTER 0
#define a_tag 0
#define b_tag 1
#define c_tag 2


void Multiply(int block_size, const int *A, const int *B, int *C) {
    for (int k = 0; k < block_size; k++) {
        for (int i = 0; i < block_size; i++) {
            int r = A[i * block_size + k];
            for (int j = 0; j < block_size; j++) {
                C[i * block_size + j] += r * B[k * block_size + j];
            }
        }
    }
}


int do_cannon(int n_threads, int blocks_per_dimension, int block_size,int*A,int*B,int*C) {

    int dims[2];
    int periods[2];
    int source_x_rank, dest_x_rank, temp_rank, dest_y_rank, source_y_rank;
    dims[0] = blocks_per_dimension;
    dims[1] = blocks_per_dimension;
    int coordinates[2];

    periods[0] = 1;
    periods[1] = 1;

    MPI_Comm MPI_COMM_GRID;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &MPI_COMM_GRID);
    MPI_Comm_rank(MPI_COMM_GRID, &temp_rank);
    MPI_Cart_coords(MPI_COMM_GRID, temp_rank, 2, coordinates);

    MPI_Cart_shift(MPI_COMM_GRID, 0, -1, &source_x_rank, &dest_x_rank);
    MPI_Cart_shift(MPI_COMM_GRID, 1, -1, &source_y_rank, &dest_y_rank);


    for (int i = 0; i < dims[0]; i++) {

        Multiply(block_size,A,B,C);
        MPI_Sendrecv_replace(A,block_size*block_size,MPI_INTEGER,dest_x_rank,1,source_x_rank,1,MPI_COMM_GRID,MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(B,block_size*block_size,MPI_INTEGER,dest_y_rank,1,source_y_rank,1,MPI_COMM_GRID,MPI_STATUS_IGNORE);

    }

    MPI_Comm_free(&MPI_COMM_GRID);

}



int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int rank;
    int matrix_size;
    int *A = NULL, *B = NULL, *C = NULL;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == MASTER) {

        FILE *filed_A = fopen(argv[1], "rb+");
        assert(filed_A != NULL);
        FILE *filed_B = fopen(argv[2], "rb+");
        assert(filed_B != NULL);
        int size_a, size_b;
        fread(&size_a, sizeof(int), 1, filed_A);
        printf("%d\n", size_a);
        fread(&size_b, sizeof(int), 1, filed_B);
        printf("%d\n", size_b);
        if (size_a != size_b) {
            return 1;
        }
        A = calloc(size_a * size_a, sizeof(int));
        B = calloc(size_b * size_b, sizeof(int));
        C = calloc(size_a*size_b,sizeof(int));
        fread(A,sizeof(int),size_a*size_a,filed_A);
        fread(B, sizeof(int),size_b*size_b,filed_B);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&matrix_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int number_of_blocks = matrix_size * matrix_size / size;
    int blocks_per_dimension = sqrt(number_of_blocks);
    int block_length = matrix_size / blocks_per_dimension;
    MPI_Datatype block;
    MPI_Type_vector(block_length, block_length, matrix_size, MPI_INTEGER, &block);
    MPI_Type_commit(&block);
    int *block_A;
    int *block_B;
    int *block_C;
    block_A = calloc(block_length * block_length, sizeof(int));
    block_B = calloc(block_length * block_length, sizeof(int));
    block_C = calloc(block_length * block_length, sizeof(int));
    assert(block_A != NULL);
    assert(block_B != NULL);
    int index;
    if (rank == MASTER) {
        for (int i = 0; i < number_of_blocks; i++) {
            index = (i/blocks_per_dimension)*block_length*matrix_size + (i%blocks_per_dimension) * block_length;

            MPI_Send(&A[index],
                     block_length * block_length, block, i, a_tag, MPI_COMM_WORLD);
            MPI_Send(&B[index],
                     block_length * block_length, block, i, b_tag, MPI_COMM_WORLD);
        }
    }
    MPI_Recv(block_A, block_length * block_length, block, MASTER, a_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(block_B, block_length * block_length, block, MASTER, b_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    do_cannon(size,blocks_per_dimension,block_length,block_A,block_B,block_C);


    MPI_Barrier(MPI_COMM_WORLD);

    if( rank != MASTER){
        MPI_Send(block_C,block_length*block_length,MPI_INTEGER,MASTER,c_tag,MPI_COMM_WORLD);

    } else {
        MPI_Sendrecv(block_C,block_length*block_length,MPI_INT,MASTER,0,C,1,block,MASTER,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        for( int i= 1;i < size;i++){
            index = (i/blocks_per_dimension)*block_length*matrix_size + (i%blocks_per_dimension) * block_length;

            MPI_Recv(&C[index],1,block,MASTER,c_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        }
    }
    if(rank == MASTER)
    {
        FILE* filed_C = fopen(argv[3], "wb+");
        fwrite(&matrix_size, sizeof(int32_t), 1, filed_C);
        fwrite(C, sizeof(int32_t), n * n, filed_C);
        free(A);
        free(B);
        free(C);
    }
    free(block_B);
    free(block_A);
    free(block_C);

    return MPI_Finalize();
}
