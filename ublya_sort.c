#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


/*int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)*/

// your code should be here.

// pseudocode from https://en.wikipedia.org/wiki/Odd%E2%80%93even_sort

/* function oddEvenSort(list) { */
/* 	function swap(list, i, j) { */
/* 		var temp = list[i]; */
/* 		list[i] = list[j]; */
/* 		list[j] = temp; */
/* 	} */

/* 	var sorted = false; */
/* 	while (!sorted) { */
/* 		sorted = true; */
/* 		for (var i = 1; i < list.length - 1; i += 2) { */
/* 			if(list[i] > list[i + 1]) { */
/* 				swap(list, i, i + 1); */
/* 				sorted = false; */
/* 			} */
/* 		} */
/* 		for (var i = 0; i < list.length - 1; i += 2) { */
/* 			if (list[i] > list[i + 1]) { */
/* 				swap(list, i, i + 1); */
/* 				sorted = false; */
/* 			} */
/* 		} */
/* 	} */
/* } */




#define MASTER 0

static int cmp(const void *_a, const void *_b)
{
    return *((int *)_a) - *((int *)_b);
}


static int seed(const int rank, const int size)
{
    int *seeds = NULL;
    if (rank == MASTER) {
        srand(time(NULL));
        seeds = malloc(size * sizeof(int));
        if (seeds == NULL) {
            return -1;
        }
        for (int i = 0; i < size; i++) {
            seeds[i] = rand();
        }
    }
    int seed;
    MPI_Scatter(seeds, 1, MPI_INTEGER, &seed, 1, MPI_INTEGER, MASTER, MPI_COMM_WORLD);
    if (rank == MASTER) {
        free(seeds);
    }
    srand(seed);
    return 0;
}


static int *gen_numbers(const int n)
{
    int *nums = malloc(n * sizeof(int));
    if (nums == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        nums[i] = rand() % 100500;
    }
    return nums;
}

static int is_sorted(const int *nums, const int n)
{
    for (int i = 0; i < n - 1; i++) {
        if (nums[i] > nums[i + 1]) {
            return 0;
        }
    }
    return 1;
}

static void min_max(const int *nums, const int n, int *min, int *max)
{
    *min = nums[0];
    *max = nums[0];
    for (int i = 1; i < n; i++) {
        if (nums[i] > *max) {
            *max = nums[i];
        }
        if (nums[i] < *min) {
            *min = nums[i];
        }
    }
}

static void merge(const int *a, const int *b, int *c,int size_a, int size_b){
    int i1=0,i2=0;
    for( int i =0;i<size_a+size_b;i++){
        if((i2 == size_b) || (i1!=size_a && a[i1] < b[i2])) {
            c[i] = a[i1++];
        }else{
            c[i] = b[i2++];
        }
    }
}


static void merge_two_parts(int *m, int *buffer, const int n) {
    int a = 0, b = n / 2;
    for (int i = 0; i < n; i++) {
        if (b == n || (a != n / 2 && m[a] < m[b])) {
            buffer[i] = m[a++];
        } else {
            buffer[i] = m[b++];
        }
    }

    for (int i = 0; i < n; i++) {
        m[i] = buffer[i];
    }
}


static void sort(const int rank, const int size, int *nums, const int n){
    int *local_buffer;
    local_buffer = malloc(sizeof(int) * n);
    for (int i=0;i<n;i++){
        local_buffer[i] = nums[i];
    }
    int *merge_buffer;
    int flag_to_change,test_element;
    merge_buffer = malloc(sizeof(int)*2*n);
    int *recv_buffer;
    recv_buffer = malloc(sizeof(int)*n);
    assert(recv_buffer);
    assert(merge_buffer);
    assert(local_buffer);
    qsort(local_buffer, n, sizeof(int), cmp); //Have sorted array on each process

    short sorted = 0;

    while() {
        if(rank % 2 == 1){
            MPI_Send(&local_buffer[n-1], 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD);
            MPI_Recv(&flag_to_change, 1, MPI_INTEGER, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(flag_to_change == 1){
                MPI_Send(local_buffer,n , MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD);
                MPI_Recv(local_buffer, n, MPI_INTEGER, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

        } else {
            MPI_Recv(&test_element, 1, MPI_INTEGER, rank, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(test_element > local_buffer[0]){
                flag_to_change = 1;
                MPI_Send(&flag_to_change, 1, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);

            } else{
                flag_to_change = 0;
                MPI_Send(&flag_to_change, 1, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);
            }

            if(flag_to_change == 1){
                MPI_Recv(recv_buffer, n, MPI_INTEGER, rank, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                merge(merge_buffer,local_buffer,recv_buffer,n,n);
                MPI_Send(merge_buffer, n, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank % 2 == 0){
            MPI_Send(&local_buffer[n-1], 1, MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD);
            MPI_Recv(&flag_to_change, 1, MPI_INTEGER, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(flag_to_change == 1){
                MPI_Send(local_buffer,n , MPI_INTEGER, rank+1, rank, MPI_COMM_WORLD);
                MPI_Recv(local_buffer, n, MPI_INTEGER, rank+1, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

        } else {

            MPI_Recv(&test_element, 1, MPI_INTEGER, rank, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(test_element > local_buffer[0]){
                flag_to_change = 1;
                MPI_Send(&flag_to_change, 1, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);

            } else{
                flag_to_change = 0;
                MPI_Send(&flag_to_change, 1, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);
            }

            if(flag_to_change == 1){
                MPI_Recv(recv_buffer, n, MPI_INTEGER, rank, rank-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                merge();
                MPI_Send(merge_buffer, n, MPI_INTEGER, rank-1, rank, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }
    MPI_Gather;


    free(recv_buffer);
    free(local_buffer);
    free(merge_buffer);
}



static int check(const int rank, const int size, const int *nums, const int n)
{
    if (!is_sorted(nums, n)) {
        return -1;
    }
    int pair[2];
    min_max(nums, n, pair, &pair[1]);
    int *buf = NULL;
    if (rank == MASTER) {
        buf = malloc(2 * size * sizeof(int));
        if (buf == NULL) {
            return -2;
        }
    }
    MPI_Gather(pair, 2, MPI_INTEGER, buf, 2, MPI_INTEGER, MASTER, MPI_COMM_WORLD);
    if (rank == MASTER) {
        int rc = is_sorted(buf, 2 * size);
        free(buf);
        return rc == 0;
    }
    return 0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 2) {
        return -1;
    }
    int n = atoi(argv[1]);
    assert(n > 0);
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rc = seed(rank, size);
    assert(rc == 0);

    int *nums = gen_numbers(n);
    assert(nums);

    sort(rank, size, nums, n);


    rc = check(rank, size, nums, n);
    assert(rc == 0);

    return MPI_Finalize();
}