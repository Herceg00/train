#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "assert.h"

int main(int argc, char **argv) {

    if (argc < 7) {
        printf("Not enough arguments");
        return -1;
    }
    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int x = atoi(argv[3]);
    int N = atoi(argv[4]);
    double p = strtod(argv[5], NULL);
    int P = atoi(argv[6]);

    srand(RAND_MAX);
    omp_set_num_threads(P);
    printf("%d\n",x);
    double born,lifetime;


#pragma omp parallel for schedule(dynamic,1) shared(a,b,x,p,N)
    for (int i = 0; i < N; i++) {
        int position  = x;
        born = omp_get_wtime();
        while((position!=b) && (position!=a)) {
            if ((rand() % 100) / 100.0 >= p) {
                position--;
            } else {
                position++;
            }
        }
        lifetime = omp_get_wtime() - born;

#pragma omp critical
        {
            x += i;
            printf("running%d______%d___%d___%d\n", omp_get_thread_num(), i, x,position);

        }
    }


}
