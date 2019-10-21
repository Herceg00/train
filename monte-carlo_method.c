#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "time.h"



int main(int argc, char **argv) {

    if (argc < 7) {
        printf("Not enough arguments");
        return -1;
    }
    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    if( a >= b){
        printf("Wrong segment");
        return -1;
    }
    int x = atoi(argv[3]);
    int N = atoi(argv[4]);
    double p = strtod(argv[5], NULL);
    int P = atoi(argv[6]);

    //srand(RAND_MAX);
    omp_set_num_threads(P);
    double born=0,single_lifetime=0,sum_lifetime=0, monte_start=0,monte_end=0;
    int sum_reached_b = 0;
    int *seeds =  (int *)calloc(P, sizeof(int));
    for(int i=0;i<P;i++){
        seeds[i] = omp_get_wtime();
    }

    monte_start = omp_get_wtime();
#pragma omp parallel for schedule(dynamic,1) shared(a,b,x,p,N)
    for (int i = 0; i < N; i++) {
        int position  = x;
        //srand(seeds[omp_get_thread_num()]);
        born = omp_get_wtime();
        while((position!=b) && (position!=a)) {


            if ((rand_r(&seeds[omp_get_thread_num()]) % 1000) / 1000.0 >= p) {
                position--;
            } else {
                position++;
            }
        }
        single_lifetime = omp_get_wtime() - born;

#pragma omp critical
        {
            sum_lifetime +=single_lifetime;
            if (position == b){
               sum_reached_b++;
            }

        }
    }
    monte_end = omp_get_wtime();
    FILE * stats = fopen("stats.txt","r+");
    if(stats == NULL){
        stats = stdout;
    }
    fprintf(stats,"PROBABILITY OF REACHING B:%f\n",((float)sum_reached_b/(float)N));
    fprintf(stats, "AVERAGE LIFETIME:%lf\n",sum_lifetime/N);
    fprintf(stats,"METHOD LIFETIME:%lf\n",monte_end - monte_start);
    fprintf(stats, "SEGMENT: [%d,%d], INITIAL POSITION: %d ",a,b,x);
    fprintf(stats,"NUMBER OF ELEMENTS: %d, ATOMIC PROBABILITY: %f, NUMBER OF THREADS: %d\n",N,p,P);
    fclose(stats);
    free(seeds);


}