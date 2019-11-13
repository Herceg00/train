#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "time.h"
#include "sys/time.h"


int main(int argc, char **argv) {
    if (argc < 7) {
        printf(" Not enough arguments");
        return -1;
    }
    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    if (a >= b) {
        printf("Wrong segment");
        return -1;
    }
    int x = atoi(argv[3]);
    int N = atoi(argv[4]);
    double p = strtod(argv[5], NULL);
    int P = atoi(argv[6]);
    struct timeval monte_start, monte_stop;
    omp_set_num_threads(P);
    double single_lifetime = 0, sum_lifetime = 0,
    time_init = time(NULL);
    int sum_reached_b = 0;
    gettimeofday(&monte_start,NULL);
    int iters = 0;
#pragma omp parallel default(none) shared(time_init) firstprivate(a, b, N, x, P, p) reduction (+:sum_reached_b,iters)
    {
        int seed_for_thread = omp_get_thread_num()+time_init;
#pragma omp for
        for (int i = 0; i < N; i++) {
            int position = x;
            while ((position != a) && (position != b)) {
                iters++;
                int num = rand_r(&seed_for_thread);
                if ((float)num / (float)RAND_MAX > p)
                    position--;
                else
                    position++;
            }
                if (position == b) {
                    sum_reached_b++;
                }
        }
    }
    gettimeofday(&monte_stop,NULL);
    FILE * stats = fopen("stats.txt","a");
    if(stats == NULL){
        stats = stdout;
    }
    fprintf(stats,"PROBABILITY OF REACHING B:%f\n",(float)sum_reached_b/(float)N);
    fprintf(stats, "AVERAGE LIFETIME:%lf\n",(float)(iters*1.0)/(float)N);
    fprintf(stats,"METHOD LIFETIME:%ld\n",(monte_stop.tv_sec - monte_start.tv_sec)*1000000 + monte_stop.tv_usec - monte_start.tv_usec);
    fprintf(stats, "SEGMENT: [%d,%d], INITIAL POSITION: %d ",a,b,x);
    fprintf(stats,"NUMBER OF ELEMENTS: %d, ATOMIC PROBABILITY: %f, NUMBER OF THREADS: %d\n\n",N,p,P);
    fclose(stats);
}
