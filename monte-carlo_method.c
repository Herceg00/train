#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

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

    srand(RAND_MAX);
    omp_set_num_threads(P);
    double born,single_lifetime,sum_lifetime, monte_start,monte_end;
    int sum_reached_b = 0;

    monte_start = omp_get_wtime();
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
    printf("PROBABILITY OF REACHING B:%f\n",((float)sum_reached_b/(float)N));
    printf("AVERAGE LIFETIME:%lf\n",sum_lifetime/N);
    printf("METHOD LIFETIME:%lf\n",monte_end - monte_start);
    printf("SEGMENT: [%d,%d], INITIAL POSITION: %d ",a,b,x);
    printf("NUMBER OF ELEMENTS: %d, ATOMIC PROBABILITY: %f, NUMBER OF THREADS: %d\n",N,p,P);


}
