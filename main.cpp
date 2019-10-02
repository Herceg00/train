#include "stdlib.h"
#include "stdio.h"
#include "string.h"


int main(int argc, char **argv) {
    int matrix_size = 4;
    double step_max = 0;

    //int matrix_size = atoi(argv[1]);
    double a[matrix_size][matrix_size] = {{3,4,1,2},
                                          {2,3,4,1},
                                          {1,2,3,4},
                                          {4,1,2,3}


    };
    int stringConseq[4] = {0, 1, 2,3};



    //OUTPUT THE MATRIX*
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }

    //ОСТАВЛЯЕМ ПЕРВЫЙ СТОЛБЕЦ КАК ЕСТЬ, ПОИЩЕМ У НЕГО МАКСИМУМ
    double first_max = a[0][0];
    int first_index = 0;
    for (int i = 1; i < matrix_size; i++) {
        if (a[i][0] > first_max) {
            first_max = a[i][0];
            first_index = i;
        }
    }
    stringConseq[0] = first_index;
    stringConseq[first_index] = 0;
    step_max = a[first_index][0];

    printf("%fkek\n", step_max);

    for (int j = 0; j < matrix_size; j++) { //ЦИКЛ СТРОКА-СТОЛБЕЦ-МАКСИМУМ

        for (int k = j + 1; k < matrix_size; k++) { //ДЕЛАЕМ СТРОКУ
            for (int l = 0; l < j; l++) {
                a[stringConseq[j]][k] -= a[stringConseq[j]][l] * a[stringConseq[l]][k];
            }
            a[stringConseq[j]][k] /= step_max;
        }

        //ДЕЛАЕМ СТОЛБЕЦ
        for (int k = j + 1; k < matrix_size; k++) {
            for (int l = 0; l < j + 1; l++) {
                a[stringConseq[k]][j + 1] -= a[stringConseq[k]][l] * a[stringConseq[l]][j+1];
                printf("%f *** %f  ",a[stringConseq[k]][l] , a[stringConseq[l]][j+1]);
            }
        }

        if (j != matrix_size-1) {
            step_max = a[stringConseq[j + 1]][j + 1];

        }
        int max_index = j + 1;
        for (int i = j + 2; i < matrix_size; i++) {
            if (a[stringConseq[i]][j + 1] > step_max) {
                step_max = a[stringConseq[i]][j + 1];
                max_index = i;
            }
        }

        //МЕНЯЕМ ЭЛЕМЕНТЫ ВЕКТОРА-УКАЗАТЕЛЯ СТРОКИ

        int c = stringConseq[j + 1];
        stringConseq[j + 1] = stringConseq[max_index];
        stringConseq[max_index] = c;

        printf("%f ", step_max);

    }






    printf("STRING CONSEQ:  ");

    for (int k = 0; k < matrix_size; k++) {
        printf("%d", stringConseq[k]);
    }

    printf("\n\n");

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }


}

