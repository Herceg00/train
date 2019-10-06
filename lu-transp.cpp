#include <stdlib.h>
#include "stdio.h"
#include "string.h"
#include <cmath>

double min(double x, double y) {
    if (x < y) {
        return x;
    }
    return y;
}

short decision(int n, float **a, float *b, float *x, int *stringConseq) {
    int matrix_size = n;

    float *y = (float *) calloc(sizeof(float), matrix_size);
    //ЗАПОЛНИЛИ STRING CONSEQ
    for (int i = 0; i < matrix_size; i++) {
        stringConseq[i] = i;
    }
    float step_max = 0;
    float first_max = a[0][0];
    int first_index = 0;
    for (int i = 1; i < matrix_size; i++) {
        if (abs(a[i][0]) > abs(first_max)) {
            first_max = a[i][0];
            first_index = i;
        }
    }
    if (first_max == 0) {
        free(y);
        return 0;
    }
    stringConseq[0] = first_index;
    stringConseq[first_index] = 0;
    step_max = a[first_index][0];
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
                a[stringConseq[k]][j + 1] -= a[stringConseq[k]][l] * a[stringConseq[l]][j + 1];
            }
        }

        if (j != matrix_size - 1) {
            step_max = a[stringConseq[j + 1]][j + 1];

        }
        int max_index = j + 1;
        for (int i = j + 2; i < matrix_size; i++) {
            if (abs(a[stringConseq[i]][j + 1]) > abs(step_max)) {
                step_max = a[stringConseq[i]][j + 1];
                max_index = i;
            }
        }

        //МЕНЯЕМ ЭЛЕМЕНТЫ ВЕКТОРА-УКАЗАТЕЛЯ СТРОКИ
        if (j != matrix_size - 1) {
            int c = stringConseq[j + 1];
            stringConseq[j + 1] = stringConseq[max_index];
            stringConseq[max_index] = c;
        }
        if (step_max == 0) {
            free(y);
            return 0;
        }
    }



    //РЕШАЕМ Ly = b//

    for (int i = 0; i < matrix_size; i++) {
        y[stringConseq[i]] = b[stringConseq[i]];
        for (int j = 0; j < i; j++) {
            y[stringConseq[i]] -= y[stringConseq[j]] * a[stringConseq[i]][j];
        }
        y[stringConseq[i]] /= a[stringConseq[i]][i];
    }
    printf("\n\n");

    //РЕШАЕМ Ux = y//

    for (int i = matrix_size - 1; i > -1; i--) {
        x[stringConseq[i]] = y[stringConseq[i]];
        for (int j = i + 1; j < matrix_size; j++) {
            x[stringConseq[i]] -= a[stringConseq[i]][j] * x[stringConseq[j]];
        }
    }
    free(y);
}

void output(float **a, int *stringConseq, int matrix_size) {
    printf("\nМАТРИЦА L:\n");
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < i + 1; j++) {
            printf("%f ", a[stringConseq[i]][j]);
        }
        for (int j = i + 1; j < matrix_size; j++) {
            printf("%f ", 0.00000);
        }
        printf("\n");
    }

    printf("\nМАТРИЦА U:\n");
    //ВЫВОД U
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < i; j++) {
            printf("%f ", 0.00000);
        }
        printf("%f ", 1.00000);
        for (int j = i + 1; j < matrix_size; j++) {
            printf("%f ", a[stringConseq[i]][j]);
        }
        printf("\n");
    }
}

void file_output(float **a, int *stringConseq, int matrix_size, FILE *filestream) {
    fprintf(filestream, "\nМАТРИЦА L:\n");
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < i + 1; j++) {
            fprintf(filestream, "%f ", a[stringConseq[i]][j]);
        }
        for (int j = i + 1; j < matrix_size; j++) {
            fprintf(filestream, "%f ", 0.00000);
        }
        fprintf(filestream, "\n");
    }

    fprintf(filestream, "\nМАТРИЦА U:\n");
    //ВЫВОД U
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < i; j++) {
            fprintf(filestream, "%f ", 0.00000);
        }
        fprintf(filestream, "%f ", 1.00000);
        for (int j = i + 1; j < matrix_size; j++) {
            fprintf(filestream, "%f ", a[stringConseq[i]][j]);
        }
        fprintf(filestream, "\n");
    }


}

float formula(int x, int y, int z) {
    return 2 * x + 4 * y - z;
}


float discrepancy(float **a, int *stringConseq, int matrix_size, float *b, float *x) {
    for (int i = matrix_size - 1; i > -1; i--) {
        for (int j = matrix_size - 1; j > -1; j--) {
            double sum = 0.00, mul_1, mul_2, m = min(i, j);
            for (int k = 0; k <= m; k++) {
                mul_1 = a[stringConseq[i]][k];
                mul_2 = a[stringConseq[k]][j];
                if (k == j) {
                    mul_2 = 1.00;
                }
                sum += mul_1 * mul_2;
            }
            a[stringConseq[i]][j] = sum;
        }
    }
    float discrepancy = 0.00, sum = 0.00;
    printf("\n\n\n\n");

    for (int i = 0; i < matrix_size; i++) {
        sum = 0.00;
        for (int j = 0; j < matrix_size; j++) {
            sum += (a[i][j] * x[stringConseq[j]]);
        }
        discrepancy += (sum - b[i]);
        printf("\n");
    }
    return discrepancy;

};


int main(int argc, char **argv) {
    FILE *file = fopen(argv[1], "r");
    char *mode = argv[2];
    int matrix_size;
    short error;
    fscanf(file, "%d", &matrix_size);
    float step_max = 0;
    float *b = (float *) calloc(matrix_size, sizeof(float));
    float **a = (float **) calloc(sizeof(matrix_size) * sizeof(float), matrix_size);

    // ВЫДЕЛИЛИ ПАМЯТЬ НА СТРОКУ
    for (int i = 0; i < matrix_size; i++) {
        a[i] = (float *) calloc(sizeof(float), matrix_size);
    }

    //ВЫДЕЛИЛИ ПАМЯТЬ НА РЕШЕНИЕ
    float *x = (float *) calloc(sizeof(float), matrix_size);
    if (strcmp(mode, "formula") == 0) {
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                a[i][j] = formula(i, j, matrix_size);
            }
        }

    } else {
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                fscanf(file, "%f", &a[i][j]);
            }
            fscanf(file, "%f", &b[i]);
        }

    }

    int *string_row = (int *) calloc(sizeof(int), matrix_size);
    error = decision(matrix_size, a, b, x, string_row);
    if (error == 0) {
        for (int k = 0; k < matrix_size; k++) {
            free(a[k]);
        }
        free(string_row);
        free(a);
        free(b);
        free(x);
        fclose(file);
        printf("ВЫРОЖДЕННАЯ МАТРИЦА");
        return 0;
    }
    printf("\n*********РЕШЕНИЕ[x]:**********\n");
    for (int k = 0; k < matrix_size; k++) {
        printf("%f  ", x[string_row[k]]);
    }
    if (matrix_size < 10) {
        output(a, string_row, matrix_size);

    } else {
        FILE *file1 = fopen("output.txt", "w+");
        file_output(a, string_row, matrix_size, file1);
        fclose(file1);
    }

    float disc = discrepancy(a, string_row, matrix_size, b, x);
    printf("\n*********НЕВЯЗКА:**********\n%f\n", disc);
    for (int k = 0; k < matrix_size; k++) {
        free(a[k]);
    }
    free(string_row);
    free(a);
    free(b);
    free(x);
    fclose(file);


}

