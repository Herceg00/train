#include <stdlib.h>
#include "stdio.h"
#include "string.h"
#include "math.h"


short sign(double num) {
    if (num > 0) return 1;
    if (num < 0) return -1;
}

double cosinus_Fk(double x, double y) {
    if (y == 0) {
        return sqrt(2) / 2;
    } else {
        return sqrt(0.5 + abs(y) / (2 * sqrt(pow(x, 2) + pow(y, 2))));
    }
}

double sinus_Fk(double x, double y) {
    if (y == 0) {
        return sqrt(2) / 2;
    } else {
        return sign(x * y) * abs(x) / (2 * cosinus_Fk(x, y) * sqrt(pow(x, 2) + pow(y, 2)));
    }
}

double out_of_diag_sum(double** A,int matrix_size) {
    double sum = 0;
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            if (i != j) {
                sum += pow(abs(A[i][j]),2);
            }
        }
    }
}


void chosen_element(double ** A, short strategy,int matrix_size, int* i,int* j,int iteration){
    switch (strategy) {
        case 0 : { // Выбираем опорный элемент среди всех элементов матрицы
            double max_elem = 0;
            for (int x = 0; x < matrix_size; x++) {
                for (int y = 0; y < matrix_size; y++) {
                    if ((x != y) && (abs(A[x][y]) > abs(max_elem))) {
                        max_elem = A[x][y];
                        *i = x;
                        *j = y;
                    }
                }
            }
            break;
        }
        case 1: {  // Пронумеровали все элементы и соотнесли с итерацией цикла
            int number = pow(matrix_size,2) - matrix_size;
            iteration = iteration % number;
            int current_string = (iteration) / (matrix_size - 1);
            *i = current_string;

            int mod;
            mod = iteration % (matrix_size - 1);
            if(mod>=current_string) {
                *j = mod + 1;
            }else{
                *j = mod;
            }
            printf("%d,%d\n",*i,*j);
            break;

        }

        case 2: {  // Выбор максимального элемента из максимальной по модулю строчки
            int max_line = 0;
            double max_sum = 0;

            for (int string  = 1; string  < matrix_size; string++) {
                double line_sum = 0;
                for (int column = 0; column < matrix_size; column++) {
                    if (string != column) {
                        line_sum += pow(A[string][column],2);
                    }
                }
                if (line_sum >= max_sum) {
                    max_line = string;
                    max_sum = line_sum;
                }

            }


            double max_elem = 0;


            for (int column = 0; column < matrix_size;column++) {
                if ((column != max_line) && (abs(A[max_line][column]) > abs(max_elem))) {
                    max_elem = abs(A[max_line][column]);
                    *j = column;
                    *i = max_line;
                }
            }
            break;
        }
        default:{

        }
    }
}



void yakobi_rotation(int matrix_size, double** A,double** Eig, short strategy, double epsilon, double * eigenvalues) {
    int iteration = 0;
    while (out_of_diag_sum(A, matrix_size) > epsilon) {
        //iteration++;
        int i = 0;
        int j = 0;
        //int position = pow(matrix_size,2) - matrix_size;
        //int current_iterarion = (iteration % position) +1;
        chosen_element(A,strategy,matrix_size,&i,&j,iteration);

        double cos  = cosinus_Fk(A[i][j]*(-2),A[i][i] - A[j][j]);
        double sin = sinus_Fk(A[i][j]*(-2),A[i][i] - A[j][j]);
        //A[i][i]
        double Aii = pow(cos,2)*A[i][i] + pow(sin,2)*A[j][j] - 2*cos*sin*A[i][j];
        //A[j][j]
        double Ajj =  pow(cos,2)*A[j][j] + pow(sin,2)*A[i][i] + 2*cos*sin*A[i][j];

        for(int k = 0;k< matrix_size;k++){
            if( (k!=i) && (k!=j)) {
                double temp1 = A[i][k]*cos - sin*A[j][k];
                double temp2 = A[k][i]*sin + cos*A[k][j];
                A[i][k] = temp1;
                A[k][i] = temp1;
                A[k][j] = temp2;
                A[j][k] = temp2;

            }

        }

        A[i][i] = Aii;
        A[j][j] = Ajj;
        A[i][j] = 0;
        A[j][i] = 0;
        iteration++;
    }

    for(int i = 0; i < matrix_size;i++){
        eigenvalues[i] = A[i][i];
    }

}


double formula(int x, int y, int matrix_size){

}



int main(int argc, char** argv) {
    FILE * file = fopen("matrix.txt","r+");
    //int matrix_size = atoi(argv[1]);
    int matrix_size = 3;
    double **A = (double**)calloc(sizeof(double*),matrix_size);
    for (int i = 0; i < matrix_size; i++) {
        A[i] = (double*)calloc(sizeof(double),matrix_size);
    }


    double **Eig = (double**)calloc(sizeof(double*),matrix_size);
    for (int i = 0; i < matrix_size; i++) {
        Eig[i] = (double*)calloc(sizeof(double),matrix_size);
        Eig[i][i] = 1;
    }

    double* eigenvalues = (double*)calloc(sizeof(double),matrix_size);
    //char* input = argv[2];
    //char input = 'k';
    //if (strcmp(input, "formula") == 0) {
    if ( 3 == 0){
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                A[i][j] = formula(i, j, matrix_size);
            }
        }
    }
    else {
        for (int i = 0; i < matrix_size; i++) {
            for (int j = 0; j < matrix_size; j++) {
                fscanf(file, "%lf", &A[i][j]);
            }
        }
    }

    printf("ИСХОДНАЯ МАТРИЦА:\n");
    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("%lf  ",A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n\n");

    //short strategy = atoi(argv[3]);

    short strategy = 0;
    yakobi_rotation(matrix_size,A, Eig, strategy, 0.1, eigenvalues);


    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            printf("%lf  ",A[i][j]);
        }
        printf("\n");
    }
    printf("\n\n\nEIGENVALUES:\n");
    for (int j = 0; j < matrix_size; j++) {
        printf("%lf  ",eigenvalues[j]);
    }
    printf("\n\n");


    fclose(file);
    for (int k = 0; k < matrix_size; k++) {
        free(A[k]);
        free(Eig[k]);
    }
    free(Eig);
    free(eigenvalues);
    free(A);
}