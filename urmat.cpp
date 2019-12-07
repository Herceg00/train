#include <iostream>
#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
#include "string.h"
#include "cublas.h"
#include "/usr/local/cuda-10.1/include/cuda_runtime.h"


double u(double x, double t) {
    return 0;
    //tutor function

}

double f(double x, double t) {
    return 0;
    //tutor function
}

double du_dx(double x, double t) {
    return 0;
    //tutor function
}

//prev_tau is used to calculate f(X_i,prev_tau)

double
get_next_step_def(double *prev_layer, double *cur_layer, double a, double tau, double h, double prev_tau, int h_steps) {
    for (int i = 1; i < h_steps - 1; i++) {

        cur_layer[i] = (tau * a * a * (prev_layer[i - 1] - 2 * prev_layer[i] + prev_layer[i + 1])) / (h * h) +
                       f(i * h, prev_tau) + prev_layer[i];
    }

}


int main(int argc, char **argv) {
    int tau_steps = atoi(argv[1]); //чему равно N, то есть всего точек, начиная с 0, N+1
    int h_steps = atoi(argv[2]); //чему равно M, то есть всего точек, начиная с 0, M+1
    //doing only test mode now


}

//h_steps + 1 - количество точек на сетке по координате OX
//tau_steps + 1 - количество точек на сетке по координате OY

//неявная схема 1_1
void indefinite_scheme_first_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *fi1, *fi2;
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = u(0, i * tau);
        fi2[i] = u(1, i * tau); //массивы из краевых условий
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u(i * h, 0); //задаем нулевой уровень
    }

    //одинаковые значения для всей сетки
    double A_i = a * a / h * h;
    double B_i = a * a / h * h;
    double C_i = 2 * a * a / h * h + 1 / tau;

    //краевые условия
    double hi_1 = 0;
    double mu_1;
    double mu_2;
    double hi_2 = 0;

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[0] = 0;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        betta[0] = fi1[i];
        for (int j = 1; j < h_steps + 1 ; j++) { //заполняем прогоночные коэффициенты
            alpha[j] = B_i / (C_i - alpha[j - 1] * A_i);
            betta[j] = (A_i * betta[j - 1] + ((grid[(i - 1) % 2][j]) / tau + f(j * h, i * tau))) /
                       (C_i - alpha[i - 1] * A_i);
        }
        bool grid_layer = i%2;
        grid[grid_layer][h_steps] = fi2[i];
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[i + 1];
        }
    }
}

void indefinite_scheme_second(int tau_steps, int h_steps, double a) {

}

void indefinite_scheme_first_second(int tau_steps, int h_steps, double a) {

}

void indefinite_scheme_second_first(int tau_steps, int h_steps, double a) {


}


void definite_scheme_first_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps;
    double h = 1.0 / h_steps;
    double *fi1, *fi2;
    fi1 = (double *) calloc(tau_steps+1, sizeof(double));
    fi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        fi1[i] = u(0, i * tau);
        fi2[i] = u(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps+1; i++) { //начальная инициализация нулевого слоя
        grid[0][i] = u(i * h, 0);
    }
    for (int i = 1; i < tau_steps+1; i++) { //делаем i - ый слой по нижележащему
        get_next_step_def(grid[(i - 1)%2], grid[i%2], a, tau, h, (i - 1) * tau, h_steps);
        grid[i%2][0] = fi1[i];
        grid[i%2][h_steps] = fi2[i];
    }
}

void definite_scheme_second_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *psi2;
    psi1 = (double *) calloc(tau_steps+1, sizeof(double));
    psi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        psi1[i] = du_dx(0, i * tau);
        psi2[i] = du_dx(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps; i++) {
        grid[0][i] = u(i * h, 0);
    }
    for (int i = 1; i < tau_steps+1; i++) {
        get_next_step_def(grid[(i - 1)%2], grid[i%2], a, tau, h, (i - 1) * tau, h_steps);
        grid[i%2][0] = (psi1[i] * 2 * h + grid[i%2][2] - 4 * grid[i%2][1]) / (-3);
        grid[i%2][h_steps - 1] = (psi2[i] * 2 * h - grid[i%2][h_steps - 3] + 4 * grid[i%2][h_steps - 2]) / 3;
    }

}

void definite_scheme_first_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *fi1, *psi2;
    fi1 = (double *) calloc(tau_steps, sizeof(double));
    psi2 = (double *) calloc(tau_steps, sizeof(double));
    for (int i = 0; i < tau_steps; i++) {
        fi1[i] = u(0, i * tau);
        psi2[i] = du_dx(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps, sizeof(double));
    grid[1] = (double *) calloc(h_steps, sizeof(double));
    for (int i = 0; i < h_steps; i++) {
        grid[0][i] = u(i * h, 0);
    }
    for (int i = 1; i < tau_steps; i++) {
        get_next_step_def(grid[i - 1], grid[i], a, tau, h, (i - 1) * tau, h_steps);
        grid[i][0] = fi1[i];
        grid[i][h_steps - 1] = (psi2[i] * 2 * h - grid[i][h_steps - 3] + 4 * grid[i][h_steps - 2]) / 3;
    }
}

void definite_scheme_second_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *fi2;
    psi1 = (double *) calloc(tau_steps, sizeof(double));
    fi2 = (double *) calloc(tau_steps, sizeof(double));
    for (int i = 0; i < tau_steps; i++) {
        psi1[i] = du_dx(0, i * tau);
        fi2[i] = u(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps, sizeof(double));
    grid[1] = (double *) calloc(h_steps, sizeof(double));
    for (int i = 0; i < h_steps; i++) {
        grid[0][i] = u(i * h, 0);
    }
    for (int i = 1; i < tau_steps; i++) {
        get_next_step_def(grid[i - 1], grid[i], a, tau, h, (i - 1) * tau, h_steps);
        grid[i][0] = (psi1[i] * 2 * h + grid[i][2] - 4 * grid[i][1]) / (-3);
        grid[i][h_steps - 1] = fi2[i];
    }
}





