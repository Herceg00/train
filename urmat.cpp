#include <iostream>
//#include "math.h"
#include <cmath>
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
#include "string.h"


double u(double x, double t) {
    return (t * t + 1) * x - 1 + 2 / (25 * M_PI * M_PI)
                                 * (1 - pow(M_E, -25 * M_PI * M_PI * t)) * cos(5 * M_PI / 2 * x);
}

double f(double x, double t) {
    return 0;
}

double u_0(double x,double t) {
    return x - 1.0;
}

double fi_1(double x,double t) {
    return 0.0;
}

double fi_2(double x,double t) {
    return t * t;
}

double psi_1(double x,double t) {
    return t * t + 1;
}

double psi_2(double x,double t) {
    return 0.0;
}


//prev_tau is used to calculate f(X_i,prev_tau)

double
get_next_step_def(double *prev_layer, double *cur_layer, double a, double tau, double h, double prev_tau, int h_steps) {
    for (int i = 1; i < h_steps ; i++) {
        double op1 = prev_layer[i] + f(i*h,prev_tau) * tau;
        double op2 = a*a*tau/(h*h);
        double op3 = prev_layer[i-1]+prev_layer[i+1] - 2*prev_layer[i];
        cur_layer[i] = op1 + op2*op3;
        printf("%d:  %f %f %f ___%f\n",i,op1,op2,op3,cur_layer[i]);
    }

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
        fi1[i] = fi_1(0,i*tau);
        fi2[i] = fi_2(1, i * tau); //массивы из краевых условий
    }

    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u_0(i * h, 0); //задаем нулевой уровень
    }

    //одинаковые значения для всей сетки
    double A_i = (a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = ((2 * tau*a * a) / (h * h)) + 1;

    //краевые условия
    double hi_1 = 0;
    double mu_1;
    double mu_2;
    double hi_2 = 0;

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[1] = 0;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        betta[1] = fi1[i];
        for (int j = 2; j < h_steps + 1 ; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] + ((grid[(i - 1) % 2][j]) + f(j * h, i * tau)*tau)) / temp;
        }
        bool grid_layer = i%2;
        grid[grid_layer][h_steps] = fi2[i];
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[i + 1];
        }
    }

    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);
}


void indefinite_scheme_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *psi1, *psi2;
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau); //массивы из краевых условий
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u(i * h, 0); //задаем нулевой уровень
    }

    //одинаковые значения для всей сетки
    double A_i = (a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = ((2 * tau*a * a) / (h * h)) + 1;

    //краевые условия
    double hi_1 = 1;
    double mu_1;
    double mu_2;
    double hi_2 = 1;

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[1] = hi_1;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        mu_1 = (psi1[i]*h*(-1));
        mu_2 = (psi2[i]*h*(-1));
        betta[1] = mu_1;
        for (int j = 2; j < h_steps + 1 ; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] + ((grid[(i - 1) % 2][j]) + f(j * h, i * tau)*tau)) / temp;
        }
        int grid_layer = i%2;
        grid[grid_layer][h_steps] = ((mu_2+hi_2*betta[h_steps])/(1 - hi_2*alpha[h_steps]));
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[i + 1];
        }
    }

    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);
}

void indefinite_scheme_first_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *fi1, *psi2;
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = u(0, i * tau);
        psi2[i] = psi_2(1, i * tau); //массивы из краевых условий
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u(i * h, 0); //задаем нулевой уровень
    }

    //одинаковые значения для всей сетки
    double A_i = (a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = ((2 * tau*a * a) / (h * h)) + 1;

    //краевые условия
    double hi_1 = 0;
    double mu_1;
    double mu_2;
    double hi_2 = 1;

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[1] = hi_1;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        mu_1 = fi1[i];
        mu_2 = (psi2[i]*h*(-1));
        betta[1] = mu_1;
        for (int j = 2; j < h_steps + 1 ; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] + ((grid[(i - 1) % 2][j]) + f(j * h, i * tau)*tau)) / temp;
        }
        int grid_layer = i%2;
        grid[grid_layer][h_steps] = ((mu_2+hi_2*betta[h_steps])/(1 - hi_2*alpha[h_steps]));
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[i + 1];
        }
    }
    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);

}

void indefinite_scheme_second_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *psi1, *fi2;
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        fi2[i] = u(1, i * tau); //массивы из краевых условий
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u(i * h, 0); //задаем нулевой уровень
    }

    //одинаковые значения для всей сетки
    double A_i = (a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = ((2 * tau*a * a) / (h * h)) + 1;

    //краевые условия
    double hi_1 = 1;
    double mu_1;
    double mu_2;
    double hi_2 = 0;

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[1] = 1;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        mu_1 = (psi1[i]*h*(-1));
        betta[1] = mu_1;
        for (int j = 2; j < h_steps + 1 ; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] + ((grid[(i - 1) % 2][j]) + f(j * h, i * tau)*tau)) / temp;
        }
        bool grid_layer = i%2;
        grid[grid_layer][h_steps] = fi2[i];
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[i + 1];
        }
    }
    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);

}

void definite_scheme_first_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps;
    double h = 1.0 / h_steps;
    double *fi1, *fi2;
    fi1 = (double *) calloc(tau_steps+1, sizeof(double));
    fi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        fi1[i] = fi_1(0, i * tau);
        fi2[i] = fi_2(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps+1; i++) { //начальная инициализация нулевого слоя
        grid[0][i] = u_0(i*h,0);
    }
    for (int i = 1; i < tau_steps+1; i++) { //делаем i - ый слой по нижележащему
        get_next_step_def(grid[(i + 1)%2], grid[i%2], a, tau, h, (i) * tau, h_steps);
        grid[i%2][0] = fi1[i];
        grid[i%2][h_steps] = fi2[i];
    }

    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);

    printf("%lf",sum);
}



void definite_scheme_second_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *psi2;
    psi1 = (double *) calloc(tau_steps+1, sizeof(double));
    psi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        psi1[i] = psi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u_0(i * h, 0);
    }
    for (int i = 1; i < tau_steps+1; i++) {
        get_next_step_def(grid[(i + 1)%2], grid[i%2], a, tau, h, (i) * tau, h_steps);
        grid[i%2][0] = (psi1[i] * 2 * h + grid[i%2][2] - 4 * grid[i%2][1]) / (-3);
        grid[i%2][h_steps] = (psi2[i] * 2 * h - grid[i%2][h_steps - 2] + 4 * grid[i%2][h_steps - 1]) / 3;
    }
    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);
}



void definite_scheme_first_second(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *fi1, *psi2;
    fi1 = (double *) calloc(tau_steps+1, sizeof(double));
    psi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        fi1[i] = fi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps+1; i++) {
        grid[0][i] = u_0(i * h, 0);
    }
    for (int i = 1; i < tau_steps+1; i++) {
        get_next_step_def(grid[(i - 1)%2], grid[i%2], a, tau, h, (i - 1) * tau, h_steps);
        grid[i%2][0] = fi1[i];
        grid[i%2][h_steps] = (psi2[i] * 2 * h - grid[i%2][h_steps - 2] + 4 * grid[i%2][h_steps - 1]) / 3;
    }
    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);
}



void definite_scheme_second_first(int tau_steps, int h_steps, double a) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *fi2;
    psi1 = (double *) calloc(tau_steps+1, sizeof(double));
    fi2 = (double *) calloc(tau_steps+1, sizeof(double));
    for (int i = 0; i < tau_steps+1; i++) {
        psi1[i] = psi_1(0, i * tau);
        fi2[i] = fi_2(1, i * tau);
    }
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps+1, sizeof(double));
    grid[1] = (double *) calloc(h_steps+1, sizeof(double));
    for (int i = 0; i < h_steps+1; i++) {
        grid[0][i] = u_0(i * h, 0);
    }
    for (int i = 1; i < tau_steps+1; i++) {
        get_next_step_def(grid[(i - 1)%2], grid[i%2], a, tau, h, (i - 1) * tau, h_steps);
        grid[i%2][0] = (psi1[i] * 2 * h + grid[i%2][2] - 4 * grid[i%2][1]) / (-3);
        grid[i%2][h_steps] = fi2[i];
    }
    double sum =0;
    for(int i = 0;i<h_steps+1;i++){
        double diff = u(i*h,1.0) - grid[tau_steps%2][i];
        sum+=(diff*diff);
    }
    sum = sum*h;
    sum = sqrt(sum);
    printf("%lf",sum);
}

int main(int argc, char **argv) {
    int tau_steps = atoi(argv[1]); //чему равно N, то есть всего точек, начиная с 0, N+1
    int h_steps = atoi(argv[2]); //чему равно M, то есть всего точек, начиная с 0, M+1
    //doing only test mode now
    //definite_scheme_first_first(tau_steps,h_steps,2);
    //definite_scheme_second_second(tau_steps,h_steps,2);
    //definite_scheme_first_second(tau_steps,h_steps,2);
    //definite_scheme_second_first(tau_steps,h_steps,2);
    //indefinite_scheme_first_first(tau_steps,h_steps,2);
    //indefinite_scheme_second(tau_steps,h_steps,2);
    //indefinite_scheme_first_second(tau_steps,h_steps,2);
    indefinite_scheme_second_first(tau_steps,h_steps,2);

}





