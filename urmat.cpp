#include <iostream>
#include <cmath>
#include "stdio.h"
#include "assert.h"
#include <fstream>
#include "stdlib.h"
#include "string.h"



double a() {
    return 1.0;
};

double u(double x, double t) {
    return (1 + (t * t)) * sin(M_PI * x) * cos(3 * M_PI * x);
    //return (2 - cos(t)) * (x * x * x - x * x + 1);
}

double f(double x, double t) {
    return sin(M_PI * x) * cos(3 * M_PI * x) * 2 * t + (1 + t * t) * (10 * M_PI * M_PI * sin(M_PI * x) * cos(3 * M_PI * x) +
    6 * M_PI * M_PI * cos(M_PI * x) * sin(3 * M_PI * x));
    //return (sin(t) * (x * x * x - x * x + 1) - (2 - cos(t)) * (6 * x - 2));
}

double u_0(double x, double t) {
    return sin(M_PI * x) * cos(3 * M_PI * x);
    //return x * x * x - x * x + 1;
}

double fi_1(double x, double t) {
    return 0.0;
    //return 2 - cos(t);
}

double fi_2(double x, double t) {
    return 0.0;
    //return 2 - cos(t);
}

double psi_1(double x, double t) {
    return M_PI * (1 + t * t);
    //return 0;
}

double psi_2(double x, double t) {
    return M_PI * (1 + t * t);
    //return 0;
}


//prev_tau is used to calculate f(X_i,prev_tau)

void
get_next_step_def(double *prev_layer, double *cur_layer, double a, double tau, double h, double prev_tau, int h_steps) {
    for (int i = 1; i < h_steps; i++) {
        double op1 = prev_layer[i] + f(i * h, prev_tau) * tau;
        double op2 = a * a * tau / (h * h);
        double op3 = prev_layer[i - 1] + prev_layer[i + 1] - 2 * prev_layer[i];
        cur_layer[i] = op1 + op2 * op3;
    }
}


//h_steps + 1 - количество точек на сетке по координате OX
//tau_steps + 1 - количество точек на сетке по координате OY

void indefinite_scheme_first_first(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *fi1, *fi2;
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = fi_1(0, i * tau);
        fi2[i] = fi_2(1, i * tau); //массивы из краевых условий
    }
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    //одинаковые значения для всей сетки
    double A_i = ((-1) * a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = (((2 * tau * a * a) / (h * h)) + 1) * (-1);
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
        for (int j = 2; j < h_steps + 1; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            double F = grid[(i - 1) % 2][j - 1] + f((j - 1) * h, i * tau) * tau;
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] - F) / temp;
        }
        int grid_layer = i % 2;
        grid[grid_layer][h_steps] = fi2[i];
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[j + 1];
        }
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }

}

void indefinite_scheme_second_second(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    double *psi1, *psi2;
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau); //массивы из краевых условий
    }
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    //одинаковые значения для всей сетки
    double A_i = ((-1) * a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = (((2 * tau * a * a) / (h * h)) + 1) * (-1);

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
        mu_1 = (psi1[i] * h * (-1));
        mu_2 = (psi2[i] * h);
        betta[1] = mu_1;
        for (int j = 2; j < h_steps + 1; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            double F = grid[(i - 1) % 2][j - 1] + f((j - 1) * h, i * tau) * tau;
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] - F) / temp;
        }
        int grid_layer = i % 2;
        grid[grid_layer][h_steps] = ((mu_2 + hi_2 * betta[h_steps]) / (1 - hi_2 * alpha[h_steps]));
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[j + 1];
        }
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

void indefinite_scheme_first_second(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    double *fi1, *psi2;
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = u(0, i * tau);
        psi2[i] = psi_2(1, i * tau); //массивы из краевых условий
    }


    //одинаковые значения для всей сетки
    double A_i = ((-1) * a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = (((2 * tau * a * a) / (h * h)) + 1) * (-1);

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
        mu_2 = (psi2[i] * h);
        betta[1] = mu_1;
        for (int j = 2; j < h_steps + 1; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            double F = grid[(i - 1) % 2][j - 1] + f((j - 1) * h, i * tau) * tau;
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] - F) / temp;
        }
        int grid_layer = i % 2;
        grid[grid_layer][h_steps] = ((mu_2 + hi_2 * betta[h_steps]) / (1 - hi_2 * alpha[h_steps]));
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[j + 1];
        }
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

void indefinite_scheme_second_first(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //шаги сетки - величина каждого кусочка
    double h = 1.0 / h_steps;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    double *psi1, *fi2;
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        fi2[i] = u(1, i * tau); //массивы из краевых условий
    }

    //одинаковые значения для всей сетки
    double A_i = ((-1) * a * a * tau) / (h * h);
    double B_i = A_i;
    double C_i = (((2 * tau * a * a) / (h * h)) + 1) * (-1);

    //массивы из прогоночных коэффициентов
    double *alpha, *betta;
    alpha = (double *) calloc(h_steps + 1, sizeof(double));
    betta = (double *) calloc(h_steps + 1, sizeof(double));

    alpha[1] = 1;
    for (int i = 1; i < tau_steps + 1; i++) { //i = заполняем этот уровень, стоим на i-1
        betta[1] = (psi1[i] * h * (-1));
        for (int j = 2; j < h_steps + 1; j++) { //заполняем прогоночные коэффициенты
            double temp = (C_i - alpha[j - 1] * A_i);
            double F = grid[(i - 1) % 2][j - 1] + f((j - 1) * h, i * tau) * tau;
            alpha[j] = B_i / temp;
            betta[j] = (A_i * betta[j - 1] - F) / temp;
        }
        int grid_layer = i % 2;
        grid[grid_layer][h_steps] = fi2[i];
        for (int j = h_steps - 1; j > -1; j--) {
            grid[grid_layer][j] = alpha[j + 1] * grid[grid_layer][j + 1] + betta[j + 1];
        }
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

void definite_scheme_first_first(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps;
    double h = 1.0 / h_steps;
    double *fi1, *fi2;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = fi_1(0, i * tau);
        fi2[i] = fi_2(1, i * tau);
    }
    for (int i = 1; i < tau_steps + 1; i++) { //делаем i - ый слой по нижележащему
        get_next_step_def(grid[(i + 1) % 2], grid[i % 2], a, tau, h, (i) * tau, h_steps);
        grid[i % 2][0] = fi1[i];
        grid[i % 2][h_steps] = fi2[i];

        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

void definite_scheme_second_second(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *psi2;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau);
    }
    for (int i = 1; i < tau_steps + 1; i++) {
        get_next_step_def(grid[(i + 1) % 2], grid[i % 2], a, tau, h, (i) * tau, h_steps);
        grid[i % 2][0] = (psi1[i] * 2 * h + grid[i % 2][2] - 4 * grid[i % 2][1]) / (-3);
        grid[i % 2][h_steps] = (psi2[i] * 2 * h - grid[i % 2][h_steps - 2] + 4 * grid[i % 2][h_steps - 1]) / 3;
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }

}

void definite_scheme_first_second(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    double *fi1, *psi2;
    fi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    psi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        fi1[i] = fi_1(0, i * tau);
        psi2[i] = psi_2(1, i * tau);
    }
    for (int i = 1; i < tau_steps + 1; i++) {
        get_next_step_def(grid[(i - 1) % 2], grid[i % 2], a, tau, h, i * tau, h_steps);
        grid[i % 2][0] = fi1[i];
        grid[i % 2][h_steps] = (psi2[i] * 2 * h - grid[i % 2][h_steps - 2] + 4 * grid[i % 2][h_steps - 1]) / 3;
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

void definite_scheme_second_first(int tau_steps, int h_steps, double a, double **grid, int mode) {
    double tau = 1.0 / tau_steps; //atomic steps for every dimension
    double h = 1.0 / h_steps;
    double *psi1, *fi2;
    std::ofstream fout;
    if (mode == 1) {
        fout.open("results.dat");
    }
    psi1 = (double *) calloc(tau_steps + 1, sizeof(double));
    fi2 = (double *) calloc(tau_steps + 1, sizeof(double));
    for (int i = 0; i < tau_steps + 1; i++) {
        psi1[i] = psi_1(0, i * tau);
        fi2[i] = fi_2(1, i * tau);
    }
    for (int i = 1; i < tau_steps + 1; i++) {
        get_next_step_def(grid[(i - 1) % 2], grid[i % 2], a, tau, h, i * tau, h_steps);
        grid[i % 2][0] = (psi1[i] * 2 * h + grid[i % 2][2] - 4 * grid[i % 2][1]) / (-3);
        grid[i % 2][h_steps] = fi2[i];
        if((mode ==1)&&(i%15==0)){
            for(int k =0 ; k<h_steps+1 ;k++){
                fout << k*h << ' '<<grid[i%2][k] << std::endl;
            }
            fout<<std::endl<<std::endl;
        }
    }
    if (mode == 1) {
        fout.close();
    }
}

int main(int argc, char **argv) {

    int tau_steps = atoi(argv[1]);   //чему равно N, то есть всего точек, начиная с 0, N+1
    int h_steps = atoi(argv[2]);     //чему равно M, то есть всего точек, начиная с 0, M+1
    char *mode = argv[3];            // для функций 0 - тестовый режим, 1 - графический
    int scheme = argv[4][0] - '0';
    double a_cof = a();
    char mode_type;
    if (strcmp(mode, "test") == 0) {
        mode_type = 0;
    }
    if (strcmp(mode, "viz") == 0) {
        mode_type = 1;
    }
    double h = 1.0 / h_steps;
    double **grid = (double **) calloc(sizeof(double *), 2);
    grid[0] = (double *) calloc(h_steps + 1, sizeof(double));
    grid[1] = (double *) calloc(h_steps + 1, sizeof(double));
    for (int i = 0; i < h_steps + 1; i++) {
        grid[0][i] = u_0(i * h, 0);
    }

    switch (scheme) {
        case 1: {
            definite_scheme_first_first(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 2: {
            definite_scheme_second_second(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 3: {
            definite_scheme_first_second(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 4: {
            definite_scheme_second_first(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 5: {
            indefinite_scheme_first_first(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 6: {
            indefinite_scheme_second_second(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 7: {
            indefinite_scheme_first_second(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
        case 8: {
            indefinite_scheme_second_first(tau_steps, h_steps, a_cof, grid, mode_type);
            break;
        }
    }
    if (mode_type == 0) {
        double sum = 0, diff_max = -0.05;
        for (int i = 0; i < h_steps + 1; i++) {
            double diff = u(i * h, 1.0) - grid[tau_steps % 2][i];
            if (abs(diff) > diff_max) {
                diff_max = abs(diff);
            }
            sum += (diff * diff);
        }
        sum = sum * h;
        sum = sqrt(sum);
        double max_value = -0.05;
        for (int i = 0; i < h_steps + 1; i++) {
            double current_value = u(i * h, 1.0);
            if (abs(current_value) > max_value) {
                max_value = abs(current_value);
            }
        }
        double sum1 = 0;
        for (int i = 0; i < h_steps + 1; i++) {
            sum1 += grid[tau_steps % 2][i] * grid[tau_steps % 2][i];
        }
        sum1 = sum1 * h;
        sum1 = sqrt(sum1);
        printf("I2h abs: %f\n", sum);
        printf("Ch abs: %f\n", diff_max);
        printf("I2h rel: %f\n", sum / sum1);
        printf("Ch rel: %f\n", diff_max / max_value);
    }
}





