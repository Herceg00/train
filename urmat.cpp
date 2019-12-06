#include <iostream>
#include "math.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"
#include "string.h"



double u(double x, double t){
    return 0;
    //tutor function
    
}

double f(double x, double t){
    return 0;
    //tutor function

}

//prev_tau is used to calculate f(X_i,prev_tau)

double get_next_step_def(double *prev_layer,double *cur_layer,double a,double tau,double h,double prev_tau,int h_steps){
    for(int i=1;i<h_steps-1;i++){

        cur_layer[i] = (tau*a*a*(prev_layer[i-1] - 2*prev_layer[i]+prev_layer[i+1]))/(h*h) + f(i*h,prev_tau) + prev_layer[i];
    }

}



int main(int argc,char** argv) {
    int tau_steps = atoi(argv[1]);
    //doing only test mode now
}

//h_steps = grid dimension on OX coordinate
//tau_steps = grid gimaension on OY coordinate

void definite_scheme_first_first(int tau_steps,int h_steps){

    //DOING INITIAIZATION

    double tau = 1.0/tau_steps; //atomic steps for every dimension
    double h = 1.0/h_steps;
    double *fi1,*fi2;
    fi1 = (double*)calloc(tau_steps, sizeof(double));
    fi2 = (double*)calloc(tau_steps,sizeof(double));
    for(int i = 0;i<tau_steps;i++){
        fi1[i] = u(0,i*tau);
        fi2[i] = u(1,i*tau);
        
    }
    double **grid = (double**) calloc(sizeof(double *),2);
    grid[0] = (double*) calloc(h_steps,sizeof(double));
    grid[1] = (double*) calloc(h_steps, sizeof(double));
    for(int i = 0;i<h_steps;i++){
        grid[0][i] = u(i*h,0);
    }



    //INITIALIZATION DONE





}

