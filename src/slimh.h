#ifndef SLIMH_H
#define SLIMH_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define delta 1e-8
#define innerIter 1000
#define outerIter 1000

void eplb(double * x, double *root, int * steps, double * v,int n, double z, double lambda0);
void slim1(double *x, double *v, int n, double rho);
void slim2(double *x, double *v, int n, double rho);
void slimInf(double *x, double * c, int * iter_step, double *v,  int n, double rho, double c0);
void zerofind(double *root, int * iterStep, double v, double p, double c, double x0);
double norm(double * v, double p, int n);
void slimO(double *x, double * cc, int * iter_step, double *v,  int n, double rho, double p);
void slim(double *x, double * c, int * iter_step, double * v, int n, double rho, double p, double c0);

#endif
