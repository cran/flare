#ifndef MYMATH_H
#define MYMATH_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

double sign(double x);

double max(double x,double y);

double max_abs_vec(double * x, int n);

double max_vec(double * x, int n);

void max_selc(double *x, double vmax, double *x_s, int n, int *n_s, double z);

double min(double x,double y);

double l1norm(double * x, int n);

void euc_proj(double * v, double z, int n);

double fun1(double lambda, double * v, double z, int n);

double mod_bisec(double * v, double z, int n);

void fabs_vc(double *v_in, double *v_out, int n);

void max_fabs_vc(double *v_in, double *v_out, double *vmax, int *n1, int n, double z);

void sort_up_bubble(double *v, int n);

void get_residual(double *r, double *y, double *A, double *x, int *xa_idx, int *nn, int *mm);

void get_dual(double *u, double *r, double *mmu, int *nn);

void get_dual1(double *u, double *r, double *mmu, int *nn);

void get_dual2(double *u, double *r, double *mmu, int *nn);

void get_grad(double *g, double *A, double *u, int *dd, int *nn);

void get_base(double *base, double *u, double *r, double *mmu, int *nn);

#endif
