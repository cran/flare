#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "mymath.h"
#include <time.h>


void lasso_ladm_scr(double *XY0, double *Xy, double *X0, double *X, double *XX0, double *XX, int *idx_scr, int num_scr, int ndata, int dim, double *beta, double *T, double rho, int *ite, double lambda, int max_ite, double prec, int intercept, int flag, int nlamb)
{
    int j,k,m,w_idx,size_a,size_a1,size_a_pre;
    double gap_ext,max_dif,beta_dif,threshold,tmpd,ratio,epsT;
    int * idx_tmp;
    
    
    //start = clock();
    double *beta0 = (double*) malloc(num_scr*sizeof(double));
    double *beta1 = (double*) malloc(num_scr*sizeof(double));
    int *idx_a = (int*) malloc(num_scr*sizeof(int)); //active set
    int *idx_i = (int*) malloc(num_scr*sizeof(int)); //inactive set
    int *idx_a1 = (int*) malloc(num_scr*sizeof(int)); //active set
    int *idx_i1 = (int*) malloc(num_scr*sizeof(int)); //inactive set
    double *beta_tild = (double*) malloc(num_scr*sizeof(double));
    double *beta_pre = (double*) malloc(num_scr*sizeof(double));
    double *X_beta = (double*) malloc(ndata*sizeof(double));
    double *y_i = (double*) malloc(ndata*sizeof(double));
    double *g = (double*) malloc(num_scr*sizeof(double));
    double *r = (double*) malloc(ndata*sizeof(double));
    ratio = 0.5;
    epsT = 1e1;
    //stop = clock();
    //time1 += stop - start;
    
    //start = clock();
    size_a = 0;
    for (j=0; j<num_scr; j++) {
        beta0[j] = beta[idx_scr[j]];
        if (beta0[j]==0) {
            idx_i[j] = 1;
        }
        else{
            idx_a[size_a] = j;
            size_a++;
            idx_i[j] = 0;
        }
    }
    
    //stop = clock();
    //time2 += stop - start;
    //start = clock();
    
    if(nlamb==0){
        for (j=0; j<num_scr; j++) {
            for (k=0; k<ndata; k++) {
                X[j*ndata+k] = X0[idx_scr[j]*ndata+k];
            }
        }
        for(j=0; j<num_scr; j++){
            Xy[j] = XY0[idx_scr[j]];
        }
    }
    else{
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                for (k=0; k<ndata; k++) {
                    X[j*ndata+k] = X0[idx_scr[j]*ndata+k];
                }
            }
            for(j=0; j<num_scr; j++){
                Xy[j] = XY0[idx_scr[j]];
            }
        }
    }
    
    if(nlamb==0){
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                for (k=j; k<num_scr; k++) {
                    XX[j*dim+k] = XX0[j*dim+k];
                    XX[k*dim+j] = XX[j*dim+k];
                }
            }
            *T = 0;
            for (j=0; j<num_scr; j++) {
                tmpd = 0;
                for (k=0; k<num_scr; k++) {
                    tmpd += fabs(XX[j*dim+k]);
                }
                *T = max(*T,tmpd);
            }
            *T = (*T)*0.8;
            for(j=0; j<num_scr; j++){
                Xy[j] = XY0[idx_scr[j]];
            }
        }
    }
    else{
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                for (k=j; k<num_scr; k++) {
                    XX[j*dim+k] = 0;
                    for (m=0; m<ndata; m++) {
                        XX[j*dim+k] += X[k*ndata+m] * X[j*ndata+m];
                    }
                    XX[k*dim+j] = XX[j*dim+k];
                }
            }
            *T = 0;
            for (j=0; j<num_scr; j++) {
                tmpd = 0;
                for (k=0; k<num_scr; k++) {
                    tmpd += fabs(XX[j*dim+k]);
                }
                *T = max(*T,tmpd);
            }
            *T = (*T)*0.8;
            for(j=0; j<num_scr; j++){
                Xy[j] = XY0[idx_scr[j]];
            }
        }
    }
    //stop = clock();
    //time3 += stop - start;
    //start = clock();
    gap_ext = 1;
    *ite = 0;
    max_dif = 1;
    threshold = lambda/(*T);
    while((gap_ext !=0 || max_dif > prec) && *ite < max_ite){
        
        // update omega
        size_a_pre = size_a;
        for(j=0; j<num_scr; j++){
            g[j] = 0;
            for(k=0; k<size_a; k++){
                w_idx = idx_a[k];
                g[j] += XX[w_idx*dim+j]*beta0[w_idx];
            }
            g[j] = g[j] - Xy[j];
        }
        lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, *T, threshold, intercept, num_scr);
        beta_dif = 0;
        for(j=0; j<num_scr; j++){
            beta_dif += fabs(beta0[j] - beta1[j]);
            beta0[j] = beta1[j];
        }
        max_dif = beta_dif;
        size_a = size_a1;
        idx_tmp = idx_a; idx_a = idx_a1; idx_a1 = idx_tmp;
        idx_tmp = idx_i; idx_i = idx_i1; idx_i1 = idx_tmp;
        gap_ext = size_a - size_a_pre;
        (*ite)++;
        //if(*ite % 1000 == 999)
            //printf("ite_ext=%d,beta_dif=%f,max_dif=%f \n",*ite,beta_dif,max_dif);
    }
    
    for (j=0; j<num_scr; j++) {
        beta[idx_scr[j]] = beta0[j];
    }
    //stop = clock();
    //time0 += stop - start;
    
    //start = clock();
    free(beta0);
    free(beta1);
    free(idx_a);
    free(idx_i);
    free(idx_a1);
    free(idx_i1);
    free(beta_tild);
    free(beta_pre);
    free(X_beta);
    free(y_i);
    free(g);
    free(r);
    //stop = clock();
    //time1 += stop - start;
}

void slim_lasso_ladm_scr(double *XY0, double *X0, double *XX, double *beta, int *n, int *d, int *ite_int, int *ite_int1, int *ite_int2, int *num_scr_1, int *num_scr_2, int *idx_scr, int * idx_scr_1, int * idx_scr_2, double * gamma, double * lambda, int * nnlambda, int *max_ite, double *prec, int * intercept)
{
    int j,k,m,ndata,dim,nlambda,ite1,ite2,ite,max_ite0,max_ite1,max_ite2,num_scr,num_scr1,num_scr2,num_scr1_tmp,num_scr2_tmp,flag,flag1,flag2;
    double T,T1,T2,rho,zero,eps,eps1,eps2,ilambda;
    
    //time0 = 0;
    //time1 = 0;
    //time2 = 0;
    //time3 = 0;
    dim = *d;
    ndata = *n;
    T = (*gamma)*0.8;
    T1 = 0;
    T2 = 0;
    nlambda = *nnlambda;
    num_scr = dim;
    num_scr1 = *num_scr_1;
    num_scr2 = *num_scr_2;
    max_ite1 = *max_ite;
    max_ite2 = *max_ite;
    max_ite0 = ceil(max_ite2/10);
    eps1 = *prec;
    eps2 = *prec;
    eps = eps2*10;
    zero = 0;

    
    double *beta0 = (double*) malloc(dim*sizeof(double));
    int *idx_scr1 = (int*) malloc(dim*sizeof(int));
    int *idx_scr2 = (int*) malloc(dim*sizeof(int));
    double *X = (double*) malloc(ndata*dim*sizeof(double));
    double *Y = (double*) malloc(ndata*sizeof(double));
    double *X1 = (double*) malloc(ndata*dim*sizeof(double));
    double *XY1 = (double*) malloc(dim*sizeof(double));
    double *X2 = (double*) malloc(ndata*dim*sizeof(double));
    double *XY2 = (double*) malloc(dim*sizeof(double));
    double *XY = (double*) malloc(dim*sizeof(double));
    double *XX1 = (double*) malloc(dim*dim*sizeof(double));
    double *XX2 = (double*) malloc(dim*dim*sizeof(double));
    
    for(j=0; j<dim; j++) {
        beta0[j] = 0;
    }
    for(j=0; j<num_scr1; j++) {
        idx_scr1[j] = idx_scr_1[j]-1;
    }
    for(j=0; j<num_scr2; j++) {
        idx_scr2[j] = idx_scr_2[j]-1;
    }
    for(j=0; j<num_scr; j++) {
        idx_scr[j] = idx_scr[j]-1;
    }
    rho = 1;
    flag = 0;
    flag1 = 1;
    flag2 = 1;
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m]*ndata;
        lasso_ladm_scr(XY0,XY1,X0,X1,XX,XX1,idx_scr1,num_scr1,ndata,dim,beta0,&T1,rho,&ite1,ilambda,max_ite1,eps1,*intercept,flag1,m);
        lasso_ladm_scr(XY0,XY2,X0,X2,XX,XX2,idx_scr2,num_scr2,ndata,dim,beta0,&T2,rho,&ite2,ilambda,max_ite2,eps2,*intercept,flag2,m);
        lasso_ladm_scr(XY0,XY,X0,X,XX,XX,idx_scr,num_scr,ndata,dim,beta0,&T,rho,&ite,ilambda,max_ite0,eps,*intercept,flag,m);
        
        num_scr1_tmp = num_scr1;
        num_scr2_tmp = num_scr2;
        for(k=0;k<dim;k++){
            beta[m*dim+k] = beta0[k];
            if(beta0[k]!=0){
                for(j=0; j<num_scr1; j++) {
                    if(idx_scr1[j] == k) break;
                }
                if(j==num_scr1){
                    idx_scr1[num_scr1] = k;
                    num_scr1++;
                }
                
                for(j=0; j<num_scr2; j++) {
                    if(idx_scr2[j] == k) break;
                }
                if(j==num_scr2){
                    idx_scr2[num_scr2] = k;
                    num_scr2++;
                }
            }
        }
        if(num_scr1>num_scr1_tmp){
            flag1=1;
        }
        else{
            flag1 = 0;
        }
        if(num_scr2>num_scr2_tmp){
            flag2=1;
        }
        else{
            flag2 = 0;
        }
        //flag = 0;
        
        ite_int[m] = ite;
        ite_int1[m] = ite1;
        ite_int2[m] = ite2;
    }
    //printf("C time0=%f,time1=%f,time2=%f,time3=%f \n",((double)(time0))/CLOCKS_PER_SEC,((double)(time1))/CLOCKS_PER_SEC,((double)(time2))/CLOCKS_PER_SEC,((double)(time3))/CLOCKS_PER_SEC);
    free(beta0);
    free(idx_scr1);
    free(idx_scr2);
    free(X);
    free(Y);
    free(X1);
    free(XY1);
    free(X2);
    free(XY2);
    free(XY);
    free(XX1);
    free(XX2);
}

