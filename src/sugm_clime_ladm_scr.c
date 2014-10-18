#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "mymath.h"
#include <time.h>

void sugm_clime_ladm_scr(double *X0, double *XX, double *beta, double *x0, int *d, double * gamma, double * lambda, int * nnlambda, double *rho0, int *col_cnz0, int * row_idx0, int *ite_int, int *ite_int1, int *ite_int2, int *num_scr_1, int *num_scr_2, int *idx_scr, int * idx_scr_1, int * idx_scr_2, int *max_ite, double *prec, int * idx0)
{
    int j,k,m,dim,nlambda,ite1,ite2,ite,max_ite0,max_ite1,max_ite2,num_scr,num_scr1,num_scr2,num_scr1_tmp,num_scr2_tmp,flag,flag1,flag2,cnz,idx,intercept;
    double T,T1,T2,rho,zero,eps,eps1,eps2,ilambda;
    
    dim = *d;
    rho = *rho0;
    idx = (*idx0)-1;
    //if(idx==0){
    //    time0 = 0;
    //    time1 = 0;
    //    time2 = 0;
    //    time3 = 0;
    //}
    T = *gamma;
    T1 = 0;
    T2 = 0;
    nlambda = *nnlambda;
    num_scr = dim;
    num_scr1 = *num_scr_1;
    num_scr2 = *num_scr_2;
    max_ite1 = *max_ite;
    max_ite2 = ceil(max_ite1/5);
    max_ite0 = ceil(max_ite1/10);
    eps1 = *prec;
    eps2 = *prec;
    eps = eps2*10;
    zero = 0;
    cnz = 0;
    intercept = 0;

    
    double *beta0 = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(dim*sizeof(double));
    double *mu = (double*) malloc(dim*sizeof(double));
    int *idx_scr1 = (int*) malloc(dim*sizeof(int));
    int *idx_scr2 = (int*) malloc(dim*sizeof(int));
    double *Y0 = (double*) malloc(dim*sizeof(double));
    double *X = (double*) malloc(dim*dim*sizeof(double));
    double *Y = (double*) malloc(dim*sizeof(double));
    double *X1 = (double*) malloc(dim*dim*sizeof(double));
    double *Y1 = (double*) malloc(dim*sizeof(double));
    double *X2 = (double*) malloc(dim*dim*sizeof(double));
    double *Y2 = (double*) malloc(dim*sizeof(double));
    double *XX1 = (double*) malloc(dim*dim*sizeof(double));
    double *XX2 = (double*) malloc(dim*dim*sizeof(double));
    
    for(j=0; j<dim; j++) {
        beta0[j] = 0;
        alp[j] = 0;
        mu[j] = 0;
        Y0[j] = 0;
    }
    Y0[idx] = 1;
    for(j=0; j<num_scr1; j++) {
        idx_scr1[j] = idx_scr_1[j]-1;
    //    printf("%d=%d,%d ",j,idx_scr1[j],idx_scr_1[j]);
    }
    //printf("num_scr1=%d,j=%d \n \n",num_scr1,j);
    for(j=0; j<num_scr2; j++) {
        idx_scr2[j] = idx_scr_2[j]-1;
    //    printf("%d=%d,%d ",j,idx_scr2[j],idx_scr_2[j]);
    }
    //printf("num_scr2=%d,j=%d \n \n",num_scr2,j);
    for(j=0; j<num_scr; j++) {
        idx_scr[j] = idx_scr[j]-1;
    //    printf("%d=%d ",j,idx_scr[j]);
    }
    //printf("num_scr=%d,j=%d \n \n",num_scr,j);
    flag = 0;
    flag1 = 1;
    flag2 = 1;
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m];//*ndata
        
        dantzig_ladm_scr(Y0,X0,Y1,X1,XX1,idx_scr1,num_scr1,dim,alp,beta0,mu,&T1,rho,&ite1,ilambda,max_ite1,eps1,intercept,flag1,m);
        //printf("m=%d,alp=%f,beta=%f,mu=%f,ite1=%d \n",m,alp[0],beta0[0],mu[0],ite1);
        dantzig_ladm_scr(Y0,X0,Y2,X2,XX2,idx_scr2,num_scr2,dim,alp,beta0,mu,&T2,rho,&ite2,ilambda,max_ite2,eps2,intercept,flag2,m);
        //printf("m=%d,alp=%f,beta=%f,mu=%f,ite2=%d \n",m,alp[0],beta0[0],mu[0],ite2);
        dantzig_ladm_scr(Y0,X0,Y,X,XX,idx_scr,num_scr,dim,alp,beta0,mu,&T,rho,&ite,ilambda,max_ite0,eps,intercept,flag,m);
        //printf("m=%d,alp=%f,beta=%f,mu=%f,ite=%d \n",m,alp[0],beta0[0],mu[0],ite);
        
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
                if(idx != k) {
                    x0[cnz] = beta0[k];
                    row_idx0[cnz] = m*dim+k;
                    cnz++;
                }
            }
        }
        if(num_scr1>num_scr1_tmp){
        //    printf("before, flag1=%d,num_scr1=%d \n",flag1,num_scr1);
        //    free(XX1);
        //    double *XX1 = (double*) malloc(num_scr1*num_scr1*sizeof(double));
            flag1=1;
        //    printf("after, flag1=%d,num_scr1=%d \n",flag1,num_scr1);
        }
        else{
        //    printf("before, flag1=%d,num_scr1=%d \n",flag1,num_scr1);
        //    free(XX1);
        //    double *XX1 = (double*) malloc(num_scr1*num_scr1*sizeof(double));
            flag1 = 0;
        //    printf("after, flag1=%d,num_scr1=%d \n",flag1,num_scr1);
        }
        if(num_scr2>num_scr2_tmp){
        //    printf("before, flag2=%d,num_scr2=%d  \n",flag2,num_scr2);
        //    free(XX2);
        //    double *XX2 = (double*) malloc(num_scr2*num_scr2*sizeof(double));
            flag2=1;
        //    printf("after, flag2=%d,num_scr2=%d  \n",flag2,num_scr2);
        }
        else{
        //    printf("before, flag2=%d,num_scr2=%d  \n",flag2,num_scr2);
        //    free(XX2);
        //    double *XX2 = (double*) malloc(num_scr2*num_scr2*sizeof(double));
            flag2 = 0;
        //    printf("after, flag2=%d,num_scr2=%d  \n",flag2,num_scr2);
        }
        //flag = 0;
        
        ite_int[m] = ite;
        ite_int1[m] = ite1;
        ite_int2[m] = ite2;
    }
    *col_cnz0 = cnz;
    //if(idx==dim-1)
        //printf("C time0=%f,time1=%f,time2=%f,time3=%f \n",((double)(time0))/CLOCKS_PER_SEC,((double)(time1))/CLOCKS_PER_SEC,((double)(time2))/CLOCKS_PER_SEC,((double)(time3))/CLOCKS_PER_SEC);
    free(beta0);
    free(alp);
    free(mu);
    free(idx_scr1);
    free(idx_scr2);
    //free(XX);
    free(Y0);
    free(X);
    free(Y);
    free(X1);
    free(Y1);
    free(X2);
    free(Y2);
    free(XX1);
    free(XX2);
}
