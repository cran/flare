#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "mymath.h"
#include <time.h>


void tiger_lasso_ladm_scr(double *Y0, double *XY0, double *Xy, double *X0, double *X, double *XX0, double *XX, int *idx_scr, int num_scr, int ndata, int dim, double *beta, double *T, double rho, int *ite, double lambda, int max_ite, double prec, int flag, int nlamb, double *tau)
{
    int j,k,m,w_idx,size_a,size_a1,size_a_pre;
    double gap_ext,max_dif,beta_dif,threshold,tmpd,ratio,epsT,tau0,tau1,tmp;
    int * idx_tmp;
    //char c;
    
    tau0 = *tau;
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
    while((gap_ext !=0 || max_dif > prec) && *ite < max_ite){
        
        threshold = lambda*tau0/(*T);
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
        lineaization_clime(beta0, beta1, beta_tild, g, idx_a, &size_a, idx_a1, idx_i1, &size_a1, *T, threshold, num_scr);
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
        
        tau1 = 0;
        for (j=0; j<ndata; j++) {
            X_beta[j]=0;
            for(k=0; k<size_a; k++){
                w_idx = idx_a[k];
                X_beta[j]+=X[w_idx*ndata+j]*beta0[w_idx];
            }
            tmp = Y0[j] - X_beta[j];
            tau1 += tmp*tmp;
            //printf(" %f ",X_beta[j]);
        }
        //printf(" \n");
        tau1=sqrt(tau1/ndata);
        max_dif = max(fabs(tau1-tau0)/tau1,beta_dif);
        tau0 = tau1;
        
        (*ite)++;
        //if(*ite % 1000 == 999)
        //printf("ite_ext=%d,beta_dif=%f,max_dif=%f \n",*ite,beta_dif,max_dif);
    }
    
    *tau = tau0;
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


void sugm_tiger_ladm_scr(double *Y0, double *X0, double *XY0, double *XX, double *beta, double *x0, int *d, int *n, double * gamma, double * lambda, int * nnlambda, double *rho0, int *col_cnz0, int * row_idx0, int *ite_int, int *ite_int1, int *ite_int2, int *ite_int3, int *ite_int4, int *ite_int5, int *num_scr_1, int *num_scr_2, int *idx_scr, int * idx_scr_1, int * idx_scr_2, int *max_ite, double *prec, int * idx0)

{
    int j,k,m,ndata,dim,dim0,nlambda,ite1,ite2,ite,max_ite0,max_ite1,max_ite2,num_scr,num_scr1,num_scr2,num_scr1_tmp,num_scr2_tmp,flag,flag1,flag2,intercept;
    int idx,y_col,cnz,mdim0;
    double T,T1,T2,rho,zero,eps,eps1,eps2,ilambda,sqrtn,nrholamb,tau0;
    
    //time0 = 0;
    //time1 = 0;
    //time2 = 0;
    //time3 = 0;
    dim0 = *d;
    dim = dim0-1;
    ndata = *n;
    rho = *rho0;
    T = *gamma;
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
    sqrtn = sqrt((double) ndata);
    intercept = 0;
    idx = (*idx0)-1;
    cnz = 0;
    
    
    double *beta0 = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(ndata*sizeof(double));
    double *mu = (double*) malloc(ndata*sizeof(double));
    int *idx_a = (int*) malloc(dim*sizeof(int)); //active set
    int *idx_i = (int*) malloc(dim*sizeof(int)); //inactive set
    int *idx_scr1 = (int*) malloc(dim*sizeof(int));
    int *idx_scr2 = (int*) malloc(dim*sizeof(int));
    double *X = (double*) malloc(ndata*dim*sizeof(double));
    double *Y = (double*) malloc(ndata*sizeof(double));
    double *X1 = (double*) malloc(ndata*dim*sizeof(double));
    double *Y1 = (double*) malloc(ndata*sizeof(double));
    double *X2 = (double*) malloc(ndata*dim*sizeof(double));
    double *Y2 = (double*) malloc(ndata*sizeof(double));
    double *XX1 = (double*) malloc(dim*dim*sizeof(double));
    double *XX2 = (double*) malloc(dim*dim*sizeof(double));
    double *XY = (double*) malloc(dim*sizeof(double));
    double *XY1 = (double*) malloc(dim*sizeof(double));
    double *XY2 = (double*) malloc(dim*sizeof(double));
    int *idx_col = (int*) malloc(dim*sizeof(int)); //column set
    
    for(j=0; j<ndata; j++) {
        alp[j] = 0;
        mu[j] = 0;
    }
    for(j=0; j<dim; j++) {
        beta0[j] = 0;
    }
    for(j=0; j<num_scr1; j++) {
        idx_scr1[j] = idx_scr_1[j]-1;
    //        printf("%d=%d,%d ",j,idx_scr1[j],idx_scr_1[j]);
    }
    //printf("num_scr1=%d,j=%d \n \n",num_scr1,j);
    for(j=0; j<num_scr2; j++) {
        idx_scr2[j] = idx_scr_2[j]-1;
    //        printf("%d=%d,%d ",j,idx_scr2[j],idx_scr_2[j]);
    }
    //printf("num_scr2=%d,j=%d \n \n",num_scr2,j);
    for(j=0; j<num_scr; j++) {
        idx_scr[j] = idx_scr[j]-1;
    //        printf("%d=%d ",j,idx_scr[j]);
    }
    //printf("num_scr=%d,j=%d \n \n",num_scr,j);
    
    y_col = 0;
    for(j=0; j<dim; j++) {
        if(y_col == idx) y_col++;
        idx_col[j] = y_col;
        y_col++;
    }
    flag = 0;
    flag1 = 1;
    flag2 = 1;
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m];//*ndata
        nrholamb = 1/(rho*sqrtn*ilambda);
        sqrt_ladm_scr(Y0,X0,X1,XX,XX1,idx_scr1,num_scr1,ndata,dim,alp,beta0,mu,&T1,rho,&ite1,ilambda,max_ite1,eps1,intercept,flag1,m,nrholamb);
        //printf("idx=%d,m=%d,alp=%f,beta=%f,mu=%f,ite1=%d \n",idx,m,alp[0],beta0[0],mu[0],ite1);
        sqrt_ladm_scr(Y0,X0,X2,XX,XX2,idx_scr2,num_scr2,ndata,dim,alp,beta0,mu,&T2,rho,&ite2,ilambda,max_ite2,eps2,intercept,flag2,m,nrholamb);
        //printf("idx=%d,m=%d,alp=%f,beta=%f,mu=%f,ite2=%d \n",idx,m,alp[0],beta0[0],mu[0],ite2);
        sqrt_ladm_scr(Y0,X0,X,XX,XX,idx_scr,num_scr,ndata,dim,alp,beta0,mu,&T,rho,&ite,ilambda,max_ite0,eps,intercept,flag,m,nrholamb);
        //printf("idx=%d,m=%d,alp=%f,beta=%f,mu=%f,ite=%d \n",idx,m,alp[0],beta0[0],mu[0],ite);
        
        ite_int[m] = ite;
        ite_int1[m] = ite1;
        ite_int2[m] = ite2;
        
        tau0=0;
        for(j=0;j<ndata;j++)
            tau0 += alp[j]*alp[j];
        tau0 = sqrt(tau0/ndata);
        
        
        num_scr1_tmp = num_scr1;
        num_scr2_tmp = num_scr2;
        mdim0 = m*dim0;
        beta[mdim0+idx] = 1/(tau0*tau0);
        for(k=0;k<dim;k++){
            beta[mdim0+idx_col[k]] = -beta0[k]*beta[mdim0+idx];
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
                x0[cnz] = beta[mdim0+idx_col[k]];
                row_idx0[cnz] = mdim0+idx_col[k];
                cnz++;
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
        //if(m==10) return;
        //flag = 0;
    }
    *col_cnz0 = cnz;
    //printf("C time0=%f,time1=%f,time2=%f,time3=%f \n",((double)(time0))/CLOCKS_PER_SEC,((double)(time1))/CLOCKS_PER_SEC,((double)(time2))/CLOCKS_PER_SEC,((double)(time3))/CLOCKS_PER_SEC);
    free(beta0);
    free(alp);
    free(mu);
    free(idx_a);
    free(idx_i);
    free(idx_scr1);
    free(idx_scr2);
    //free(XX);
    free(X);
    free(Y);
    free(X1);
    free(Y1);
    free(X2);
    free(Y2);
    free(XX1);
    free(XX2);
    free(XY);
    free(XY1);
    free(XY2);
    free(idx_col);
}

