#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "mymath.h"
#include <time.h>


void lad_ladm_scr_btr(double *Y0, double *X0, double *X, double *XX0, double *XX, int *idx_scr, int num_scr, int ndata, int dim, double *alp, double *beta, double * mu, double *T, double rho, int *ite, double lambda, int max_ite, double prec, int intercept, int flag, int nlamb, double nrholamb)
{
    int j,k,m,w_idx,size_a,size_a1,size_a_pre,gap_track,ite2;
    double gap_ext,max_dif,alp_dif,beta_dif,mu_dif,threshold,tmpd,alp_abs,alp_tmp,T0,ratio,epsT,Q,Q0,F;
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
    double *alp_pre = (double*) malloc(ndata*sizeof(double));
    double *alp_tild = (double*) malloc(ndata*sizeof(double));
    double *mu_grad = (double*) malloc(ndata*sizeof(double));
    double *X_beta = (double*) malloc(ndata*sizeof(double));
    double *y_i = (double*) malloc(ndata*sizeof(double));
    double *Xy = (double*) malloc(num_scr*sizeof(double));
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
    }
    else{
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                for (k=0; k<ndata; k++) {
                    X[j*ndata+k] = X0[idx_scr[j]*ndata+k];
                }
            }
        }
    }
    for (j=0; j<ndata; j++) {
        X_beta[j] = 0;
        for (k=0; k<num_scr; k++) {
            X_beta[j] += X[k*ndata+j]*beta0[k];
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
        }
    }
    //stop = clock();
    //time3 += stop - start;
    //start = clock();
    gap_ext = 1;
    *ite = 0;
    max_dif = 1;
    while((gap_ext !=0 || max_dif > prec) && *ite < max_ite){
        // update alpha
        alp_dif = 0;
        for(j=0; j<ndata; j++){
            alp_pre[j] = alp[j];
            alp_tild[j] = Y0[j]-X_beta[j]-mu[j];
            alp_abs = fabs(alp_tild[j]);
            alp_tmp = alp_abs>nrholamb ? (alp_abs-nrholamb) : 0;
            alp[j] = (alp_tild[j]>0 ? 1 : -1)*alp_tmp;
            alp_dif += fabs(alp_pre[j] - alp[j]);
        }
        max_dif = alp_dif;
        
        // update omega
        for(j=0; j<ndata; j++){
            y_i[j] = Y0[j]-alp[j]-mu[j];
        }
        for(j=0; j<num_scr; j++){
            Xy[j] = 0;
            for(k=0; k<ndata; k++){
                Xy[j] += X[j*ndata+k]*y_i[k];
            }
        }
        size_a_pre = size_a;
        for(j=0; j<num_scr; j++){
            g[j] = 0;
            for(k=0; k<size_a; k++){
                w_idx = idx_a[k];
                g[j] += XX[w_idx*dim+j]*beta0[w_idx];
            }
            g[j] = g[j] - Xy[j];
        }
        
        T0 = *T;
        Q0 = dif_l2norm(r,y_i,X,beta0,ndata,dim,size_a,idx_a)/2;
        ite2=0;
        if((*ite)<1e4){
            gap_track = 1;
            while(gap_track == 1 && T0>epsT){
                threshold = 1/((rho)*T0);
                lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, T0, threshold, intercept, num_scr);
                Q = Q0+inner_prod2(g,beta1,beta0,num_scr)+dif_vec_l2norm(beta1,beta0,num_scr)*T0/2;
                F = dif_l2norm(r,y_i,X,beta1,ndata,num_scr,size_a1,idx_a1)/2;
                ite2++;
                if(F<Q) T0 = T0*ratio;
                else {
                    if(ite2==1){T0 = *T; gap_track = 0;}
                    else{T0 = T0/ratio; gap_track = 0;}
                }
            }
        }
        
        threshold = 1/(rho*(T0));
        lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, T0, threshold, intercept, num_scr);
        beta_dif = 0;
        for(j=0; j<num_scr; j++){
            beta_dif += fabs(beta0[j] - beta1[j]);
            beta0[j] = beta1[j];
        }
        max_dif = max_dif>beta_dif ? max_dif : beta_dif;
        size_a = size_a1;
        idx_tmp = idx_a; idx_a = idx_a1; idx_a1 = idx_tmp;
        idx_tmp = idx_i; idx_i = idx_i1; idx_i1 = idx_tmp;
        gap_ext = size_a - size_a_pre;
        
        // update mu
        mu_dif = 0;
        for(j=0; j<ndata; j++){
            X_beta[j]=0;
            for(k=0; k<size_a; k++){
                w_idx = idx_a[k];
                X_beta[j]+=X[w_idx*ndata+j]*beta1[w_idx];
            }
            mu_grad[j]=(alp[j]+X_beta[j]-Y0[j]);
            mu[j] += mu_grad[j];
            mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
        }
        max_dif = max_dif>mu_dif ? max_dif : mu_dif;
        (*ite)++;
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
    free(alp_pre);
    free(alp_tild);
    free(mu_grad);
    free(X_beta);
    free(y_i);
    free(Xy);
    free(g);
    free(r);
    //stop = clock();
    //time1 += stop - start;
}

void slim_lad_ladm_scr_btr(double *Y0, double *X0, double *XX, double *beta, int *n, int *d, double *rho0, int *ite_int, int *ite_int1, int *ite_int2, int *num_scr_1, int *num_scr_2, int *idx_scr, int * idx_scr_1, int * idx_scr_2, double * gamma, double * lambda, int * nnlambda, int *max_ite, double *prec, int * intercept)
{
    int j,k,m,ndata,dim,nlambda,ite1,ite2,ite,max_ite0,max_ite1,max_ite2,num_scr,num_scr1,num_scr2,num_scr1_tmp,num_scr2_tmp,flag,flag1,flag2;
    double T,T1,T2,rho,zero,eps,eps1,eps2,ilambda,nrholamb;
    
    //time0 = 0;
    //time1 = 0;
    //time2 = 0;
    //time3 = 0;
    dim = *d;
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

    
    double *beta0 = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(ndata*sizeof(double));
    double *mu = (double*) malloc(ndata*sizeof(double));
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
    
    for(j=0; j<ndata; j++) {
        alp[j] = 0;
        mu[j] = 0;
    }
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
    flag = 0;
    flag1 = 1;
    flag2 = 1;
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m];//*ndata
        nrholamb = 1/(rho*ndata*ilambda);
        lad_ladm_scr_btr(Y0,X0,X1,XX,XX1,idx_scr1,num_scr1,ndata,dim,alp,beta0,mu,&T1,rho,&ite1,ilambda,max_ite1,eps1,*intercept,flag1,m,nrholamb);
        lad_ladm_scr_btr(Y0,X0,X2,XX,XX2,idx_scr2,num_scr2,ndata,dim,alp,beta0,mu,&T2,rho,&ite2,ilambda,max_ite2,eps2,*intercept,flag2,m,nrholamb);
        lad_ladm_scr_btr(Y0,X0,X,XX,XX,idx_scr,num_scr,ndata,dim,alp,beta0,mu,&T,rho,&ite,ilambda,max_ite0,eps,*intercept,flag,m,nrholamb);
        
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
    free(beta0);
    free(alp);
    free(mu);
    free(idx_scr1);
    free(idx_scr2);
    free(X);
    free(Y);
    free(X1);
    free(Y1);
    free(X2);
    free(Y2);
    free(XX1);
    free(XX2);
}

