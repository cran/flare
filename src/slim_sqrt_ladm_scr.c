#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "mymath.h"
#include <time.h>


void sqrt_ladm_scr(double *Y0, double *X0, double *X, double *XX0, double *XX, int *idx_scr, int num_scr, int ndata, int dim, double *alp, double *beta, double * mu, double *T, double rho, int *ite, double lambda, int max_ite, double prec, int intercept, int flag, int nlamb, double nrholamb)
{
    int j,k,m,w_idx,size_a,size_a1,size_a_pre;
    double gap_ext,max_dif,alp_dif,beta_dif,mu_dif,threshold,tmpd,alp_tild_sq,alp_th,ratio,epsT;
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
    //double *alp0 = (double*) malloc(ndata*sizeof(double));
    double *alp_pre = (double*) malloc(ndata*sizeof(double));
    double *alp_tild = (double*) malloc(ndata*sizeof(double));
    //double *mu0 = (double*) malloc(ndata*sizeof(double));
    double *mu_grad = (double*) malloc(ndata*sizeof(double));
    double *X_beta = (double*) malloc(ndata*sizeof(double));
    double *y_i = (double*) malloc(ndata*sizeof(double));
    double *Xy = (double*) malloc(num_scr*sizeof(double));
    double *g = (double*) malloc(num_scr*sizeof(double));
    double *r = (double*) malloc(ndata*sizeof(double));
    //double *X = (double*) malloc(num_scr*num_scr*sizeof(double));
    //double *Y = (double*) malloc(num_scr*sizeof(double));
    ratio = 0.5;
    epsT = 1e1;
    //stop = clock();
    //time1 += stop - start;
    
    //start = clock();
    size_a = 0;
    for (j=0; j<num_scr; j++) {
        //alp0[j] = alp[idx_scr[j]];
        beta0[j] = beta[idx_scr[j]];
        //mu0[j] = mu[idx_scr[j]];
        if (beta0[j]==0) {
            idx_i[j] = 1;
        }
        else{
            idx_a[size_a] = j;
            size_a++;
            idx_i[j] = 0;
        }
        //printf("j=%d,alp=%f,beta=%f,mu=%f \n",j,alp[j],beta0[j],mu[j]);
    }
    //printf("size_a=%d \n",size_a);
    
    //stop = clock();
    //time2 += stop - start;
    //start = clock();
    
    if(nlamb==0){
        for (j=0; j<num_scr; j++) {
            //if(j==1) printf("X row%d ",j);
            for (k=0; k<ndata; k++) {
                X[j*ndata+k] = X0[idx_scr[j]*ndata+k];
            //        if(j==1) printf(" %d,%f,%f ",k,X[j*ndata+k],X0[j*ndata+k]);
            }
            //if(j==1) printf("\n");
        }
    }
    else{
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                //if(j==0) printf("row%d ",j);
                for (k=0; k<ndata; k++) {
                    X[j*ndata+k] = X0[idx_scr[j]*ndata+k];
                    //    if(j==0) printf(" %f ",X[j*num_scr+k]);
                }
            }
        }
    }
    for (j=0; j<ndata; j++) {
        X_beta[j] = 0;
        //printf("X_beta ");
        for (k=0; k<num_scr; k++) {
            X_beta[j] += X[k*ndata+j]*beta0[k];
        }
        //printf(" %f ",X_beta[j]);
    }
    //printf(" \n");
    
    if(nlamb==0){
        if(flag==1){
            for (j=0; j<num_scr; j++) {
                //if(j==1) printf("XX row%d ",j);
                for (k=j; k<num_scr; k++) {
                    XX[j*dim+k] = XX0[j*dim+k];
                    //XX[j*dim+k] = 0;
                    //for (m=0; m<ndata; m++) {
                    //    XX[j*dim+k] += X[k*ndata+m] * X[j*ndata+m];
                    //    if(j==1 && k==1) printf(" %d,%f,%f ",k,X[k*ndata+m],X[j*ndata+m]);
                    //}
                    //if(j==1) printf(" %f ",XX[j*dim+k]);
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
    //printf("T=%f \n",*T);
    //stop = clock();
    //time3 += stop - start;
    //start = clock();
    //printf("T=%f \n",T);
    gap_ext = 1;
    *ite = 0;
    max_dif = 1;
    while((gap_ext !=0 || max_dif > prec) && *ite < max_ite){
        // update alpha
        alp_dif = 0;
        alp_tild_sq=0;
        for(j=0; j<ndata; j++){
            alp_pre[j] = alp[j];
            alp_tild[j] = Y0[j]-X_beta[j]-mu[j];
            alp_tild_sq += alp_tild[j]*alp_tild[j];
        }
        alp_th = 1-nrholamb/sqrt(alp_tild_sq);
        if(alp_th<0) alp_th=0;
        for(j=0; j<ndata; j++){
            alp[j] = alp_tild[j]*alp_th;
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
        //Q0 = dif_vec_l2norm(Y0,X_beta,ndata);
        //T0 = *T;
        //ite2=0;
        //if(*ite<1e4){
        //    gap_track = 1;
        //    while(gap_track == 1 && T0>epsT){
        //        threshold = 4/(rho*T0);
        //        lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, T0, threshold, intercept, num_scr);
        //        Q = Q0+inner_prod2(g,beta1,beta0,num_scr)+dif_vec_l2norm(beta1,beta0,num_scr)*T0/2;
        //        F = dif_l2norm(r,Y0,X,beta1,ndata,num_scr,size_a1,idx_a1);
        //        ite2++;
        //        if(F<Q) T0 = T0*ratio;
        //        else {
        //            if(ite2==1){T0 = *T; gap_track = 0;}
        //            else{T0 = T0/ratio; gap_track = 0;}
        //        }
        //    }
        //}
        //printf("Y=%f,y=%f,Xy=%f,g=%f \n",Y[0],y_i[0],Xy[0],g[0]);
        //threshold = 4/(rho*(T0));
        //lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, T0, threshold, intercept, num_scr);
        threshold = 1/(rho*(*T));
        lineaization_lasso_dantzig(XX, Xy, beta0, beta1, beta_tild, g, idx_a, idx_scr, &size_a, idx_a1, idx_i1, &size_a1, *T, threshold, intercept, num_scr);
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
                X_beta[j]+=X[w_idx*ndata+j]*beta0[w_idx];
            }
            mu_grad[j]=(alp[j]+X_beta[j]-Y0[j]);
            mu[j] += mu_grad[j];
            mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
        }
        max_dif = max_dif>mu_dif ? max_dif : mu_dif;
        //if(ite_ext == 0)
        //    printf("te_ext=%d,alp=%f,beta=%f,mu=%f \n",ite_ext,alp[0],beta0[0],mu[0]);
        //if(ite_ext % 1000 == 999)
        //printf("ite_ext=%d,alp_dif=%f,beta_dif=%f,mu_dif=%f \n",ite_ext,alp_dif,beta_dif,mu_dif);
        (*ite)++;
    }
    
    //for (j=0; j<dim; j++) {
    //    alp0[j] = 0;
    //    beta[j] = 0;
    //    mu0[j] = 0;
    //}
    
    for (j=0; j<num_scr; j++) {
        //alp[idx_scr[j]] = alp0[j];
        beta[idx_scr[j]] = beta0[j];
        //mu[idx_scr[j]] = mu0[j];
    }
    //stop = clock();
    //time0 += stop - start;
    //printf("ite=%d,alp=%f,beta=%f,mu=%f \n\n",*ite,alp[0],beta0[0],mu[0]);
    //printf("ite=%d,alp=%f,beta=%f,mu=%f \n\n",*ite,alp[1],beta0[1],mu[1]);
    //printf("ite=%d,alp=%f,beta=%f,mu=%f \n\n",*ite,alp[2],beta0[2],mu[2]);
    
    //start = clock();
    free(beta0);
    free(beta1);
    free(idx_a);
    free(idx_i);
    free(idx_a1);
    free(idx_i1);
    free(beta_tild);
    free(beta_pre);
    //free(alp0);
    free(alp_pre);
    free(alp_tild);
    //free(mu0);
    free(mu_grad);
    free(X_beta);
    free(y_i);
    free(Xy);
    free(g);
    free(r);
    //free(X);
    //free(Y);
    //stop = clock();
    //time1 += stop - start;
}

void slim_sqrt_ladm_scr(double *Y0, double *X0, double *XX, double *beta, int *n, int *d, double *rho0, int *ite_int, int *ite_int1, int *ite_int2, int *num_scr_1, int *num_scr_2, int *idx_scr, int * idx_scr_1, int * idx_scr_2, double * gamma, double * lambda, int * nnlambda, int *max_ite, double *prec, int * intercept)
{
    int j,k,m,ndata,dim,nlambda,ite1,ite2,ite,max_ite0,max_ite1,max_ite2,num_scr,num_scr1,num_scr2,num_scr1_tmp,num_scr2_tmp,flag,flag1,flag2;
    double T,T1,T2,rho,zero,eps,eps1,eps2,ilambda,sqrtn,nrholamb;
    
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
    sqrtn = sqrt((double) ndata);

    
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
        nrholamb = 1/(rho*sqrtn*ilambda);
        sqrt_ladm_scr(Y0,X0,X1,XX,XX1,idx_scr1,num_scr1,ndata,dim,alp,beta0,mu,&T1,rho,&ite1,ilambda,max_ite1,eps1,*intercept,flag1,m,nrholamb);
        //printf("m=%d,alp=%f,beta=%f,mu=%f,ite1=%d \n",m,alp[0],beta0[0],mu[0],ite1);
        sqrt_ladm_scr(Y0,X0,X2,XX,XX2,idx_scr2,num_scr2,ndata,dim,alp,beta0,mu,&T2,rho,&ite2,ilambda,max_ite2,eps2,*intercept,flag2,m,nrholamb);
        //printf("m=%d,alp=%f,beta=%f,mu=%f,ite2=%d \n",m,alp[0],beta0[0],mu[0],ite2);
        sqrt_ladm_scr(Y0,X0,X,XX,XX,idx_scr,num_scr,ndata,dim,alp,beta0,mu,&T,rho,&ite,ilambda,max_ite0,eps,*intercept,flag,m,nrholamb);
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
    //printf("C time0=%f,time1=%f,time2=%f,time3=%f \n",((double)(time0))/CLOCKS_PER_SEC,((double)(time1))/CLOCKS_PER_SEC,((double)(time2))/CLOCKS_PER_SEC,((double)(time3))/CLOCKS_PER_SEC);
    free(beta0);
    free(alp);
    free(mu);
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
}

