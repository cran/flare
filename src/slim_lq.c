#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "slimh.h"
#include "mymath.h"

void slim_lq(double *Y, double *X, double *XX, double *beta, int *n, int *d, double *rho, int *ite_cnt_ext, int *ite_cnt_int1, int *ite_cnt_int2, double *q, double *lambda, int * nnlambda, int *max_ite, double *prec, int * intercept, double * obj, double * runt)
{
    int j,k,m,ndata,dim,ileft,iright,junk_a,size_a,size_a_pre,w_idx,rs_idx,nlambda;
    int ite1,ite2,ite_ext,gap_ext,max_ite1,max_ite2,max_ite3;
    double gamma,cc,qrt_n, beta_sum1,beta_dif_sum1,beta_sum2,gap_int,F,beta_norm1;
    double ilambda,ilambda_tmp,tmp1,tmp2,beta_dif,eps,eps1,eps2,max_dif,beta_out1,mu_dif,rho0;
    clock_t start, end;

    dim = *d;
    ndata = *n;
    nlambda = *nnlambda;
    double *beta0 = (double*) malloc(dim*sizeof(double));
    double *beta1 = (double*) malloc(dim*sizeof(double));
    int *idx_a = (int*) malloc(dim*sizeof(int)); //active set
    int *idx_i = (int*) malloc(dim*sizeof(int)); //inactive set
    double *XX_diag = (double*) malloc(dim*sizeof(double));
    double *gamma_col = (double*) malloc(dim*sizeof(double));
    double *beta_tild = (double*) malloc(dim*sizeof(double));
    double *alp_tild = (double*) malloc(ndata*sizeof(double));
    double *mu_grad = (double*) malloc(ndata*sizeof(double));
    double *beta_pre = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(ndata*sizeof(double));
    double *mu = (double*) malloc(ndata*sizeof(double));
    double *X_beta = (double*) malloc(ndata*sizeof(double));
    double *y_i = (double*) malloc(ndata*sizeof(double));
    double *Xy = (double*) malloc(dim*sizeof(double));
    double *yXbeta = (double*) malloc(ndata*sizeof(double));
    int *iter_step = (int*) malloc(2*sizeof(int));


    max_ite1 = *max_ite;
    max_ite2 = 1e2;
    max_ite3 = 1e2;
    eps1 = *prec;
    eps2 = 1e-3;
    rho0 = (*rho);
    qrt_n = pow(ndata,1/(*q));
    if(*q < 1e6)
        gamma = 1/(rho0*qrt_n);
    else
        gamma = sqrt(ndata)/(rho0);

    for(j=0; j<dim; j++) {
        beta0[j] = 0;
        idx_i[j] = 1;
        XX_diag[j] = XX[j*dim+j];
    }
    for(j=0; j<ndata; j++) {
        alp[j] = 0;
        mu[j] = 0;
    }
    size_a = 0;
            
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m];
        start = clock();
        max_dif = 1;

        if(*q < 1e6)
            for(j=0;j<dim;j++)
                gamma_col[j] = ilambda/(rho0*XX_diag[j]);
        else
            for(j=0;j<dim;j++)
                gamma_col[j] = ilambda/(rho0*XX_diag[j]);

        ite_ext=0;
        while(max_dif > eps1 && ite_ext < max_ite1){
            // update alpha
            for(j=0; j<ndata; j++){
                X_beta[j]=0;
                for(k=0; k<size_a; k++){
                    w_idx = idx_a[k];
                    X_beta[j]+=X[w_idx*ndata+j]*beta0[w_idx];
                }
                alp_tild[j] = Y[j]-X_beta[j]-mu[j];
            }
            slim(alp, &cc, iter_step, alp_tild, ndata, gamma, *q, 0);
            // update omega
            //ilambda_tmp = ilambda*ndata;//*ndata
            for(j=0; j<ndata; j++){
                y_i[j] = Y[j]-alp[j]-mu[j];
            }
            for(j=0; j<dim; j++){
                Xy[j] = 0;
                beta_pre[j] = beta0[j];
                for(k=0; k<ndata; k++){
                    Xy[j] += X[j*ndata+k]*y_i[k];
                }
            }
            gap_ext = 1;
            ite1 = 0;
            while(gap_ext !=0 && ite1<max_ite2){
                size_a_pre = size_a;
                for(j=0; j<dim; j++){
                    if(idx_i[j] == 1){
                        beta_tild[j] = 0;
                        for(k=0; k<size_a; k++){
                            w_idx = idx_a[k];
                            beta_tild[j] += XX[w_idx*dim+j]*beta0[w_idx];
                        }
                        beta_tild[j] = (Xy[j]-beta_tild[j]+XX_diag[j]*beta0[j])/XX_diag[j];
                        if(j==0 && *intercept==1){
                            beta1[j] = beta_tild[j];
                            idx_a[size_a] = j;
                            size_a++;
                            idx_i[j] = 0;
                        }
                        else {
                            if(fabs(beta_tild[j])<=gamma_col[j]) {
                                beta1[j] = 0;
                            }
                            else{
                                if(beta_tild[j]>gamma_col[j])
                                    beta1[j] = beta_tild[j] - gamma_col[j];
                                else
                                    beta1[j] = beta_tild[j] + gamma_col[j];
                                idx_a[size_a] = j;
                                size_a++;
                                idx_i[j] = 0;
                            }    
                        }
                        beta0[j] = beta1[j];
                    }
                }
                gap_ext = size_a - size_a_pre;
                gap_int = 1;
                ite2 = 0;
                while(gap_int>eps2 && ite2<max_ite3){
                    tmp1 = 0;
                    tmp2 = 0;
                    for(j=0; j<size_a; j++){
                        w_idx = idx_a[j];
                        beta_tild[w_idx] = 0;
                        for(k=0; k<size_a; k++){
                            rs_idx = idx_a[k];
                            beta_tild[w_idx] += XX[rs_idx*dim+w_idx]*beta0[rs_idx];
                        }
                        beta_tild[w_idx] = (Xy[w_idx]-beta_tild[w_idx]+XX_diag[w_idx]*beta0[w_idx])/XX_diag[w_idx];
                        if(w_idx==0 && *intercept==1){
                            beta1[w_idx] = beta_tild[w_idx];
                            tmp2 = tmp2+fabs(beta1[w_idx]);
                        }
                        else {
                            if (fabs(beta_tild[w_idx]) <= gamma_col[w_idx]) {
                                beta1[w_idx] = 0;
                            }
                            else {
                                if (beta_tild[w_idx]>gamma_col[w_idx])
                                    beta1[w_idx] = beta_tild[w_idx] - gamma_col[w_idx];
                                else
                                    beta1[w_idx] = beta_tild[w_idx] + gamma_col[w_idx];
                                tmp2 = tmp2+fabs(beta1[w_idx]);
                            }
                        }
                        beta_dif = beta1[w_idx]-beta0[w_idx];
                        tmp1 = tmp1+fabs(beta_dif);
                        beta0[w_idx] = beta1[w_idx];
                    }
                    gap_int = tmp1/tmp2;
                    ite2++;  
                }
                ite_cnt_int2[m] += ite2;

                junk_a = 0;
                for(j=0; j<size_a; j++){
                    w_idx = idx_a[j];
                    if (beta1[w_idx] == 0){
                        junk_a++;
                        idx_i[w_idx] = 1;
                    }
                    else
                        idx_a[j-junk_a] = w_idx;
                }
                size_a = size_a - junk_a;
                ite1++;
            }
            ite_cnt_int1[m] += ite1;
			
            beta_norm1 = 0;//
            for(j=0;j<size_a;j++){
                w_idx = idx_a[j];
                beta_norm1 += fabs(beta1[w_idx]);//
            }
            F = 0;//
            for(j=0; j<ndata; j++){
                X_beta[j]=0;
                for(k=0; k<size_a; k++){
                    w_idx = idx_a[k];
                    X_beta[j]+=X[w_idx*ndata+j]*beta1[w_idx];
                }
                yXbeta[j] = y_i[j]-X_beta[j];//
                F += yXbeta[j]*yXbeta[j];//
            }
            F = ilambda*beta_norm1+rho0*F/2;//
            
            end = clock();
            runt[m*max_ite1+ite_ext] = (end-start)/ (double)CLOCKS_PER_SEC;
            obj[m*max_ite1+ite_ext] = F;

            beta_sum1 = 0;
            beta_dif_sum1 = 0;
            beta_sum2 = 0;
            for(j=0; j<dim; j++){
                beta_sum1 += fabs(beta0[j]);
                beta_sum2 += fabs(beta_pre[j]);
            }
            beta_out1 = fabs(beta_sum2-beta_sum1)/beta_sum1;

            // update mu
            mu_dif = 0;
            for(j=0; j<ndata; j++){
                //X_beta[j]=0;
                //for(k=0; k<size_a; k++){
                //    w_idx = idx_a[k];
                //    X_beta[j]+=X[w_idx*ndata+j]*beta1[w_idx];
                //}
                mu_grad[j]=alp[j]+X_beta[j]-Y[j];
                mu[j] += mu_grad[j];
                mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
            }
            max_dif = beta_out1>mu_dif ? beta_out1 : mu_dif;
            ite_ext++;
        }
        ite_cnt_ext[m] = ite_ext;

        for(j=0;j<size_a;j++){
            w_idx = idx_a[j];
            beta[m*dim+w_idx] = beta1[w_idx];
        }
    }

    free(beta0);
    free(beta1);
    free(idx_a);
    free(idx_i);
    free(XX_diag);
    free(gamma_col);
    free(beta_tild);
    free(alp_tild);
    free(mu_grad);
    free(beta_pre);
    free(alp);
    free(mu);
    free(X_beta);
    free(yXbeta);
    free(y_i);
    free(Xy);
}
