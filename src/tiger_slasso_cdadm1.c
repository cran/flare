#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "R.h"

void tiger_slasso_cdadm1(double *tau0, int col_i, double lambda, int dim, int ndata, int *idx_col, int *size_a, int *idx_i, int *idx_a, double *omg0, double *Z, double *Z_i, double *SS, double *S_col)
{
    int i,j,k,m,junk_a,size_a_pre,w_idx,rs_idx;
    int ite_ext,ite_int1,ite_int2,gap_ext,max_ite1,max_ite2,max_ite3, cnz, tmp_m;
    double rho,alp_tild_sq,sq_n,alp_th, gap_int,tmp1,tmp2,err1,err2,omg_2norm,mu_2norm,omg_dif,eps1,eps2,omg_temp,alp_dif,mu_dif,max_dif;
    double omg_sum1, omg_sum2, omg_dif_sum1, omg_out1;
    char c;

    rho = sqrt(ndata);
    sq_n = sqrt(ndata);
    double *omg1 = (double*) malloc(dim*sizeof(double));
    double *alp_tild = (double*) malloc(ndata*sizeof(double));
    double *mu_grad = (double*) malloc(ndata*sizeof(double));
    double *omg_grad = (double*) malloc(dim*sizeof(double));
    double *omg_pre = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(ndata*sizeof(double));
    double *mu = (double*) malloc(ndata*sizeof(double));
    double *gamma_col = (double*) malloc(dim*sizeof(double));
    double *Z_omg = (double*) malloc(ndata*sizeof(double));
    double *omg_tild = (double*) malloc(dim*sizeof(double));
    double *y_i = (double*) malloc(ndata*sizeof(double));
    double *Zy = (double*) malloc(dim*sizeof(double));

    max_ite1 = 1e2;
    max_ite2 = 1e2;
    max_ite3 = 1e2;
    eps1 = 1e-3;
    eps2 = 1e-3;
    cnz = 0;

    for(i=0; i<dim; i++){
        gamma_col[i] = (lambda/rho)/S_col[i];
    }
    //for(i=0; i<dim; i++){
        for(j=0; j<ndata; j++) {
            alp[j] = 0;
            mu[j] = 0;
        }
        // idx_i[i] = 0;
            
        //for(m=0; m<nlambda; m++) {
            ite_ext = 0;
            max_dif = 1;
            while(max_dif > eps1 && ite_ext < max_ite1){
                // update alpha

                alp_tild_sq=0;
                for(j=0; j<ndata; j++){
                    Z_omg[j]=0;
                    for(k=0; k<*size_a; k++){
                        w_idx = idx_a[k];
                        Z_omg[j]+=Z[idx_col[w_idx]*ndata+j]*omg0[w_idx];
                    }
                    alp_tild[j] = Z_i[j]-Z_omg[j]-mu[j];
                    alp_tild_sq += alp_tild[j]*alp_tild[j];
                }
                alp_th = 1-1/(rho*sq_n*sqrt(alp_tild_sq));
                if(alp_th<0) alp_th=0;
                for(j=0; j<ndata; j++){
                    alp[j] = alp_tild[j]*alp_th;
                }

                // update omega
                for(j=0; j<ndata; j++){
                    y_i[j] = Z_i[j]-alp[j]-mu[j];
                }
                for(j=0; j<dim; j++){
                    Zy[j] = 0;
                    omg_pre[j] = omg0[j];
                    for(k=0; k<ndata; k++){
                        Zy[j] += Z[idx_col[j]*ndata+k]*y_i[k];
                    }
                }

                gap_ext = 1;
                ite_int1 = 0;
                while(gap_ext !=0 && ite_int1<max_ite2){
                    size_a_pre = *size_a;
                    for(j=0; j<dim; j++){
                        if(idx_i[j] == 1){
                            omg_tild[j] = 0;
                            for(k=0; k<*size_a; k++){
                                w_idx = idx_a[k];
                                omg_tild[j] += SS[w_idx*dim+j]*omg0[w_idx];
                            }
                            omg_tild[j] = (Zy[j]-omg_tild[j]+S_col[j]*omg0[j])/S_col[j];
                            if(fabs(omg_tild[j])<=gamma_col[j]) {
                                omg1[j] = 0;
                            }
                            else{
                                if(omg_tild[j]>gamma_col[j])
                                    omg1[j] = omg_tild[j] - gamma_col[j];
                                else
                                    omg1[j] = omg_tild[j] + gamma_col[j];
                                idx_a[*size_a] = j;
                                (*size_a)++;
                                idx_i[j] = 0;
                            }
                            omg0[j] = omg1[j];
                        }
                    }
                    gap_ext = *size_a - size_a_pre;
                    gap_int = 1;
                    ite_int2 = 0;
                    while(gap_int>eps2 && ite_int2<max_ite3){
                        tmp1 = 0;
                        tmp2 = 0;
                        for(j=0; j<*size_a; j++){
                            w_idx = idx_a[j];
                            omg_tild[w_idx] = 0;
                            for(k=0; k<*size_a; k++){
                                rs_idx = idx_a[k];
                                omg_tild[w_idx] += SS[rs_idx*dim+w_idx]*omg0[rs_idx];
                            }
                            omg_tild[w_idx] = (Zy[w_idx]-omg_tild[w_idx]+S_col[w_idx]*omg0[w_idx])/S_col[w_idx];
                            if (fabs(omg_tild[w_idx]) <= gamma_col[w_idx]) {
                                omg1[w_idx] = 0;
                            }
                            else {
                                if (omg_tild[w_idx]>gamma_col[w_idx])
                                    omg1[w_idx] = omg_tild[w_idx] - gamma_col[w_idx];
                                else
                                    omg1[w_idx] = omg_tild[w_idx] + gamma_col[w_idx];
                                tmp2 = tmp2+fabs(omg1[w_idx]);
                            }
                            omg_dif = omg1[w_idx]-omg0[w_idx];
                            tmp1 = tmp1+fabs(omg_dif);
                            omg0[w_idx] = omg1[w_idx];
                        }
                        gap_int = tmp1/tmp2;
                        ite_int2++;  
                    }
                
                    junk_a = 0;
                    for(j=0; j<*size_a; j++){
                        w_idx = idx_a[j];
                        if (omg1[w_idx] == 0){
                            junk_a++;
                            idx_i[w_idx] = 1;
                        }
                        else
                            idx_a[j-junk_a] = w_idx;
                    }
                    *size_a = *size_a - junk_a;
                    ite_int1++;
                }

                omg_sum1 = 0;
                omg_dif_sum1 = 0;
                omg_sum2 = 0;
                for(j=0; j<ndata; j++){
                    //omg_dif_sum1 += fabs(omg_pre[j] - omg0[j]);
                    omg_sum1 += fabs(omg0[j]);
                    omg_sum2 += fabs(omg_pre[j]);
                    Z_omg[j]=0;
                    for(k=0; k<*size_a; k++){
                        w_idx = idx_a[k];
                        Z_omg[j]+=Z[idx_col[w_idx]*ndata+j]*omg1[w_idx];
                    }
                }
                //omg_out1 = omg_dif_sum1/omg_sum1;
                omg_out1 = fabs(omg_sum2-omg_sum1)/omg_sum1;

                // update mu
                mu_dif = 0;
                for(j=0; j<ndata; j++){
                    mu_grad[j]=alp[j]+Z_omg[j]-Z_i[j];
                    mu[j] += mu_grad[j];
                    mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
                }
                max_dif = omg_out1>mu_dif ? omg_out1 : mu_dif;
                ite_ext++;
            }
        //}
    //}
    *tau0=0;
    for(i=0;i<ndata;i++)
        *tau0 += alp[i]*alp[i];
    *tau0 = sqrt(*tau0/ndata);

    free(omg1);
    free(alp_tild);
    free(mu_grad);
    free(omg_grad);
    free(omg_pre);
    free(alp);
    free(mu);
    free(gamma_col);
    free(Z_omg);
    free(omg_tild);
    free(y_i);
    free(Zy);
}
