#include "R.h"
#include "math.h"

void slim_lasso(double *X, double *XX, double *XY, double *beta, int *n, int *d, int *ite_cnt_int1, int *ite_cnt_int2, double *lambda, int * nnlambda, int *max_ite, double *prec, int * intercept)
{
    int j,k,m,ndata,dim,ileft,iright,junk_a,size_a,size_a_pre,w_idx,rs_idx,nlambda;
    int ite1,ite2,gap_ext,max_ite1,max_ite2;
    double gap_int,ilambda,ilambda_tmp,tmp1,tmp2,beta_dif,eps,eps1,eps2,max_dif;

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

    max_ite1 = *max_ite;
    max_ite2 = 1e2;
    eps1 = *prec;
    eps2 = 1e-3;

    for(j=0; j<dim; j++) {
        beta0[j] = 0;
        idx_i[j] = 1;
        XX_diag[j] = XX[j*dim+j];
    }
    size_a = 0;
            
    for(m=0; m<nlambda; m++) {
        ilambda = lambda[m];
        max_dif = 1;

        // update omega
        ilambda_tmp = ilambda*ndata;//*ndata
        for(j=0; j<dim; j++) {
            gamma_col[j] = ilambda_tmp/XX_diag[j];
        }
        gap_ext = 1;
        ite1 = 0;
        while(gap_ext !=0 && ite1<max_ite1){
            size_a_pre = size_a;
            for(j=0; j<dim; j++){
                if(idx_i[j] == 1){
                    beta_tild[j] = 0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        beta_tild[j] += XX[w_idx*dim+j]*beta0[w_idx];
                    }
                    beta_tild[j] = (XY[j]-beta_tild[j]+XX_diag[j]*beta0[j])/XX_diag[j];
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
            while(gap_int>eps2 && ite2<max_ite2){
                tmp1 = 0;
                tmp2 = 0;
                for(j=0; j<size_a; j++){
                    w_idx = idx_a[j];
                    beta_tild[w_idx] = 0;
                    for(k=0; k<size_a; k++){
                        rs_idx = idx_a[k];
                        beta_tild[w_idx] += XX[rs_idx*dim+w_idx]*beta0[rs_idx];
                    }
                    beta_tild[w_idx] = (XY[w_idx]-beta_tild[w_idx]+XX_diag[w_idx]*beta0[w_idx])/XX_diag[w_idx];
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
        ite_cnt_int1[m] = ite1;
            
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
}
