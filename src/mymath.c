#include "mymath.h"

double sign(double x){
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

double max(double x,double y){
    return (x > y) ? x : y;
}

double max_abs_vec(double * x, int n){
    int i;
    double tmp = fabs(x[0]);
    
    for(i=1; i<n ; i++){
        tmp = max(tmp, fabs(x[i]));
    }
    return tmp;
}

double max_vec(double * x, int n){
    int i;
    double tmp = x[0];
    
    for(i=1; i<n ; i++){
        tmp = max(tmp, x[i]);
    }
    return tmp;
}

void max_selc(double *x, double vmax, double *x_s, int n, int *n_s, double z){
    int i,tmp;
    double thresh = vmax-z;

    tmp = 0;
    for(i=0; i<n ; i++){
        if(x[i]>thresh){
            x_s[tmp] = x[i];
            tmp++;
        }
    }
    *n_s = tmp;
}

double min(double x,double y){
    return (x < y) ? x : y;
}

double l1norm(double * x, int n){
    int i;
    double tmp=0;

    for(i=0; i<n; i++) 
        tmp += fabs(x[i]);
    return tmp;
}

void fabs_vc(double *v_in, double *v_out, int n){
    int i;

    for(i=0; i<n; i++)
        v_out[i] = fabs(v_in[i]);
}

void max_fabs_vc(double *v_in, double *v_out, double *vmax, int *n1, int n, double z){
    int i,cnt;
    double tmp, v_abs;

    tmp = 0;
    cnt = 0;
    for(i=0; i<n; i++){
        v_abs = fabs(v_in[i]);
        v_out[i] = v_abs;
        tmp = max(tmp, v_abs);
    }
    *vmax = tmp;
    *n1 = n;
}

void sort_up_bubble(double *v, int n){
    int i,j;
    double tmp;
    int ischanged;

    for(i=n-1; i>=0; i--){ 
        ischanged = 0;
        for(j=0; j<i; j++){
            if(v[j]>v[j+1]){
                tmp = v[j];
                v[j] = v[j+1];
                v[j+1] = tmp;
                ischanged = 1;
            }
        }
        if(ischanged==0) 
            break;
    }
}


void get_residual(double *r, double *y, double *A, double *x, int *xa_idx, int *nn, int *mm)
{
    int i,j,b_idx;
    int n,m;
    double tmp;
    n = *nn;
    m = *mm;
    
    for(i=0;i<n;i++){
        tmp=0;
        for(j=0;j<m;j++){
            b_idx = xa_idx[j];
            tmp+=A[b_idx*n+i]*x[b_idx];
        }
        r[i] = y[i]-tmp;
    }
}

void get_dual(double *u, double *r, double *mmu, int *nn)
{
    int i,n;
    double mu, zv;
    mu = *mmu;
    n = *nn;
    zv = 1;
    for(i=0;i<n;i++){
        u[i] = r[i]/mu;
    }
    euc_proj(u, zv, n); //euclidean projection
}

void get_dual1(double *u, double *r, double *mmu, int *nn)
{
    int i,n;
    double mu, zv;
    mu = *mmu;
    n = *nn;
    zv = 1;
    for(i=0;i<n;i++){
        u[i] = r[i]/mu;
        if(u[i]>zv)
            u[i] = zv;
        if(u[i]<-zv)
            u[i] = -zv;
    }
}

void get_dual2(double *u, double *r, double *mmu, int *nn)
{
    int i,n;
    double mu, zv, tmp_sum;
    mu = *mmu;
    n = *nn;
    zv = 1;
    tmp_sum = 0;
    for(i=0;i<n;i++){
        u[i] = r[i]/mu;
        tmp_sum += u[i]*u[i];
    }
    tmp_sum = sqrt(tmp_sum);
    if(tmp_sum>=zv){
        for(i=0;i<n;i++){
            u[i] = u[i]/tmp_sum;
        }
    }
}

void get_grad(double *g, double *A, double *u, int *dd, int *nn)
{
    int i,j;
    int d,n;
    
    d = *dd;
    n = *nn;
    
    for(i=0;i<d;i++){
        g[i]=0;
        for(j=0;j<n;j++){
            g[i] -= A[i*n+j]*u[j];
        }
    }
}

void get_base(double *base, double *u, double *r, double *mmu, int *nn)
{
    int i,n;
    double mu,tmp;
    mu = *mmu;
    n = *nn; 
    tmp = 0;
    for(i=0;i<n;i++){
        tmp += u[i]*u[i];
    }

    *base = 0;
    for(i=0;i<n;i++){
        *base += u[i]*r[i];
    }
    *base -= mu*tmp/2;
}
