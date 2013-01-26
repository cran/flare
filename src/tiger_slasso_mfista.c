#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "R.h"
#include "mymath.h"

void tiger_slasso_mfista(double *Z, double *icov, double *x, int *d, int *n, double *mu, int *ite_cnt, double *lambda, int * nnlambda, int *max_ite, int *col_cnz, int *row_idx, double *prec, double *L, double *obj, double *runt)
{
    int i,j,k,m,dim,dim1,dim_sq,ndata,nlambda,max_ite1,ite1,ite2,gap_track,y_col,tmp_m,cnz;
    int size_x0, size_y1, size_z1, w_idx;
    double T,T0,imu,cpuTime,u_norm_y,u_norm_z,u_norm_x,y1_norm1_pre,obj_base,munru2,ilambda,Q,Fx,Fz;
    double u_norm2y,u_norm2z,u_norm2x,u_norm1x,x0_norm1,x1_norm1,y1_norm1,z1_norm1;
    double norm_dif,z_dif,x_dif,y1_dif,obj_tmp,eps1,t1,t2,ratio,epsT,tau0;
    clock_t start, end;

    dim1 = *d;
    dim = dim1-1;
    dim_sq = dim1*dim1;
    ndata = *n;
    imu = *mu;
    max_ite1 = *max_ite;
    nlambda = *nnlambda;
    eps1 = *prec;
    ratio = 0.8;
    epsT = 1e-3;
    y1_norm1 = 0;
    cnz = 0;
    double *b = (double*) malloc(ndata*sizeof(double));
    double *A = (double*) malloc(ndata*dim*sizeof(double));
    double *x1 = (double*) malloc(dim*sizeof(double));
    double *x0 = (double*) malloc(dim*sizeof(double));
    double *y1 = (double*) malloc(dim*sizeof(double));
    double *y2 = (double*) malloc(dim*sizeof(double));
    double *z1 = (double*) malloc(dim*sizeof(double));
    double *z0 = (double*) malloc(dim*sizeof(double));
    double *u_y = (double*) malloc(ndata*sizeof(double));
    double *u_z = (double*) malloc(ndata*sizeof(double));
    double *u_x = (double*) malloc(ndata*sizeof(double));
    double *Ax0 = (double*) malloc(ndata*sizeof(double));
    double *bAx0 = (double*) malloc(ndata*sizeof(double));
    double *Ay1 = (double*) malloc(ndata*sizeof(double));
    double *bAy1 = (double*) malloc(ndata*sizeof(double));
    double *Az1 = (double*) malloc(ndata*sizeof(double));
    double *bAz1 = (double*) malloc(ndata*sizeof(double));
    double *g = (double*) malloc(dim*sizeof(double));
    int *idx_col = (int*) malloc(dim*sizeof(int));
    double *idx_x0 = (double*) malloc(dim*sizeof(double));
    double *idx_y1 = (double*) malloc(dim*sizeof(double));
    double *idx_z1 = (double*) malloc(dim*sizeof(double));
    
    for(i=0;i<dim;i++){
        x1[i]=0;
        x0[i]=0;
        y1[i]=0;
        z1[i]=0;
        idx_x0[i]=0;
        idx_y1[i]=0;
        idx_z1[i]=0;
    }
    size_x0 = 0;
    size_y1 = 0;
    size_z1 = 0;
    for(k=0;k<dim1;k++){
        for(i=0;i<ndata;i++){
            b[i] = Z[k*ndata+i];
        }
        y_col = 0;
        for(i=0;i<dim;i++){
            if(y_col == k) y_col++;
            for(j=0;j<ndata;j++){
                A[i*ndata+j] = Z[y_col*ndata+j];
            }
            idx_col[i] = y_col;
            y_col++;
        }
        
        for(m=0;m<nlambda;m++){
            T = (*L)/imu;
            T0 = T;

            t1 = 1;
            ilambda = lambda[m];
    
            start = clock();
            u_norm_y = 0;
            for(i=0;i<ndata;i++){
                Ay1[i]=0;
                for(j=0;j<size_y1;j++){
                    w_idx = idx_y1[j];
                    Ay1[i]+=A[w_idx*ndata+i]*y1[w_idx];
                }
                bAy1[i] = b[i]-Ay1[i];
                u_y[i]=bAy1[i]/imu;
                u_norm_y += u_y[i]*u_y[i];
            }
            u_norm2y = u_norm_y;
            u_norm_y = sqrt(u_norm_y);
            if(u_norm_y>=1){
                u_norm2y = 0;
                for(i=0;i<ndata;i++){
                    u_y[i] = u_y[i]/u_norm_y;
                    u_norm2y += u_y[i]*u_y[i];
                }
            }
            for(i=0;i<dim;i++){
                g[i]=0;
                for(j=0;j<ndata;j++){
                    g[i]-=A[i*ndata+j]*u_y[j];
                }
            }
            obj_base = 0;
            for(i=0;i<ndata;i++){
                obj_base += u_y[i]*bAy1[i];
            }
            obj_base += imu*u_norm2y;
            gap_track = 1;
            while(gap_track == 1 && T>epsT){
                Q = obj_base;
                z1_norm1 = 0;
                size_z1 = 0;
                for(i=0;i<dim;i++){
                    z1[i] = y1[i]-g[i]/T;
                    z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T,0);
                    z_dif = z1[i] - y1[i];
                    Q += g[i]*z_dif + T*z_dif*z_dif/2;
                    if(z1[i] != 0){
                        idx_z1[size_z1]=i;
                        size_z1++;
                        z1_norm1 += fabs(z1[i]);
                    }
                }
                u_norm_z = 0;
                for(i=0;i<ndata;i++){
                    Az1[i]=0;
                   for(j=0;j<size_z1;j++){
                        w_idx=idx_z1[j];
                        Az1[i]+=A[w_idx*ndata+i]*z1[w_idx];
                    }
                    bAz1[i] = b[i]-Az1[i];
                    u_z[i]=bAz1[i]/imu;
                    u_norm_z += u_z[i]*u_z[i];
                }
                Q += ilambda*z1_norm1;
                u_norm2z = u_norm_z;
                u_norm_z = sqrt(u_norm_z);
                if(u_norm_z>=1){
                    u_norm2z = 0;
                    for(i=0;i<ndata;i++){
                        u_z[i] = u_z[i]/u_norm_z;
                        u_norm2z += u_z[i]*u_z[i];
                    }
                }
                Fz = 0;
                for(i=0;i<ndata;i++){
                    Fz += u_z[i]*bAz1[i];
                }
                Fz += imu*u_norm2z+ilambda*z1_norm1;
                if(Fz<Q) T = T*ratio;
                else {
                    T = T/ratio;
                    gap_track = 0;
                }
            }

            z1_norm1 = 0;
            size_z1 = 0;
            for(i=0;i<dim;i++){
                z1[i] = y1[i]-g[i]/(T/ratio);
                z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/(T/ratio),0);
                if(z1[i] != 0){
                    idx_z1[size_z1]=i;
                    size_z1++;
                    z1_norm1 += fabs(z1[i]);
                }
            }
            t2 = (1+sqrt(1+4*t1*t1))/2;
            u_norm_z = 0;
            for(i=0;i<ndata;i++){
                Az1[i]=0;
                for(j=0;j<size_z1;j++){
                    w_idx=idx_z1[j];
                    Az1[i]+=A[w_idx*ndata+i]*z1[w_idx];
                }
                bAz1[i] = b[i]-Az1[i];
                u_z[i]=bAz1[i]/imu;
                u_norm_z += u_z[i]*u_z[i];
            }
            u_norm2z = u_norm_z;
            u_norm_z = sqrt(u_norm_z);
            if(u_norm_z>=1){
                u_norm2z = 0;
                for(i=0;i<ndata;i++){
                    u_z[i] = u_z[i]/u_norm_z;
                    u_norm2z += u_z[i]*u_z[i];
                }
            }
            Fz = 0;
            for(i=0;i<ndata;i++){
                Fz += u_z[i]*bAz1[i];
            }
            Fz += imu*u_norm2z+ilambda*z1_norm1;

            u_norm_x = 0;
            for(i=0;i<ndata;i++){
                Ax0[i]=0;
                for(j=0;j<size_x0;j++){
                    w_idx=idx_x0[j];
                    Ax0[i]+=A[w_idx*ndata+i]*x0[w_idx];
                }
                bAx0[i] = b[i]-Ax0[i];
                u_x[i]=bAx0[i]/imu;
                u_norm_x += u_x[i]*u_x[i];
            }
            u_norm2x = u_norm_x;
            u_norm_x = sqrt(u_norm_x);
            if(u_norm_x>=1){
                u_norm2x = 0;
                for(i=0;i<ndata;i++){
                    u_x[i] = u_x[i]/u_norm_x;
                    u_norm2x += u_x[i]*u_x[i];
                }
            }
            Fx = 0;
            for(i=0;i<ndata;i++){
                Fx += u_x[i]*bAx0[i];
            }
            x0_norm1 = 0;
            for(i=0;i<dim;i++){
                x0_norm1 += fabs(x0[i]);
            }
            Fx += imu*u_norm2x+ilambda*x0_norm1;
        
            if(Fx>Fz){
                for(i=0;i<dim;i++)
                    x1[i] = z1[i];
            }
            else {
                for(i=0;i<dim;i++)
                    x1[i] = x0[i];
            }
        
            size_y1=0;
            size_x0=0;
            for(i=0;i<dim;i++){
                y2[i] = x1[i]+(x1[i]-x0[i])*(t1-1)/t2+(z1[i]-x1[i])*t1/t2;
                z0[i] = z1[i];
                x0[i] = x1[i];
                if(x0[i]!=0){
                    idx_x0[size_x0]=i;
                    size_x0++;
                }
                y1[i] = y2[i];
                if(y1[i]!=0){
                    idx_y1[size_y1]=i;
                    size_y1++;
                }
            }
            t1 = t2;
//if(m==0)
//printf("Q=%f,F=%f,T=%f,T0=%f \n",Q,F,T,T0);
            ite1=0;
            y1_dif = 1;
            while(y1_dif>eps1 && ite1<max_ite1){
                y1_norm1_pre = y1_norm1;
                u_norm_y = 0;
                for(i=0;i<ndata;i++){
                    Ay1[i]=0;
                    for(j=0;j<size_y1;j++){
                        w_idx = idx_y1[j];
                        Ay1[i]+=A[w_idx*ndata+i]*y1[w_idx];
                    }
                    bAy1[i] = b[i]-Ay1[i];
                    u_y[i]=bAy1[i]/imu;
                    u_norm_y += u_y[i]*u_y[i];
                }
                u_norm2y = u_norm_y;
                u_norm_y = sqrt(u_norm_y);
                if(u_norm_y>=1){
                    u_norm2y = 0;
                    for(i=0;i<ndata;i++){
                        u_y[i] = u_y[i]/u_norm_y;
                        u_norm2y += u_y[i]*u_y[i];
                    }
                }
                for(i=0;i<dim;i++){
                    g[i]=0;
                    for(j=0;j<ndata;j++){
                        g[i]-=A[i*ndata+j]*u_y[j];
                    }
                }
        
                ite2=0;
                if(T<T0){
                    obj_base = 0;
                    for(i=0;i<ndata;i++){
                        obj_base += u_y[i]*bAy1[i];
                    }
                    obj_base += imu*u_norm2y;
                    gap_track = 1;
                    while(gap_track == 1){
                        Q = obj_base;
                        z1_norm1 = 0;
                        size_z1 = 0;
                        for(i=0;i<dim;i++){
                            z1[i] = y1[i]-g[i]/T;
                            z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T,0);
                            z_dif = z1[i] - y1[i];
                            Q += g[i]*z_dif + T*z_dif*z_dif/2;
                            if(z1[i] != 0){
                                idx_z1[size_z1]=i;
                                size_z1++;
                                z1_norm1 += fabs(z1[i]);
                            }
                        }
                        Q += ilambda*z1_norm1;
                        u_norm_z = 0;
                        for(i=0;i<ndata;i++){
                            Az1[i]=0;
                            for(j=0;j<size_z1;j++){
                                w_idx=idx_z1[j];
                                Az1[i]+=A[w_idx*ndata+i]*z1[w_idx];
                            }
                            bAz1[i] = b[i]-Az1[i];
                            u_z[i]=bAz1[i]/imu;
                            u_norm_z += u_z[i]*u_z[i];
                        }
                        u_norm2z = u_norm_z;
                        u_norm_z = sqrt(u_norm_z);
                        if(u_norm_z>=1){
                            u_norm2z = 0;
                            for(i=0;i<ndata;i++){
                                u_z[i] = u_z[i]/u_norm_z;
                                u_norm2z += u_z[i]*u_z[i];
                            }
                        }
                        Fz = 0;
                        for(i=0;i<ndata;i++){
                            Fz += u_z[i]*bAz1[i];
                        }
                        Fz += ilambda*z1_norm1 + imu*u_norm2z;
                        if(Fz>Q) T = T/ratio;
                        else gap_track = 0;
                        ite2++;
                    }
                }
                else {
                    z1_norm1 = 0;
                    size_z1 = 0;
                    for(i=0;i<dim;i++){
                        z1[i] = y1[i]-g[i]/T0;
                        z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T0,0);
                        if(z1[i] != 0){
                            idx_z1[size_z1]=i;
                            size_z1++;
                            z1_norm1 += fabs(z1[i]);
                        }
                    }
                }

                t2 = (1+sqrt(1+4*t1*t1))/2;
                u_norm_z = 0;
                for(i=0;i<ndata;i++){
                    Az1[i]=0;
                    for(j=0;j<size_z1;j++){
                        w_idx=idx_z1[j];
                        Az1[i]+=A[w_idx*ndata+i]*z1[w_idx];
                    }
                    bAz1[i] = b[i]-Az1[i];
                    u_z[i]=bAz1[i]/imu;
                    u_norm_z += u_z[i]*u_z[i];
                }
                u_norm2z = u_norm_z;
                u_norm_z = sqrt(u_norm_z);
                if(u_norm_z>=1){
                    u_norm2z = 0;
                    for(i=0;i<ndata;i++){
                        u_z[i] = u_z[i]/u_norm_z;
                        u_norm2z += u_z[i]*u_z[i];
                    }
                }
                Fz = 0;
                for(i=0;i<ndata;i++){
                    Fz += u_z[i]*bAz1[i];
                }
                Fz += imu*u_norm2z + ilambda*z1_norm1;

                u_norm_x = 0;
                for(i=0;i<ndata;i++){
                    Ax0[i]=0;
                    for(j=0;j<size_x0;j++){
                        w_idx=idx_x0[j];
                        Ax0[i]+=A[w_idx*ndata+i]*x0[w_idx];
                    }
                    bAx0[i] = b[i]-Ax0[i];
                    u_x[i]=bAx0[i]/imu;
                    u_norm_x += u_x[i]*u_x[i];
                }
                u_norm2x = u_norm_x;
                u_norm_x = sqrt(u_norm_x);
                if(u_norm_x>=1){
                    u_norm2x = 0;
                    for(i=0;i<ndata;i++){
                        u_x[i] = u_x[i]/u_norm_x;
                        u_norm2x += u_x[i]*u_x[i];
                    }
                }
                Fx = 0;
                for(i=0;i<ndata;i++){
                    Fx += u_x[i]*bAx0[i];
                }
                x0_norm1 = 0;
                for(i=0;i<dim;i++){
                    x0_norm1 += fabs(x0[i]);
                }
                Fx += imu*u_norm2x + ilambda*x0_norm1;
        
                if(Fx>Fz){
                    for(i=0;i<dim;i++)
                        x1[i] = z1[i];
                }
                else {
                    for(i=0;i<dim;i++)
                        x1[i] = x0[i];
                }
            
                y1_norm1 = 0;
                x1_norm1 = 0;
                size_x0 = 0;
                size_y1 = 0;
                for(i=0;i<dim;i++){
                    y2[i] = x1[i]+(x1[i]-x0[i])*(t1-1)/t2+(z1[i]-x1[i])*t1/t2;
                    z0[i] = z1[i];
                    x0[i] = x1[i];
                    if(x0[i]!=0){
                        idx_x0[size_x0]=i;
                        size_x0++;
                    }
                    y1[i] = y2[i];
                    if(y1[i]!=0){
                        idx_y1[size_y1]=i;
                        size_y1++;
                    }
                    y1_norm1 += fabs(y1[i]);
                    x1_norm1 += fabs(x1[i]);
                }
                y1_dif = fabs(y1_norm1 - y1_norm1_pre);
                t1 = t2;

                obj_tmp = 0;
                for(i=0;i<ndata;i++){
                    Ax0[i]=0;
                    for(j=0;j<size_x0;j++){
                        w_idx=idx_x0[j];
                        Ax0[i]+=A[w_idx*ndata+i]*x0[w_idx];
                    }
                    bAx0[i] = b[i]-Ax0[i];
                    obj_tmp += bAx0[i]*bAx0[i];
                }
                if(k==0){
                    end = clock();
                    runt[m*max_ite1+ite1] = (end-start)/ (double)CLOCKS_PER_SEC;
                    obj[m*max_ite1+ite1] = sqrt(obj_tmp)+ilambda*x1_norm1;
                }

//if(m==0 && k==0 && ite1%1==0)// && ite1%20==0
//printf("ite_ext=%d,Q=%f,F=%f,obj2=%f,T=%f,T0=%f \n",ite1,Q,Fz,obj[m*max_ite1+ite1],T,T0);
                ite1++;
            }
            //tau0=sqrt(obj_tmp/ndata);
            //omg[tmp_m+i] = 1/(tau0*tau0);
            tmp_m = m*dim_sq+k*dim1;
            icov[tmp_m+k] = ndata/obj_tmp;
            for(j=0;j<size_x0;j++){
                w_idx=idx_x0[j];
                icov[tmp_m+idx_col[w_idx]] = -x0[w_idx]*icov[tmp_m+k];
                x[cnz] = icov[tmp_m+idx_col[w_idx]];
                row_idx[cnz] = m*dim1+idx_col[w_idx];
                cnz++;
            }
            
            ite_cnt[m*dim1+k] = ite1;
        }
        col_cnz[k+1]=cnz;
    }
    //while (getchar() != '\n');
    free(b);
    free(A);
    free(x0);
    free(x1);
    free(y1);
    free(y2);
    free(z0);
    free(z1);
    free(u_x);
    free(u_y);
    free(u_z);
    free(Ax0);
    free(bAx0);
    free(Ay1);
    free(bAy1);
    free(Az1);
    free(bAz1);
    free(g);
    free(idx_col);
    free(idx_x0);
    free(idx_y1);
    free(idx_z1);
}

