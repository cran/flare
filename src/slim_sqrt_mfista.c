#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "R.h"
#include "mymath.h"

void slim_sqrt_mfista(double *b, double *A, double *beta, int *n, int *d, double *mu, int *ite_cnt_init, int *ite_cnt_ex, int *ite_cnt_in, double *lambda, int * nnlambda, int *max_ite, double *prec, double *L, int *intercept)
{
    int i,j,m,dim,ndata,nlambda,max_ite1,ite0,ite1,ite2,ite21,gap_track;
    int size_x0, size_y1, size_z1, w_idx;
    double T,T0,imu,cpuTime,opt_dif;
    double x1_norm1,y1_norm1,z1_norm1,y1_norm1_pre,obj_base,ilambda,Q,Fx,Fz;
    double norm_dif,z_dif,x_dif,y1_dif,obj_tmp,eps1,t1,t2,ratio,epsT;
    clock_t start, end;

    dim = *d;
    ndata = *n;
    imu = *mu;
    max_ite1 = *max_ite;
    nlambda = *nnlambda;
    eps1 = *prec;
    ratio = 0.8;
    epsT = 1e-3;
    y1_norm1 = 0;
    double *x1 = (double*) malloc(dim*sizeof(double));
    double *x0 = (double*) malloc(dim*sizeof(double));
    double *y1 = (double*) malloc(dim*sizeof(double));
    double *y2 = (double*) malloc(dim*sizeof(double));
    double *z1 = (double*) malloc(dim*sizeof(double));
    double *z0 = (double*) malloc(dim*sizeof(double));
    double *u_y = (double*) malloc(ndata*sizeof(double));
    double *u_z = (double*) malloc(ndata*sizeof(double));
    double *u_x = (double*) malloc(ndata*sizeof(double));
    double *bAx0 = (double*) malloc(ndata*sizeof(double));
    double *bAy1 = (double*) malloc(ndata*sizeof(double));
    double *bAz1 = (double*) malloc(ndata*sizeof(double));
    double *g = (double*) malloc(dim*sizeof(double));
    int *idx_x0 = (int*) malloc(dim*sizeof(int));
    int *idx_y1 = (int*) malloc(dim*sizeof(int));
    int *idx_z1 = (int*) malloc(dim*sizeof(int));
    
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
    for(m=0;m<nlambda;m++){
        T = (*L)/imu;
        T0 = T;

        t1 = 1;
        ilambda = lambda[m];
    
        get_residual(bAy1, b, A, y1, idx_y1, &ndata, &size_y1); // bAy1=b-A*y1
        get_dual2(u_y, bAy1, &imu, &ndata); //u_y=proj(bAy1)
        get_grad(g, A, u_y, &dim, &ndata); //g=-A*u_y
        get_base(&obj_base, u_y, bAy1, &imu, &ndata); //obj_base=u_y*bAy1-imu*||u_y||_2^2/2

        gap_track = 1;
        ite0 = 0;
        while(gap_track == 1 && T>epsT){
            Q = obj_base;
            z1_norm1 = 0;
            size_z1 = 0;
            for(i=0;i<dim;i++){
                z1[i] = y1[i]-g[i]/T;
                if(i>0 && *intercept==1 || *intercept==0){
                    z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T,0);
                }
                z_dif = z1[i] - y1[i];
                Q += g[i]*z_dif + T*z_dif*z_dif/2;
                if(z1[i] != 0){
                    idx_z1[size_z1]=i;
                    size_z1++;
                    z1_norm1 += fabs(z1[i]);
                }
            }
            Q += ilambda*z1_norm1;
            
            get_residual(bAz1, b, A, z1, idx_z1, &ndata, &size_z1); // bAz1=b-A*z1
            get_dual2(u_z, bAz1, &imu, &ndata); //u_z=proj(bAz1)
            get_base(&Fz, u_z, bAz1, &imu, &ndata); //obj_base=u_z*bAz1-imu*||u_z||_2^2/2
            
            Fz += z1_norm1*ilambda;
            ite0++;
            if(Fz<Q) T = T*ratio;
            else {
                T = T/ratio;
                if(ite0>1)
                    gap_track = 0;
            }
        }

        z1_norm1 = 0;
        size_z1 = 0;
        for(i=0;i<dim;i++){
            z1[i] = y1[i]-g[i]/(T/ratio);
            if(i>0 && *intercept==1 || *intercept==0){
                z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/(T/ratio),0);
            }
            if(z1[i] != 0){
                idx_z1[size_z1]=i;
                size_z1++;
                z1_norm1 += fabs(z1[i]);
            }
        }
        t2 = (1+sqrt(1+4*t1*t1))/2;
        get_residual(bAz1, b, A, z1, idx_z1, &ndata, &size_z1); // bAz1=b-A*z1
        get_dual2(u_z, bAz1, &imu, &ndata); //u_z=proj(bAz1)
        get_base(&Fz, u_z, bAz1, &imu, &ndata); //obj_base=u_z*bAz1-imu*||u_z||_2^2/2
        Fz += z1_norm1*ilambda;

        get_residual(bAx0, b, A, x0, idx_x0, &ndata, &size_x0); // bAx0=b-A*x0
        get_dual2(u_x, bAx0, &imu, &ndata); //u_x=proj(bAx0)
        get_base(&Fx, u_x, bAx0, &imu, &ndata); //obj_base=u_x*bAx0-imu*||u_x||_2^2/2
        Fx += l1norm(x0, dim)*ilambda;
        
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
        ite21=0;
        y1_dif = 1;
        opt_dif = 1;
        get_residual(bAy1, b, A, y1, idx_y1, &ndata, &size_y1); // bAy1=b-A*y1
        get_dual2(u_y, bAy1, &imu, &ndata); //u_y=proj(bAy1)
        get_grad(g, A, u_y, &dim, &ndata); //g=-A*u_y
        while(y1_dif>eps1 && ite1<max_ite1){
        //while(opt_dif>eps1 && ite1<max_ite1){
            y1_norm1_pre = y1_norm1;
            ite2=0;
            if(T<T0){
                get_base(&obj_base, u_y, bAy1, &imu, &ndata); //obj_base=u_y*bAy1-imu*||u_y||_2^2/2
                gap_track = 1;
                while(gap_track == 1){
                    Q = obj_base;
                    z1_norm1 = 0;
                    size_z1 = 0;
                    for(i=0;i<dim;i++){
                        z1[i] = y1[i]-g[i]/T;
                        if(i>0 && *intercept==1 || *intercept==0){
                            z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T,0);
                        }
                        z_dif = z1[i] - y1[i];
                        Q += g[i]*z_dif + T*z_dif*z_dif/2;
                        if(z1[i] != 0){
                            idx_z1[size_z1]=i;
                            size_z1++;
                            z1_norm1 += fabs(z1[i]);
                        }
                    }
                    Q += ilambda*z1_norm1;
                    
                    get_residual(bAz1, b, A, z1, idx_z1, &ndata, &size_z1); // bAz1=b-A*z1
                    get_dual2(u_z, bAz1, &imu, &ndata); //u_z=proj(bAz1)
                    get_base(&Fz, u_z, bAz1, &imu, &ndata); //obj_base=u_z*bAz1-imu*||u_z||_2^2/2
                    Fz += z1_norm1*ilambda;
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
                    if(i>0 && *intercept==1 || *intercept==0){
                        z1[i] = sign(z1[i])*max(fabs(z1[i])-ilambda/T0,0);
                    }
                    if(z1[i] != 0){
                        idx_z1[size_z1]=i;
                        size_z1++;
                        z1_norm1 += fabs(z1[i]);
                    }
                }
            }

            t2 = (1+sqrt(1+4*t1*t1))/2;
            
            get_residual(bAz1, b, A, z1, idx_z1, &ndata, &size_z1); // bAz1=b-A*z1
            get_dual2(u_z, bAz1, &imu, &ndata); //u_z=proj(bAz1)
            get_base(&Fz, u_z, bAz1, &imu, &ndata); //obj_base=u_z*bAz1-imu*||u_z||_2^2/2
            Fz += z1_norm1*ilambda;

            get_residual(bAx0, b, A, x0, idx_x0, &ndata, &size_x0); // bAx0=b-A*x0
            get_dual2(u_x, bAx0, &imu, &ndata); //u_x=proj(bAx0)
            get_base(&Fx, u_x, bAx0, &imu, &ndata); //obj_base=u_x*bAx0-imu*||u_x||_2^2/2
            Fx += l1norm(x0, dim)*ilambda;
        
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
            y1_dif = 0;
            for(i=0;i<dim;i++){
                y2[i] = x1[i]+(x1[i]-x0[i])*(t1-1)/t2+(z1[i]-x1[i])*t1/t2;
                z0[i] = z1[i];
                x0[i] = x1[i];
                if(x0[i]!=0){
                    idx_x0[size_x0]=i;
                    size_x0++;
                }
                //y1_dif += fabs(y1[i] - y2[i]);
                y1[i] = y2[i];
                if(y1[i]!=0){
                    idx_y1[size_y1]=i;
                    size_y1++;
                }
                y1_norm1 += fabs(y1[i]);
            }
            y1_dif = fabs(y1_norm1 - y1_norm1_pre);
            t1 = t2;

            get_residual(bAy1, b, A, y1, idx_y1, &ndata, &size_y1); // bAy1=b-A*y1
            get_dual2(u_y, bAy1, &imu, &ndata); //u_y=proj(bAy1)
            get_grad(g, A, u_y, &dim, &ndata); //g=-A*u_y
            ite1++;
            ite21 += ite2;
        }
        for(i=0;i<size_x0;i++){
            w_idx=idx_x0[i];
            beta[m*dim+w_idx] = x0[w_idx];
        }
        ite_cnt_init[m] = ite0;
        ite_cnt_ex[m] = ite1;
        ite_cnt_in[m] = ite21;
    }
    //while (getchar() != '\n');
    free(x0);
    free(x1);
    free(y1);
    free(y2);
    free(z0);
    free(z1);
    free(u_x);
    free(u_y);
    free(u_z);
    free(bAx0);
    free(bAy1);
    free(bAz1);
    free(g);
    free(idx_x0);
    free(idx_y1);
    free(idx_z1);
}
