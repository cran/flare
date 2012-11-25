#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include <time.h>

//extern "C" void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);

void svd1(int *mm,     // number of rows in matrix
          int *nn,     // number of columns in matrix
          int *lldu,   // leading dimension of matrix
          int *lldvt,   // leading dimension of matrix
          double *a,
          double *u,
          double *vt,
          double *s) // pointer to top-left corner
{
    char jobu, jobvt;
    // Setup a buffer to hold the singular values:
    int sv_n,i,m,n,ldu,ldvt;

    // Workspace and status variables:
    double workSize;
    int lwork = -1;
    time_t t1,t2,t3,t4;
    double dft1, dft2;

    m=*mm;
    n=*nn;
    sv_n = m < n ? m : n;
    int *iwork = (int*) malloc(8*sv_n*sizeof(int));
    int info = 0;

    double *work = (double*) malloc(1*sizeof(double));
    
    jobu = 'S';
    jobvt = 'S';
    ldu=*lldu;
    ldvt=*lldvt;
    // Call dgesdd_ with lwork = -1 to query optimal workspace size:
    //dgesvd(A, m, n, s, u, vt);
//time(&t1);
    //dgesvd_(&jobu, &jobvt, &m, &n, a, &ldu, s, u,&ldu, vt, &ldvt, work, &lwork, &info);
    dgesdd_(&jobu, &m, &n, a, &ldu, s, u,&ldu, vt, &ldvt, work, &lwork, iwork, &info);
//time(&t2);
//dft1=difftime(t2,t1);
//printf("runt=%f sec \n",dft1);
    if (info != 0) //{printf("error1\n");}// handle error conditions here

    // Optimal workspace size is returned in work[0].
    lwork = work[0];
    //lwork = 5*sv_n;
    work = (double*) malloc(lwork * sizeof(double));
    // Call dgesdd_ to do the actual computation:
    //dgesvd(A, m, n, s, u, vt);
//time(&t3);
    //dgesvd_(&jobu, &jobvt, &m, &n, a, &ldu, s, u,&ldu, vt, &ldvt, work, &lwork, &info);
    dgesdd_(&jobu, &m, &n, a, &ldu, s, u,&ldu, vt, &ldvt, work, &lwork, iwork, &info);
//for(i=0;i<1e9;i++)
//time(&t4);
//dft2=difftime(t4,t3);
//printf("runt=%f sec \n",dft2);
    if (info) //{printf("error2\n");} // handle error conditions here

    // Cleanup workspace:
    free(work);
    free(iwork);
}
