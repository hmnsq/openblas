/*****************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the OpenBLAS project nor the names of 
      its contributors may be used to endorse or promote products 
      derived from this software without specific prior written 
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "common.h"

void openblas_warning(số_nguyên verbose, const char * msg);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#không_thì
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#không_thì
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#không_thì
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
static số_nguyên (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#nếu defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#không_thì
#define USE_SMALL_MATRIX_OPT 0
#endif

#nếu USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(table[idx]))
#không_thì
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(*(uintptr_t *)((char *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
static size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


static size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#không_thì

static size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

static size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (số_nguyên (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(char *transa, char *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           char ta = *transa;
   	   char tb = *transb;	   
	   blasint count = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#nếu !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#không_thì
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#không_thì

void CNAME(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE transa, enum CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint count) {
#không_thì
	   void * valpha,
	   void * va, blasint lda, blasint stridea,
	   void * vb, blasint ldb, blasint strideb,
	   void * vbeta,
	   void * vc, blasint ldc, blasint stridec, blasint count) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=NULL;

  số_nguyên mode=0, group_mode=0;

  blasint i=0;

  số_nguyên group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  void * group_routine=NULL;
#ifdef SMALL_MATRIX_OPT
  void * group_small_matrix_opt_routine=NULL;
#endif
  
#nếu defined (SMP) || defined(SMALL_MATRIX_OPT)
  double MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(count * sizeof(blas_arg_t));
  
  nếu(args_array == NULL){
    openblas_warning(0, "memory alloc failed!\n");
    trả_về;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#không_thì
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#không_thì
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#không_thì
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  cho(i=0; i<count; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#nếu defined(CBLAS)
    nếu (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      nếu (transa == CblasNoTrans)     group_transa = 0;
      nếu (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      nếu (transa == CblasConjNoTrans) group_transa = 0;
      nếu (transa == CblasConjTrans)   group_transa = 1;
#không_thì
      nếu (transa == CblasConjNoTrans) group_transa = 2;
      nếu (transa == CblasConjTrans)   group_transa = 3;
#endif
      nếu (transb == CblasNoTrans)     group_transb = 0;
      nếu (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      nếu (transb == CblasConjNoTrans) group_transb = 0;
      nếu (transb == CblasConjTrans)   group_transb = 1;
#không_thì
      nếu (transb == CblasConjNoTrans) group_transb = 2;
      nếu (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      nếu (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      nếu (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      nếu (group_ldc < group_m) info = 13;
      nếu (group_ldb < group_nrowb)  info = 10;
      nếu (group_lda < group_nrowa)  info =  8;
      nếu (group_k < 0)        info =  5;
      nếu (group_n < 0)        info =  4;
      nếu (group_m < 0)        info =  3;
      nếu (group_transb < 0)        info =  2;
      nếu (group_transa < 0)        info =  1;

    }không_thì nếu (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      nếu (transb == CblasNoTrans)     group_transa = 0;
      nếu (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      nếu (transb == CblasConjNoTrans) group_transa = 0;
      nếu (transb == CblasConjTrans)   group_transa = 1;
#không_thì
      nếu (transb == CblasConjNoTrans) group_transa = 2;
      nếu (transb == CblasConjTrans)   group_transa = 3;
#endif
      nếu (transa == CblasNoTrans)     group_transb = 0;
      nếu (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      nếu (transa == CblasConjNoTrans) group_transb = 0;
      nếu (transa == CblasConjTrans)   group_transb = 1;
#không_thì
      nếu (transa == CblasConjNoTrans) group_transb = 2;
      nếu (transa == CblasConjTrans)   group_transb = 3;
#endif

#không_thì
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      nếu (ta == 'N')   group_transa = 0;
      nếu (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      nếu (ta == 'C')   group_transa = 1;
#không_thì
      nếu (ta == 'C')   group_transa = 3;
#endif
      nếu (tb == 'N')   group_transb = 0;
      nếu (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      nếu (tb == 'C')   group_transb = 1;
#không_thì
      nếu (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      nếu (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      nếu (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      nếu (group_ldc < group_m) info = 13;
      nếu (group_ldb < group_nrowb)  info = 10;
      nếu (group_lda < group_nrowa)  info =  8;
      nếu (group_k < 0)        info =  5;
      nếu (group_n < 0)        info =  4;
      nếu (group_m < 0)        info =  3;
      nếu (group_transb < 0)        info =  2;
      nếu (group_transa < 0)        info =  1;      
#nếu defined(CBLAS)
    }
#endif

    nếu (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
      free(args_array);
      trả_về;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    nếu (group_m == 0 || group_n == 0) { free(args_array); trả_về; }

    group_mode=mode;

#nếu defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (double) group_m * (double) group_n * (double) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#nếu !defined(COMPLEX)
    nếu (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#không_thì
    nếu (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=NULL;
#nếu !defined(COMPLEX)
      nếu(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }không_thì{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#không_thì
      nếu(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }không_thì{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }không_thì{
#endif
      group_routine=(void*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
      args_array[i].m=group_m;
      args_array[i].n=group_n;
      args_array[i].k=group_k;
      args_array[i].lda=group_lda;
      args_array[i].ldb=group_ldb;
      args_array[i].ldc=group_ldc;
#ifdef COMPLEX
      args_array[i].alpha=alpha;
      args_array[i].beta=beta;
#không_thì
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#nếu defined(CBLAS)
      nếu (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }không_thì nếu(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#không_thì
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      nếu (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  nếu(count>0) {
    GEMM_BATCH_THREAD(args_array,count);
  }

  free(args_array);
}
