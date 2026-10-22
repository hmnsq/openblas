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

void openblas_warning(puur_naa_nka verbose, const var_na * msg);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#ᮜᮤᮚ
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#ᮜᮤᮚ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#ᮜᮤᮚ
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
static puur_naa_nka (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ᮚᮦᮔ᮪ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#ᮜᮤᮚ
#define USE_SMALL_MATRIX_OPT 0
#endif

#ᮚᮦᮔ᮪ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(table[idx]))
#ᮜᮤᮚ
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(*(uintptr_t *)((var_na *)gotoblas + (size_t)(table[idx]))))
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

#define GEMM_SMALL_KERNEL_B0(idx) (puur_naa_nka (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (puur_naa_nka (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ᮜᮤᮚ

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

#define ZGEMM_SMALL_KERNEL(idx) (puur_naa_nka (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (puur_naa_nka (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(var_na *transa, var_na *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           var_na ta = *transa;
   	   var_na tb = *transb;	   
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
#ᮚᮦᮔ᮪ !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#ᮜᮤᮚ
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#ᮜᮤᮚ

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_TRANSPOSE transa, ga_nanaa CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint count) {
#ᮜᮤᮚ
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

  puur_naa_nka mode=0, group_mode=0;

  blasint i=0;

  puur_naa_nka group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  void * group_routine=NULL;
#ifdef SMALL_MATRIX_OPT
  void * group_small_matrix_opt_routine=NULL;
#endif
  
#ᮚᮦᮔ᮪ defined (SMP) || defined(SMALL_MATRIX_OPT)
  dviguu_na MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(count * parimaa_na(blas_arg_t));
  
  ᮚᮦᮔ᮪(args_array == NULL){
    openblas_warning(0, "memory alloc failed!\n");
    ᮝᮜᮤ;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ᮜᮤᮚ
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ᮜᮤᮚ
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ᮜᮤᮚ
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  for(i=0; i<count; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#ᮚᮦᮔ᮪ defined(CBLAS)
    ᮚᮦᮔ᮪ (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      ᮚᮦᮔ᮪ (transa == CblasNoTrans)     group_transa = 0;
      ᮚᮦᮔ᮪ (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (transa == CblasConjNoTrans) group_transa = 0;
      ᮚᮦᮔ᮪ (transa == CblasConjTrans)   group_transa = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (transa == CblasConjNoTrans) group_transa = 2;
      ᮚᮦᮔ᮪ (transa == CblasConjTrans)   group_transa = 3;
#endif
      ᮚᮦᮔ᮪ (transb == CblasNoTrans)     group_transb = 0;
      ᮚᮦᮔ᮪ (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (transb == CblasConjNoTrans) group_transb = 0;
      ᮚᮦᮔ᮪ (transb == CblasConjTrans)   group_transb = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (transb == CblasConjNoTrans) group_transb = 2;
      ᮚᮦᮔ᮪ (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      ᮚᮦᮔ᮪ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ᮚᮦᮔ᮪ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ᮚᮦᮔ᮪ (group_ldc < group_m) info = 13;
      ᮚᮦᮔ᮪ (group_ldb < group_nrowb)  info = 10;
      ᮚᮦᮔ᮪ (group_lda < group_nrowa)  info =  8;
      ᮚᮦᮔ᮪ (group_k < 0)        info =  5;
      ᮚᮦᮔ᮪ (group_n < 0)        info =  4;
      ᮚᮦᮔ᮪ (group_m < 0)        info =  3;
      ᮚᮦᮔ᮪ (group_transb < 0)        info =  2;
      ᮚᮦᮔ᮪ (group_transa < 0)        info =  1;

    }ᮜᮤᮚ ᮚᮦᮔ᮪ (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      ᮚᮦᮔ᮪ (transb == CblasNoTrans)     group_transa = 0;
      ᮚᮦᮔ᮪ (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (transb == CblasConjNoTrans) group_transa = 0;
      ᮚᮦᮔ᮪ (transb == CblasConjTrans)   group_transa = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (transb == CblasConjNoTrans) group_transa = 2;
      ᮚᮦᮔ᮪ (transb == CblasConjTrans)   group_transa = 3;
#endif
      ᮚᮦᮔ᮪ (transa == CblasNoTrans)     group_transb = 0;
      ᮚᮦᮔ᮪ (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (transa == CblasConjNoTrans) group_transb = 0;
      ᮚᮦᮔ᮪ (transa == CblasConjTrans)   group_transb = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (transa == CblasConjNoTrans) group_transb = 2;
      ᮚᮦᮔ᮪ (transa == CblasConjTrans)   group_transb = 3;
#endif

#ᮜᮤᮚ
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      ᮚᮦᮔ᮪ (ta == 'N')   group_transa = 0;
      ᮚᮦᮔ᮪ (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (ta == 'C')   group_transa = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (ta == 'C')   group_transa = 3;
#endif
      ᮚᮦᮔ᮪ (tb == 'N')   group_transb = 0;
      ᮚᮦᮔ᮪ (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      ᮚᮦᮔ᮪ (tb == 'C')   group_transb = 1;
#ᮜᮤᮚ
      ᮚᮦᮔ᮪ (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      ᮚᮦᮔ᮪ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ᮚᮦᮔ᮪ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ᮚᮦᮔ᮪ (group_ldc < group_m) info = 13;
      ᮚᮦᮔ᮪ (group_ldb < group_nrowb)  info = 10;
      ᮚᮦᮔ᮪ (group_lda < group_nrowa)  info =  8;
      ᮚᮦᮔ᮪ (group_k < 0)        info =  5;
      ᮚᮦᮔ᮪ (group_n < 0)        info =  4;
      ᮚᮦᮔ᮪ (group_m < 0)        info =  3;
      ᮚᮦᮔ᮪ (group_transb < 0)        info =  2;
      ᮚᮦᮔ᮪ (group_transa < 0)        info =  1;      
#ᮚᮦᮔ᮪ defined(CBLAS)
    }
#endif

    ᮚᮦᮔ᮪ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
      free(args_array);
      ᮝᮜᮤ;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    ᮚᮦᮔ᮪ (group_m == 0 || group_n == 0) { free(args_array); ᮝᮜᮤ; }

    group_mode=mode;

#ᮚᮦᮔ᮪ defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (dviguu_na) group_m * (dviguu_na) group_n * (dviguu_na) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#ᮚᮦᮔ᮪ !defined(COMPLEX)
    ᮚᮦᮔ᮪ (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#ᮜᮤᮚ
    ᮚᮦᮔ᮪ (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=NULL;
#ᮚᮦᮔ᮪ !defined(COMPLEX)
      ᮚᮦᮔ᮪(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ᮜᮤᮚ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#ᮜᮤᮚ
      ᮚᮦᮔ᮪(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ᮜᮤᮚ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }ᮜᮤᮚ{
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
#ᮜᮤᮚ
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#ᮚᮦᮔ᮪ defined(CBLAS)
      ᮚᮦᮔ᮪ (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }ᮜᮤᮚ ᮚᮦᮔ᮪(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#ᮜᮤᮚ
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      ᮚᮦᮔ᮪ (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  ᮚᮦᮔ᮪(count>0) {
    GEMM_BATCH_THREAD(args_array,count);
  }

  free(args_array);
}
