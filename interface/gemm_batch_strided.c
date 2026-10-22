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

#അന്തര്ഭാവ <stdio.h>
#അന്തര്ഭാവ <stdlib.h>
#അന്തര്ഭാവ "common.h"

ശൂന്യം openblas_warning(പൂർണ്ണസംഖ്യ verbose, സ്ഥിരം അക്ഷരം * സന്ദേശം);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#അല്ലെങ്കിൽ
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#അല്ലെങ്കിൽ
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
സ്ഥിതം പൂർണ്ണസംഖ്യ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#എങ്കിൽ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#അല്ലെങ്കിൽ
#define USE_SMALL_MATRIX_OPT 0
#endif

#എങ്കിൽ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ശൂന്യം *)(table[idx]))
#അല്ലെങ്കിൽ
#define SMALL_KERNEL_ADDR(table, idx) ((ശൂന്യം *)(*(uintptr_t *)((അക്ഷരം *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
സ്ഥിതം size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


സ്ഥിതം size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#അല്ലെങ്കിൽ

സ്ഥിതം size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

സ്ഥിതം size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (പൂർണ്ണസംഖ്യ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
ശൂന്യം NAME(അക്ഷരം *transa, അക്ഷരം *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           അക്ഷരം ta = *transa;
   	   അക്ഷരം tb = *transb;	   
	   blasint എണ്ണം = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#എങ്കിൽ !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#അല്ലെങ്കിൽ
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#അല്ലെങ്കിൽ

ശൂന്യം CNAME(എണ്ണൽ CBLAS_ORDER order, എണ്ണൽ CBLAS_TRANSPOSE transa, എണ്ണൽ CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint എണ്ണം) {
#അല്ലെങ്കിൽ
	   ശൂന്യം * valpha,
	   ശൂന്യം * va, blasint lda, blasint stridea,
	   ശൂന്യം * vb, blasint ldb, blasint strideb,
	   ശൂന്യം * vbeta,
	   ശൂന്യം * vc, blasint ldc, blasint stridec, blasint എണ്ണം) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ശൂന്യ;

  പൂർണ്ണസംഖ്യ mode=0, group_mode=0;

  blasint i=0;

  പൂർണ്ണസംഖ്യ group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  ശൂന്യം * group_routine=ശൂന്യ;
#ifdef SMALL_MATRIX_OPT
  ശൂന്യം * group_small_matrix_opt_routine=ശൂന്യ;
#endif
  
#എങ്കിൽ defined (SMP) || defined(SMALL_MATRIX_OPT)
  ഇരട്ട MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(എണ്ണം * വലുപ്പം(blas_arg_t));
  
  എങ്കിൽ(args_array == ശൂന്യ){
    openblas_warning(0, "memory alloc failed!\n");
    തിരികെ;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#അല്ലെങ്കിൽ
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#അല്ലെങ്കിൽ
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#അല്ലെങ്കിലെങ്കിൽ defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#അല്ലെങ്കിൽ
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  വേണ്ടി(i=0; i<എണ്ണം; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#എങ്കിൽ defined(CBLAS)
    എങ്കിൽ (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      എങ്കിൽ (transa == CblasNoTrans)     group_transa = 0;
      എങ്കിൽ (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      എങ്കിൽ (transa == CblasConjNoTrans) group_transa = 0;
      എങ്കിൽ (transa == CblasConjTrans)   group_transa = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (transa == CblasConjNoTrans) group_transa = 2;
      എങ്കിൽ (transa == CblasConjTrans)   group_transa = 3;
#endif
      എങ്കിൽ (transb == CblasNoTrans)     group_transb = 0;
      എങ്കിൽ (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      എങ്കിൽ (transb == CblasConjNoTrans) group_transb = 0;
      എങ്കിൽ (transb == CblasConjTrans)   group_transb = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (transb == CblasConjNoTrans) group_transb = 2;
      എങ്കിൽ (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      എങ്കിൽ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      എങ്കിൽ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      എങ്കിൽ (group_ldc < group_m) info = 13;
      എങ്കിൽ (group_ldb < group_nrowb)  info = 10;
      എങ്കിൽ (group_lda < group_nrowa)  info =  8;
      എങ്കിൽ (group_k < 0)        info =  5;
      എങ്കിൽ (group_n < 0)        info =  4;
      എങ്കിൽ (group_m < 0)        info =  3;
      എങ്കിൽ (group_transb < 0)        info =  2;
      എങ്കിൽ (group_transa < 0)        info =  1;

    }അല്ലെങ്കിൽ എങ്കിൽ (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      എങ്കിൽ (transb == CblasNoTrans)     group_transa = 0;
      എങ്കിൽ (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      എങ്കിൽ (transb == CblasConjNoTrans) group_transa = 0;
      എങ്കിൽ (transb == CblasConjTrans)   group_transa = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (transb == CblasConjNoTrans) group_transa = 2;
      എങ്കിൽ (transb == CblasConjTrans)   group_transa = 3;
#endif
      എങ്കിൽ (transa == CblasNoTrans)     group_transb = 0;
      എങ്കിൽ (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      എങ്കിൽ (transa == CblasConjNoTrans) group_transb = 0;
      എങ്കിൽ (transa == CblasConjTrans)   group_transb = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (transa == CblasConjNoTrans) group_transb = 2;
      എങ്കിൽ (transa == CblasConjTrans)   group_transb = 3;
#endif

#അല്ലെങ്കിൽ
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      എങ്കിൽ (ta == 'N')   group_transa = 0;
      എങ്കിൽ (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      എങ്കിൽ (ta == 'C')   group_transa = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (ta == 'C')   group_transa = 3;
#endif
      എങ്കിൽ (tb == 'N')   group_transb = 0;
      എങ്കിൽ (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      എങ്കിൽ (tb == 'C')   group_transb = 1;
#അല്ലെങ്കിൽ
      എങ്കിൽ (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      എങ്കിൽ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      എങ്കിൽ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      എങ്കിൽ (group_ldc < group_m) info = 13;
      എങ്കിൽ (group_ldb < group_nrowb)  info = 10;
      എങ്കിൽ (group_lda < group_nrowa)  info =  8;
      എങ്കിൽ (group_k < 0)        info =  5;
      എങ്കിൽ (group_n < 0)        info =  4;
      എങ്കിൽ (group_m < 0)        info =  3;
      എങ്കിൽ (group_transb < 0)        info =  2;
      എങ്കിൽ (group_transa < 0)        info =  1;      
#എങ്കിൽ defined(CBLAS)
    }
#endif

    എങ്കിൽ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, വലുപ്പം(ERROR_NAME));
      സ്വതന്ത്രം(args_array);
      തിരികെ;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    എങ്കിൽ (group_m == 0 || group_n == 0) { സ്വതന്ത്രം(args_array); തിരികെ; }

    group_mode=mode;

#എങ്കിൽ defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (ഇരട്ട) group_m * (ഇരട്ട) group_n * (ഇരട്ട) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#എങ്കിൽ !defined(COMPLEX)
    എങ്കിൽ (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#അല്ലെങ്കിൽ
    എങ്കിൽ (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=ശൂന്യ;
#എങ്കിൽ !defined(COMPLEX)
      എങ്കിൽ(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }അല്ലെങ്കിൽ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#അല്ലെങ്കിൽ
      എങ്കിൽ(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }അല്ലെങ്കിൽ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }അല്ലെങ്കിൽ{
#endif
      group_routine=(ശൂന്യം*)(gemm[(group_transb<<2)|group_transa]);
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
#അല്ലെങ്കിൽ
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#എങ്കിൽ defined(CBLAS)
      എങ്കിൽ (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }അല്ലെങ്കിൽ എങ്കിൽ(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#അല്ലെങ്കിൽ
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      എങ്കിൽ (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  എങ്കിൽ(എണ്ണം>0) {
    GEMM_BATCH_THREAD(args_array,എണ്ണം);
  }

  സ്വതന്ത്രം(args_array);
}
