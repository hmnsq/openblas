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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ "common.h"

ଶୂନ୍ୟ openblas_warning(ପୂର୍ଣ୍ଣସଂଖ୍ୟା verbose, ଅଚଳ ଅକ୍ଷର * ସନ୍ଦେଶ);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
ସ୍ଥିର ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ଯଦି defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#ଅନ୍ୟଥା
#define USE_SMALL_MATRIX_OPT 0
#endif

#ଯଦି USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ଶୂନ୍ୟ *)(table[idx]))
#ଅନ୍ୟଥା
#define SMALL_KERNEL_ADDR(table, idx) ((ଶୂନ୍ୟ *)(*(uintptr_t *)((ଅକ୍ଷର *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ସ୍ଥିର size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ସ୍ଥିର size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ଅନ୍ୟଥା

ସ୍ଥିର size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ସ୍ଥିର size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (ପୂର୍ଣ୍ଣସଂଖ୍ୟା (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
ଶୂନ୍ୟ NAME(ଅକ୍ଷର *transa, ଅକ୍ଷର *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           ଅକ୍ଷର ta = *transa;
   	   ଅକ୍ଷର tb = *transb;	   
	   blasint ଗଣନା_2 = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#ଯଦି !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#ଅନ୍ୟଥା
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_TRANSPOSE transa, ଗଣନା CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint ଗଣନା_2) {
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ * valpha,
	   ଶୂନ୍ୟ * va, blasint lda, blasint stridea,
	   ଶୂନ୍ୟ * vb, blasint ldb, blasint strideb,
	   ଶୂନ୍ୟ * vbeta,
	   ଶୂନ୍ୟ * vc, blasint ldc, blasint stridec, blasint ଗଣନା_2) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ଖାଲି;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode=0, group_mode=0;

  blasint i=0;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  ଶୂନ୍ୟ * group_routine=ଖାଲି;
#ifdef SMALL_MATRIX_OPT
  ଶୂନ୍ୟ * group_small_matrix_opt_routine=ଖାଲି;
#endif
  
#ଯଦି defined (SMP) || defined(SMALL_MATRIX_OPT)
  ଦ୍ୱିଗୁଣ MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(ଗଣନା_2 * ଆକାର(blas_arg_t));
  
  ଯଦି(args_array == ଖାଲି){
    openblas_warning(0, "memory alloc failed!\n");
    ଫେରନ୍ତୁ;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ନଚେତଯଦି defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ଅନ୍ୟଥା
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ନଚେତଯଦି defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ଅନ୍ୟଥା
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  ପାଇଁ(i=0; i<ଗଣନା_2; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#ଯଦି defined(CBLAS)
    ଯଦି (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      ଯଦି (transa == CblasNoTrans)     group_transa = 0;
      ଯଦି (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ଯଦି (transa == CblasConjNoTrans) group_transa = 0;
      ଯଦି (transa == CblasConjTrans)   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (transa == CblasConjNoTrans) group_transa = 2;
      ଯଦି (transa == CblasConjTrans)   group_transa = 3;
#endif
      ଯଦି (transb == CblasNoTrans)     group_transb = 0;
      ଯଦି (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ଯଦି (transb == CblasConjNoTrans) group_transb = 0;
      ଯଦି (transb == CblasConjTrans)   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (transb == CblasConjNoTrans) group_transb = 2;
      ଯଦି (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      ଯଦି (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ଯଦି (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ଯଦି (group_ldc < group_m) info = 13;
      ଯଦି (group_ldb < group_nrowb)  info = 10;
      ଯଦି (group_lda < group_nrowa)  info =  8;
      ଯଦି (group_k < 0)        info =  5;
      ଯଦି (group_n < 0)        info =  4;
      ଯଦି (group_m < 0)        info =  3;
      ଯଦି (group_transb < 0)        info =  2;
      ଯଦି (group_transa < 0)        info =  1;

    }ଅନ୍ୟଥା ଯଦି (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      ଯଦି (transb == CblasNoTrans)     group_transa = 0;
      ଯଦି (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ଯଦି (transb == CblasConjNoTrans) group_transa = 0;
      ଯଦି (transb == CblasConjTrans)   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (transb == CblasConjNoTrans) group_transa = 2;
      ଯଦି (transb == CblasConjTrans)   group_transa = 3;
#endif
      ଯଦି (transa == CblasNoTrans)     group_transb = 0;
      ଯଦି (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ଯଦି (transa == CblasConjNoTrans) group_transb = 0;
      ଯଦି (transa == CblasConjTrans)   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (transa == CblasConjNoTrans) group_transb = 2;
      ଯଦି (transa == CblasConjTrans)   group_transb = 3;
#endif

#ଅନ୍ୟଥା
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      ଯଦି (ta == 'N')   group_transa = 0;
      ଯଦି (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      ଯଦି (ta == 'C')   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (ta == 'C')   group_transa = 3;
#endif
      ଯଦି (tb == 'N')   group_transb = 0;
      ଯଦି (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      ଯଦି (tb == 'C')   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      ଯଦି (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ଯଦି (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ଯଦି (group_ldc < group_m) info = 13;
      ଯଦି (group_ldb < group_nrowb)  info = 10;
      ଯଦି (group_lda < group_nrowa)  info =  8;
      ଯଦି (group_k < 0)        info =  5;
      ଯଦି (group_n < 0)        info =  4;
      ଯଦି (group_m < 0)        info =  3;
      ଯଦି (group_transb < 0)        info =  2;
      ଯଦି (group_transa < 0)        info =  1;      
#ଯଦି defined(CBLAS)
    }
#endif

    ଯଦି (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, ଆକାର(ERROR_NAME));
      ମୁକ୍ତ(args_array);
      ଫେରନ୍ତୁ;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    ଯଦି (group_m == 0 || group_n == 0) { ମୁକ୍ତ(args_array); ଫେରନ୍ତୁ; }

    group_mode=mode;

#ଯଦି defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (ଦ୍ୱିଗୁଣ) group_m * (ଦ୍ୱିଗୁଣ) group_n * (ଦ୍ୱିଗୁଣ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#ଯଦି !defined(COMPLEX)
    ଯଦି (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#ଅନ୍ୟଥା
    ଯଦି (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=ଖାଲି;
#ଯଦି !defined(COMPLEX)
      ଯଦି(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ଅନ୍ୟଥା{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#ଅନ୍ୟଥା
      ଯଦି(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ଅନ୍ୟଥା{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }ଅନ୍ୟଥା{
#endif
      group_routine=(ଶୂନ୍ୟ*)(gemm[(group_transb<<2)|group_transa]);
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
#ଅନ୍ୟଥା
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#ଯଦି defined(CBLAS)
      ଯଦି (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }ଅନ୍ୟଥା ଯଦି(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#ଅନ୍ୟଥା
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      ଯଦି (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  ଯଦି(ଗଣନା_2>0) {
    GEMM_BATCH_THREAD(args_array,ଗଣନା_2);
  }

  ମୁକ୍ତ(args_array);
}
