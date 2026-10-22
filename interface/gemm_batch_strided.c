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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ "common.h"

ಶೂನ್ಯ openblas_warning(ಪೂರ್ಣಾಙ್ಕ verbose, ಧ್ರುವ ವರ್ಣ * ಸಂದೇಶ);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
ಸ್ಥಿರ ಪೂರ್ಣಾಙ್ಕ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ಯದಿ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#ಅನ್ಯಥಾ
#define USE_SMALL_MATRIX_OPT 0
#endif

#ಯದಿ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((ಶೂನ್ಯ *)(table[idx]))
#ಅನ್ಯಥಾ
#define SMALL_KERNEL_ADDR(table, idx) ((ಶೂನ್ಯ *)(*(uintptr_t *)((ವರ್ಣ *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ಸ್ಥಿರ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ಸ್ಥಿರ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ಅನ್ಯಥಾ

ಸ್ಥಿರ size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ಸ್ಥಿರ size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಙ್ಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
ಶೂನ್ಯ NAME(ವರ್ಣ *transa, ವರ್ಣ *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           ವರ್ಣ ta = *transa;
   	   ವರ್ಣ tb = *transb;	   
	   blasint ಎಣಿಕೆ = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#ಯದಿ !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#ಅನ್ಯಥಾ
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_TRANSPOSE transa, ಗಣನಾ CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint ಎಣಿಕೆ) {
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ * valpha,
	   ಶೂನ್ಯ * va, blasint lda, blasint stridea,
	   ಶೂನ್ಯ * vb, blasint ldb, blasint strideb,
	   ಶೂನ್ಯ * vbeta,
	   ಶೂನ್ಯ * vc, blasint ldc, blasint stridec, blasint ಎಣಿಕೆ) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ರಿಕ್ತ;

  ಪೂರ್ಣಾಙ್ಕ mode=0, group_mode=0;

  blasint i=0;

  ಪೂರ್ಣಾಙ್ಕ group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  ಶೂನ್ಯ * group_routine=ರಿಕ್ತ;
#ifdef SMALL_MATRIX_OPT
  ಶೂನ್ಯ * group_small_matrix_opt_routine=ರಿಕ್ತ;
#endif
  
#ಯದಿ defined (SMP) || defined(SMALL_MATRIX_OPT)
  ದ್ವಿಗುಣ MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(ಎಣಿಕೆ * ಪರಿಮಾಣ(blas_arg_t));
  
  ಯದಿ(args_array == ರಿಕ್ತ){
    openblas_warning(0, "memory alloc failed!\n");
    ಪ್ರತ್ಯಾಗಮ;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#ಅನ್ಯಥಾ
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಅನ್ಯಥಾ
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  ಪ್ರತಿ(i=0; i<ಎಣಿಕೆ; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#ಯದಿ defined(CBLAS)
    ಯದಿ (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      ಯದಿ (transa == CblasNoTrans)     group_transa = 0;
      ಯದಿ (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (transa == CblasConjNoTrans) group_transa = 0;
      ಯದಿ (transa == CblasConjTrans)   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transa == CblasConjNoTrans) group_transa = 2;
      ಯದಿ (transa == CblasConjTrans)   group_transa = 3;
#endif
      ಯದಿ (transb == CblasNoTrans)     group_transb = 0;
      ಯದಿ (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (transb == CblasConjNoTrans) group_transb = 0;
      ಯದಿ (transb == CblasConjTrans)   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transb == CblasConjNoTrans) group_transb = 2;
      ಯದಿ (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      ಯದಿ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ಯದಿ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ಯದಿ (group_ldc < group_m) info = 13;
      ಯದಿ (group_ldb < group_nrowb)  info = 10;
      ಯದಿ (group_lda < group_nrowa)  info =  8;
      ಯದಿ (group_k < 0)        info =  5;
      ಯದಿ (group_n < 0)        info =  4;
      ಯದಿ (group_m < 0)        info =  3;
      ಯದಿ (group_transb < 0)        info =  2;
      ಯದಿ (group_transa < 0)        info =  1;

    }ಅನ್ಯಥಾ ಯದಿ (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      ಯದಿ (transb == CblasNoTrans)     group_transa = 0;
      ಯದಿ (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (transb == CblasConjNoTrans) group_transa = 0;
      ಯದಿ (transb == CblasConjTrans)   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transb == CblasConjNoTrans) group_transa = 2;
      ಯದಿ (transb == CblasConjTrans)   group_transa = 3;
#endif
      ಯದಿ (transa == CblasNoTrans)     group_transb = 0;
      ಯದಿ (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (transa == CblasConjNoTrans) group_transb = 0;
      ಯದಿ (transa == CblasConjTrans)   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transa == CblasConjNoTrans) group_transb = 2;
      ಯದಿ (transa == CblasConjTrans)   group_transb = 3;
#endif

#ಅನ್ಯಥಾ
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      ಯದಿ (ta == 'N')   group_transa = 0;
      ಯದಿ (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (ta == 'C')   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (ta == 'C')   group_transa = 3;
#endif
      ಯದಿ (tb == 'N')   group_transb = 0;
      ಯದಿ (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (tb == 'C')   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      ಯದಿ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      ಯದಿ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      ಯದಿ (group_ldc < group_m) info = 13;
      ಯದಿ (group_ldb < group_nrowb)  info = 10;
      ಯದಿ (group_lda < group_nrowa)  info =  8;
      ಯದಿ (group_k < 0)        info =  5;
      ಯದಿ (group_n < 0)        info =  4;
      ಯದಿ (group_m < 0)        info =  3;
      ಯದಿ (group_transb < 0)        info =  2;
      ಯದಿ (group_transa < 0)        info =  1;      
#ಯದಿ defined(CBLAS)
    }
#endif

    ಯದಿ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, ಪರಿಮಾಣ(ERROR_NAME));
      ಮುಕ್ತ(args_array);
      ಪ್ರತ್ಯಾಗಮ;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    ಯದಿ (group_m == 0 || group_n == 0) { ಮುಕ್ತ(args_array); ಪ್ರತ್ಯಾಗಮ; }

    group_mode=mode;

#ಯದಿ defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (ದ್ವಿಗುಣ) group_m * (ದ್ವಿಗುಣ) group_n * (ದ್ವಿಗುಣ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#ಯದಿ !defined(COMPLEX)
    ಯದಿ (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#ಅನ್ಯಥಾ
    ಯದಿ (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=ರಿಕ್ತ;
#ಯದಿ !defined(COMPLEX)
      ಯದಿ(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ಅನ್ಯಥಾ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#ಅನ್ಯಥಾ
      ಯದಿ(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ಅನ್ಯಥಾ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }ಅನ್ಯಥಾ{
#endif
      group_routine=(ಶೂನ್ಯ*)(gemm[(group_transb<<2)|group_transa]);
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
#ಅನ್ಯಥಾ
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#ಯದಿ defined(CBLAS)
      ಯದಿ (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }ಅನ್ಯಥಾ ಯದಿ(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#ಅನ್ಯಥಾ
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      ಯದಿ (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  ಯದಿ(ಎಣಿಕೆ>0) {
    GEMM_BATCH_THREAD(args_array,ಎಣಿಕೆ);
  }

  ಮುಕ್ತ(args_array);
}
