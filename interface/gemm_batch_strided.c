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

#અન્તર્ભાવ <stdio.h>
#અન્તર્ભાવ <stdlib.h>
#અન્તર્ભાવ "common.h"

શૂન્ય openblas_warning(પૂર્ણાંક verbose, અચલ અક્ષર * સંદેશ);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#નહીંતર
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#નહીંતર
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
સ્થિર પૂર્ણાંક (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#જો defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#નહીંતર
#define USE_SMALL_MATRIX_OPT 0
#endif

#જો USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((શૂન્ય *)(table[idx]))
#નહીંતર
#define SMALL_KERNEL_ADDR(table, idx) ((શૂન્ય *)(*(uintptr_t *)((અક્ષર *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
સ્થિર size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


સ્થિર size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (પૂર્ણાંક (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (પૂર્ણાંક (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#નહીંતર

સ્થિર size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

સ્થિર size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (પૂર્ણાંક (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (પૂર્ણાંક (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
શૂન્ય NAME(અક્ષર *transa, અક્ષર *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           અક્ષર ta = *transa;
   	   અક્ષર tb = *transb;	   
	   blasint ગણતરી = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#જો !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#નહીંતર
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_TRANSPOSE transa, ગણના CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint ગણતરી) {
#નહીંતર
	   શૂન્ય * valpha,
	   શૂન્ય * va, blasint lda, blasint stridea,
	   શૂન્ય * vb, blasint ldb, blasint strideb,
	   શૂન્ય * vbeta,
	   શૂન્ય * vc, blasint ldc, blasint stridec, blasint ગણતરી) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ખાલી;

  પૂર્ણાંક mode=0, group_mode=0;

  blasint i=0;

  પૂર્ણાંક group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  શૂન્ય * group_routine=ખાલી;
#ifdef SMALL_MATRIX_OPT
  શૂન્ય * group_small_matrix_opt_routine=ખાલી;
#endif
  
#જો defined (SMP) || defined(SMALL_MATRIX_OPT)
  દ્વિગુણ MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(ગણતરી * કદ(blas_arg_t));
  
  જો(args_array == ખાલી){
    openblas_warning(0, "memory alloc failed!\n");
    પાછા;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#અથવાજો defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#નહીંતર
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#નહીંતર
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#અથવાજો defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#નહીંતર
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  માટે(i=0; i<ગણતરી; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#જો defined(CBLAS)
    જો (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      જો (transa == CblasNoTrans)     group_transa = 0;
      જો (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      જો (transa == CblasConjNoTrans) group_transa = 0;
      જો (transa == CblasConjTrans)   group_transa = 1;
#નહીંતર
      જો (transa == CblasConjNoTrans) group_transa = 2;
      જો (transa == CblasConjTrans)   group_transa = 3;
#endif
      જો (transb == CblasNoTrans)     group_transb = 0;
      જો (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      જો (transb == CblasConjNoTrans) group_transb = 0;
      જો (transb == CblasConjTrans)   group_transb = 1;
#નહીંતર
      જો (transb == CblasConjNoTrans) group_transb = 2;
      જો (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      જો (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      જો (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      જો (group_ldc < group_m) info = 13;
      જો (group_ldb < group_nrowb)  info = 10;
      જો (group_lda < group_nrowa)  info =  8;
      જો (group_k < 0)        info =  5;
      જો (group_n < 0)        info =  4;
      જો (group_m < 0)        info =  3;
      જો (group_transb < 0)        info =  2;
      જો (group_transa < 0)        info =  1;

    }નહીંતર જો (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      જો (transb == CblasNoTrans)     group_transa = 0;
      જો (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      જો (transb == CblasConjNoTrans) group_transa = 0;
      જો (transb == CblasConjTrans)   group_transa = 1;
#નહીંતર
      જો (transb == CblasConjNoTrans) group_transa = 2;
      જો (transb == CblasConjTrans)   group_transa = 3;
#endif
      જો (transa == CblasNoTrans)     group_transb = 0;
      જો (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      જો (transa == CblasConjNoTrans) group_transb = 0;
      જો (transa == CblasConjTrans)   group_transb = 1;
#નહીંતર
      જો (transa == CblasConjNoTrans) group_transb = 2;
      જો (transa == CblasConjTrans)   group_transb = 3;
#endif

#નહીંતર
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      જો (ta == 'N')   group_transa = 0;
      જો (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      જો (ta == 'C')   group_transa = 1;
#નહીંતર
      જો (ta == 'C')   group_transa = 3;
#endif
      જો (tb == 'N')   group_transb = 0;
      જો (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      જો (tb == 'C')   group_transb = 1;
#નહીંતર
      જો (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      જો (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      જો (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      જો (group_ldc < group_m) info = 13;
      જો (group_ldb < group_nrowb)  info = 10;
      જો (group_lda < group_nrowa)  info =  8;
      જો (group_k < 0)        info =  5;
      જો (group_n < 0)        info =  4;
      જો (group_m < 0)        info =  3;
      જો (group_transb < 0)        info =  2;
      જો (group_transa < 0)        info =  1;      
#જો defined(CBLAS)
    }
#endif

    જો (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, કદ(ERROR_NAME));
      મુક્ત(args_array);
      પાછા;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    જો (group_m == 0 || group_n == 0) { મુક્ત(args_array); પાછા; }

    group_mode=mode;

#જો defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (દ્વિગુણ) group_m * (દ્વિગુણ) group_n * (દ્વિગુણ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#જો !defined(COMPLEX)
    જો (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#નહીંતર
    જો (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=ખાલી;
#જો !defined(COMPLEX)
      જો(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }નહીંતર{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#નહીંતર
      જો(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }નહીંતર{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }નહીંતર{
#endif
      group_routine=(શૂન્ય*)(gemm[(group_transb<<2)|group_transa]);
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
#નહીંતર
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#જો defined(CBLAS)
      જો (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }નહીંતર જો(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#નહીંતર
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      જો (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  જો(ગણતરી>0) {
    GEMM_BATCH_THREAD(args_array,ગણતરી);
  }

  મુક્ત(args_array);
}
