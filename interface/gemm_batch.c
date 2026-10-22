/*****************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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
#define ERROR_NAME "QGEMM_BATCH "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#નહીંતર
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#નહીંતર
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#અથવાજો defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#નહીંતર
#define ERROR_NAME "CGEMM_BATCH "
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
શૂન્ય NAME(અક્ષર *transa_array, અક્ષર *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_ગણતરી = *gcount;
#નહીંતર

શૂન્ય CNAME(ગણના CBLAS_ORDER order, ગણના CBLAS_TRANSPOSE *  transa_array, ગણના CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_ગણતરી, blasint * group_size) {
#નહીંતર
	   શૂન્ય * valpha_array,
	   શૂન્ય ** va_array, blasint * lda_array,
	   શૂન્ય ** vb_array, blasint * ldb_array,
	   શૂન્ય * vbeta_array,
	   શૂન્ય ** vc_array, blasint * ldc_array, blasint group_ગણતરી, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ખાલી;

  પૂર્ણાંક mode=0, group_mode=0;
  blasint total_અંક=0;

  blasint i=0, j=0, matrix_idx=0, ગણતરી=0;

  પૂર્ણાંક group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  શૂન્ય * group_alpha, * group_beta;
  શૂન્ય * group_routine=ખાલી;
#ifdef SMALL_MATRIX_OPT
  શૂન્ય * group_small_matrix_opt_routine=ખાલી;
#endif
  
#જો defined (SMP) || defined(SMALL_MATRIX_OPT)
  દ્વિગુણ MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  માટે(i=0; i<group_ગણતરી; i++){
    total_અંક+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_અંક * કદ(blas_arg_t));
  
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
  
  માટે(i=0; i<group_ગણતરી; matrix_idx+=group_size[i], i++){
    group_alpha = (શૂન્ય *)&alpha_array[i * COMPSIZE];
    group_beta = (શૂન્ય *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#જો defined(CBLAS) 
    જો (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      જો (transa_array[i] == CblasNoTrans)     group_transa = 0;
      જો (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      જો (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      જો (transa_array[i] == CblasConjTrans)   group_transa = 1;
#નહીંતર
      જો (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      જો (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      જો (transb_array[i] == CblasNoTrans)     group_transb = 0;
      જો (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      જો (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      જો (transb_array[i] == CblasConjTrans)   group_transb = 1;
#નહીંતર
      જો (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      જો (transb_array[i] == CblasConjTrans)   group_transb = 3;
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
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      જો (transb_array[i] == CblasNoTrans)     group_transa = 0;
      જો (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      જો (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      જો (transb_array[i] == CblasConjTrans)   group_transa = 1;
#નહીંતર
      જો (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      જો (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      જો (transa_array[i] == CblasNoTrans)     group_transb = 0;
      જો (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      જો (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      જો (transa_array[i] == CblasConjTrans)   group_transb = 1;
#નહીંતર
      જો (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      જો (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#નહીંતર
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      જો (transa_array[i] == 'N')   group_transa = 0;
      જો (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      જો (transa_array[i] == 'C')   group_transa = 1;
#નહીંતર
      જો (transa_array[i] == 'C')   group_transa = 3;
#endif
      જો (transb_array[i] == 'N')   group_transb = 0;
      જો (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      જો (transb_array[i] == 'C')   group_transb = 1;
#નહીંતર
      જો (transb_array[i] == 'C')   group_transb = 3;
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

    જો (group_m == 0 || group_n == 0) ચાલુ;

    group_mode=mode;

#જો defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (દ્વિગુણ) group_m * (દ્વિગુણ) group_n * (દ્વિગુણ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#જો !defined(COMPLEX)
    જો (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#નહીંતર
    જો (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=ખાલી;
#જો !defined(COMPLEX)
      જો(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }નહીંતર{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#નહીંતર
      જો(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
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

    
    માટે(j=0; j<group_size[i]; j++){
      args_array[ગણતરી].m=group_m;
      args_array[ગણતરી].n=group_n;
      args_array[ગણતરી].k=group_k;
      args_array[ગણતરી].lda=group_lda;
      args_array[ગણતરી].ldb=group_ldb;
      args_array[ગણતરી].ldc=group_ldc;
      args_array[ગણતરી].alpha=group_alpha;
      args_array[ગણતરી].beta=group_beta;

#જો defined(CBLAS)
      જો (order == CblasRowMajor) {
	args_array[ગણતરી].a=(b_array[matrix_idx+j]);
	args_array[ગણતરી].b=(a_array[matrix_idx+j]);
      } નહીંતર
#endif
      {
	args_array[ગણતરી].a=(a_array[matrix_idx+j]);
	args_array[ગણતરી].b=(b_array[matrix_idx+j]);
      }
      
      args_array[ગણતરી].c=(c_array[matrix_idx+j]);
      
      args_array[ગણતરી].routine_mode=group_mode;
      args_array[ગણતરી].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      જો (!group_routine)
        args_array[ગણતરી].routine=group_small_matrix_opt_routine;
#endif      
      ગણતરી++;
    }
  }

  જો(ગણતરી>0){
    GEMM_BATCH_THREAD(args_array,ગણતરી);
  }

  મુક્ત(args_array);
}
