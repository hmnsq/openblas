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

#ଅନ୍ତର୍ଭାଵ <stdio.h>
#ଅନ୍ତର୍ଭାଵ <stdlib.h>
#ଅନ୍ତର୍ଭାଵ "common.h"

ଶୂନ୍ୟ openblas_warning(ପୂର୍ଣ୍ଣସଂଖ୍ୟା verbose, ଅଚଳ ଅକ୍ଷର * ସନ୍ଦେଶ);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#ଅନ୍ୟଥା
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#ଅନ୍ୟଥା
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#ନଚେତଯଦି defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#ଅନ୍ୟଥା
#define ERROR_NAME "CGEMM_BATCH "
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
ଶୂନ୍ୟ NAME(ଅକ୍ଷର *transa_array, ଅକ୍ଷର *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_ଗଣନା = *gcount;
#ଅନ୍ୟଥା

ଶୂନ୍ୟ CNAME(ଗଣନା CBLAS_ORDER order, ଗଣନା CBLAS_TRANSPOSE *  transa_array, ଗଣନା CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_ଗଣନା, blasint * group_size) {
#ଅନ୍ୟଥା
	   ଶୂନ୍ୟ * valpha_array,
	   ଶୂନ୍ୟ ** va_array, blasint * lda_array,
	   ଶୂନ୍ୟ ** vb_array, blasint * ldb_array,
	   ଶୂନ୍ୟ * vbeta_array,
	   ଶୂନ୍ୟ ** vc_array, blasint * ldc_array, blasint group_ଗଣନା, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ଖାଲି;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା mode=0, group_mode=0;
  blasint total_ଅଙ୍କ=0;

  blasint i=0, j=0, matrix_idx=0, ଗଣନା_2=0;

  ପୂର୍ଣ୍ଣସଂଖ୍ୟା group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  ଶୂନ୍ୟ * group_alpha, * group_beta;
  ଶୂନ୍ୟ * group_routine=ଖାଲି;
#ifdef SMALL_MATRIX_OPT
  ଶୂନ୍ୟ * group_small_matrix_opt_routine=ଖାଲି;
#endif
  
#ଯଦି defined (SMP) || defined(SMALL_MATRIX_OPT)
  ଦ୍ୱିଗୁଣ MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  ପାଇଁ(i=0; i<group_ଗଣନା; i++){
    total_ଅଙ୍କ+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_ଅଙ୍କ * ଆକାର(blas_arg_t));
  
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
  
  ପାଇଁ(i=0; i<group_ଗଣନା; matrix_idx+=group_size[i], i++){
    group_alpha = (ଶୂନ୍ୟ *)&alpha_array[i * COMPSIZE];
    group_beta = (ଶୂନ୍ୟ *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#ଯଦି defined(CBLAS) 
    ଯଦି (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      ଯଦି (transa_array[i] == CblasNoTrans)     group_transa = 0;
      ଯଦି (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ଯଦି (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      ଯଦି (transa_array[i] == CblasConjTrans)   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      ଯଦି (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      ଯଦି (transb_array[i] == CblasNoTrans)     group_transb = 0;
      ଯଦି (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ଯଦି (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      ଯଦି (transb_array[i] == CblasConjTrans)   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      ଯଦି (transb_array[i] == CblasConjTrans)   group_transb = 3;
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
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      ଯଦି (transb_array[i] == CblasNoTrans)     group_transa = 0;
      ଯଦି (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ଯଦି (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      ଯଦି (transb_array[i] == CblasConjTrans)   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      ଯଦି (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      ଯଦି (transa_array[i] == CblasNoTrans)     group_transb = 0;
      ଯଦି (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ଯଦି (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      ଯଦି (transa_array[i] == CblasConjTrans)   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      ଯଦି (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#ଅନ୍ୟଥା
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      ଯଦି (transa_array[i] == 'N')   group_transa = 0;
      ଯଦି (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      ଯଦି (transa_array[i] == 'C')   group_transa = 1;
#ଅନ୍ୟଥା
      ଯଦି (transa_array[i] == 'C')   group_transa = 3;
#endif
      ଯଦି (transb_array[i] == 'N')   group_transb = 0;
      ଯଦି (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      ଯଦି (transb_array[i] == 'C')   group_transb = 1;
#ଅନ୍ୟଥା
      ଯଦି (transb_array[i] == 'C')   group_transb = 3;
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

    ଯଦି (group_m == 0 || group_n == 0) ଜାରି;

    group_mode=mode;

#ଯଦି defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (ଦ୍ୱିଗୁଣ) group_m * (ଦ୍ୱିଗୁଣ) group_n * (ଦ୍ୱିଗୁଣ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#ଯଦି !defined(COMPLEX)
    ଯଦି (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#ଅନ୍ୟଥା
    ଯଦି (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=ଖାଲି;
#ଯଦି !defined(COMPLEX)
      ଯଦି(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ଅନ୍ୟଥା{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#ଅନ୍ୟଥା
      ଯଦି(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
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

    
    ପାଇଁ(j=0; j<group_size[i]; j++){
      args_array[ଗଣନା_2].m=group_m;
      args_array[ଗଣନା_2].n=group_n;
      args_array[ଗଣନା_2].k=group_k;
      args_array[ଗଣନା_2].lda=group_lda;
      args_array[ଗଣନା_2].ldb=group_ldb;
      args_array[ଗଣନା_2].ldc=group_ldc;
      args_array[ଗଣନା_2].alpha=group_alpha;
      args_array[ଗଣନା_2].beta=group_beta;

#ଯଦି defined(CBLAS)
      ଯଦି (order == CblasRowMajor) {
	args_array[ଗଣନା_2].a=(b_array[matrix_idx+j]);
	args_array[ଗଣନା_2].b=(a_array[matrix_idx+j]);
      } ଅନ୍ୟଥା
#endif
      {
	args_array[ଗଣନା_2].a=(a_array[matrix_idx+j]);
	args_array[ଗଣନା_2].b=(b_array[matrix_idx+j]);
      }
      
      args_array[ଗଣନା_2].c=(c_array[matrix_idx+j]);
      
      args_array[ଗଣନା_2].routine_mode=group_mode;
      args_array[ଗଣନା_2].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      ଯଦି (!group_routine)
        args_array[ଗଣନା_2].routine=group_small_matrix_opt_routine;
#endif      
      ଗଣନା_2++;
    }
  }

  ଯଦି(ଗଣନା_2>0){
    GEMM_BATCH_THREAD(args_array,ଗଣନା_2);
  }

  ମୁକ୍ତ(args_array);
}
