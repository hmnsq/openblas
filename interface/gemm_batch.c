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

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ "common.h"

ಶೂನ್ಯ openblas_warning(ಪೂರ್ಣಾಙ್ಕ verbose, ಧ್ರುವ ವರ್ಣ * ಸಂದೇಶ);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#ಅನ್ಯಥಾ
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#ಅನ್ಯಥಾ
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#ನೋಚೇದ್ಯದಿ defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#ಅನ್ಯಥಾ
#define ERROR_NAME "CGEMM_BATCH "
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
ಶೂನ್ಯ NAME(ವರ್ಣ *transa_array, ವರ್ಣ *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_ಎಣಿಕೆ = *gcount;
#ಅನ್ಯಥಾ

ಶೂನ್ಯ CNAME(ಗಣನಾ CBLAS_ORDER order, ಗಣನಾ CBLAS_TRANSPOSE *  transa_array, ಗಣನಾ CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_ಎಣಿಕೆ, blasint * group_size) {
#ಅನ್ಯಥಾ
	   ಶೂನ್ಯ * valpha_array,
	   ಶೂನ್ಯ ** va_array, blasint * lda_array,
	   ಶೂನ್ಯ ** vb_array, blasint * ldb_array,
	   ಶೂನ್ಯ * vbeta_array,
	   ಶೂನ್ಯ ** vc_array, blasint * ldc_array, blasint group_ಎಣಿಕೆ, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=ರಿಕ್ತ;

  ಪೂರ್ಣಾಙ್ಕ mode=0, group_mode=0;
  blasint total_ಅಂಕೆ=0;

  blasint i=0, j=0, matrix_idx=0, ಎಣಿಕೆ=0;

  ಪೂರ್ಣಾಙ್ಕ group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  ಶೂನ್ಯ * group_alpha, * group_beta;
  ಶೂನ್ಯ * group_routine=ರಿಕ್ತ;
#ifdef SMALL_MATRIX_OPT
  ಶೂನ್ಯ * group_small_matrix_opt_routine=ರಿಕ್ತ;
#endif
  
#ಯದಿ defined (SMP) || defined(SMALL_MATRIX_OPT)
  ದ್ವಿಗುಣ MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  ಪ್ರತಿ(i=0; i<group_ಎಣಿಕೆ; i++){
    total_ಅಂಕೆ+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_ಅಂಕೆ * ಪರಿಮಾಣ(blas_arg_t));
  
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
  
  ಪ್ರತಿ(i=0; i<group_ಎಣಿಕೆ; matrix_idx+=group_size[i], i++){
    group_alpha = (ಶೂನ್ಯ *)&alpha_array[i * COMPSIZE];
    group_beta = (ಶೂನ್ಯ *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#ಯದಿ defined(CBLAS) 
    ಯದಿ (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      ಯದಿ (transa_array[i] == CblasNoTrans)     group_transa = 0;
      ಯದಿ (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      ಯದಿ (transa_array[i] == CblasConjTrans)   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      ಯದಿ (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      ಯದಿ (transb_array[i] == CblasNoTrans)     group_transb = 0;
      ಯದಿ (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      ಯದಿ (transb_array[i] == CblasConjTrans)   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      ಯದಿ (transb_array[i] == CblasConjTrans)   group_transb = 3;
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
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      ಯದಿ (transb_array[i] == CblasNoTrans)     group_transa = 0;
      ಯದಿ (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      ಯದಿ (transb_array[i] == CblasConjTrans)   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      ಯದಿ (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      ಯದಿ (transa_array[i] == CblasNoTrans)     group_transb = 0;
      ಯದಿ (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      ಯದಿ (transa_array[i] == CblasConjTrans)   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      ಯದಿ (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#ಅನ್ಯಥಾ
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      ಯದಿ (transa_array[i] == 'N')   group_transa = 0;
      ಯದಿ (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      ಯದಿ (transa_array[i] == 'C')   group_transa = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transa_array[i] == 'C')   group_transa = 3;
#endif
      ಯದಿ (transb_array[i] == 'N')   group_transb = 0;
      ಯದಿ (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      ಯದಿ (transb_array[i] == 'C')   group_transb = 1;
#ಅನ್ಯಥಾ
      ಯದಿ (transb_array[i] == 'C')   group_transb = 3;
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

    ಯದಿ (group_m == 0 || group_n == 0) ಅನುವರ್ತ;

    group_mode=mode;

#ಯದಿ defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (ದ್ವಿಗುಣ) group_m * (ದ್ವಿಗುಣ) group_n * (ದ್ವಿಗುಣ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#ಯದಿ !defined(COMPLEX)
    ಯದಿ (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#ಅನ್ಯಥಾ
    ಯದಿ (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=ರಿಕ್ತ;
#ಯದಿ !defined(COMPLEX)
      ಯದಿ(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }ಅನ್ಯಥಾ{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#ಅನ್ಯಥಾ
      ಯದಿ(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
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

    
    ಪ್ರತಿ(j=0; j<group_size[i]; j++){
      args_array[ಎಣಿಕೆ].m=group_m;
      args_array[ಎಣಿಕೆ].n=group_n;
      args_array[ಎಣಿಕೆ].k=group_k;
      args_array[ಎಣಿಕೆ].lda=group_lda;
      args_array[ಎಣಿಕೆ].ldb=group_ldb;
      args_array[ಎಣಿಕೆ].ldc=group_ldc;
      args_array[ಎಣಿಕೆ].alpha=group_alpha;
      args_array[ಎಣಿಕೆ].beta=group_beta;

#ಯದಿ defined(CBLAS)
      ಯದಿ (order == CblasRowMajor) {
	args_array[ಎಣಿಕೆ].a=(b_array[matrix_idx+j]);
	args_array[ಎಣಿಕೆ].b=(a_array[matrix_idx+j]);
      } ಅನ್ಯಥಾ
#endif
      {
	args_array[ಎಣಿಕೆ].a=(a_array[matrix_idx+j]);
	args_array[ಎಣಿಕೆ].b=(b_array[matrix_idx+j]);
      }
      
      args_array[ಎಣಿಕೆ].c=(c_array[matrix_idx+j]);
      
      args_array[ಎಣಿಕೆ].routine_mode=group_mode;
      args_array[ಎಣಿಕೆ].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      ಯದಿ (!group_routine)
        args_array[ಎಣಿಕೆ].routine=group_small_matrix_opt_routine;
#endif      
      ಎಣಿಕೆ++;
    }
  }

  ಯದಿ(ಎಣಿಕೆ>0){
    GEMM_BATCH_THREAD(args_array,ಎಣಿಕೆ);
  }

  ಮುಕ್ತ(args_array);
}
