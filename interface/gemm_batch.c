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

ಶೂನ್ಯ openblas_warning(ಪೂರ್ಣಾಂಕ verbose, ಸ್ಥಿರಾಂಕ ಅಕ್ಷರ * ಸಂದೇಶ);

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "QGEMM_BATCH "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "DGEMM_BATCH "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_BATCH_THREAD dgemm_batch_thread
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "SGEMM_BATCH "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_BATCH_THREAD sgemm_batch_thread
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "XGEMM_BATCH "
#ಅಥವಾಆದರೆ defined(DOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "ZGEMM_BATCH "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_BATCH_THREAD zgemm_batch_thread
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು ERROR_NAME "CGEMM_BATCH "
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_BATCH_THREAD cgemm_batch_thread
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
ಸ್ಥಿರ ಪೂರ್ಣಾಂಕ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ಆದರೆ defined(SMALL_ಮಾತೃಕೆ_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#ವ್ಯಾಖ್ಯಾನಿಸು USE_SMALL_ಮಾತೃಕೆ_OPT 1
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು USE_SMALL_ಮಾತೃಕೆ_OPT 0
#ಆದರೆ_ಅಂತ್ಯ

#ಆದರೆ USE_SMALL_ಮಾತೃಕೆ_OPT
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ DYNAMIC_ARCH
#ವ್ಯಾಖ್ಯಾನಿಸು SMALL_KERNEL_ADDR(ಕೋಷ್ಟಕ, idx) ((ಶೂನ್ಯ *)(ಕೋಷ್ಟಕ[idx]))
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿಸು SMALL_KERNEL_ADDR(ಕೋಷ್ಟಕ, idx) ((ಶೂನ್ಯ *)(*(uintptr_t *)((ಅಕ್ಷರ *)gotoblas + (size_t)(ಕೋಷ್ಟಕ[idx]))))
#ಆದರೆ_ಅಂತ್ಯ


#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
ಸ್ಥಿರ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ಸ್ಥಿರ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#ವ್ಯಾಖ್ಯಾನಿಸು GEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ಇಲ್ಲದಿದ್ದರೆ

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

#ವ್ಯಾಖ್ಯಾನಿಸು ZGEMM_SMALL_KERNEL(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#ವ್ಯಾಖ್ಯಾನಿಸು ZGEMM_SMALL_KERNEL_B0(idx) (ಪೂರ್ಣಾಂಕ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

OPENBLAS_EXPORT

#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ CBLAS
ಶೂನ್ಯ NAME(ಅಕ್ಷರ *transa_array, ಅಕ್ಷರ *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * ಗುಂಪು_ಗಾತ್ರ) {
       blasint ಗುಂಪು_ಎಣಿಕೆ = *gcount;
#ಇಲ್ಲದಿದ್ದರೆ

ಶೂನ್ಯ CNAME(ಎಣಿಕೆ CBLAS_ಕ್ರಮ ಕ್ರಮ, ಎಣಿಕೆ CBLAS_TRANSPOSE *  transa_array, ಎಣಿಕೆ CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint ಗುಂಪು_ಎಣಿಕೆ, blasint * ಗುಂಪು_ಗಾತ್ರ) {
#ಇಲ್ಲದಿದ್ದರೆ
	   ಶೂನ್ಯ * valpha_array,
	   ಶೂನ್ಯ ** va_array, blasint * lda_array,
	   ಶೂನ್ಯ ** vb_array, blasint * ldb_array,
	   ಶೂನ್ಯ * vbeta_array,
	   ಶೂನ್ಯ ** vc_array, blasint * ldc_array, blasint ಗುಂಪು_ಎಣಿಕೆ, blasint * ಗುಂಪು_ಗಾತ್ರ) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
  BLASLONG ಗುಂಪು_m, ಗುಂಪು_n, ಗುಂಪು_k;
  BLASLONG ಗುಂಪು_lda, ಗುಂಪು_ldb, ಗುಂಪು_ldc;

  blas_arg_t * ವಾದಗಳು_array=ಖಾಲಿ;

  ಪೂರ್ಣಾಂಕ ಬಹುಲಕ=0, ಗುಂಪು_ಬಹುಲಕ=0;
  blasint ಒಟ್ಟು_ಅಂಕೆ=0;

  blasint i=0, j=0, ಮಾತೃಕೆ_idx=0, ಎಣಿಕೆ_2=0;

  ಪೂರ್ಣಾಂಕ ಗುಂಪು_transa, ಗುಂಪು_transb;
  BLASLONG ಗುಂಪು_nrowa, ಗುಂಪು_nrowb;
  blasint info;

  ಶೂನ್ಯ * ಗುಂಪು_alpha, * ಗುಂಪು_beta;
  ಶೂನ್ಯ * ಗುಂಪು_routine=ಖಾಲಿ;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMALL_ಮಾತೃಕೆ_OPT
  ಶೂನ್ಯ * ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine=ಖಾಲಿ;
#ಆದರೆ_ಅಂತ್ಯ
  
#ಆದರೆ defined (SMP) || defined(SMALL_ಮಾತೃಕೆ_OPT)
  ದ್ವಿಗುಣ MNK;
#ಆದರೆ_ಅಂತ್ಯ

  PRINT_DEBUG_CNAME;
    
  ಗಾಗಿ(i=0; i<ಗುಂಪು_ಎಣಿಕೆ; i++){
    ಒಟ್ಟು_ಅಂಕೆ+=ಗುಂಪು_ಗಾತ್ರ[i];
  }

  ವಾದಗಳು_array=(blas_arg_t *)ಸ್ಮರಣೆ_ಕೊಡು(ಒಟ್ಟು_ಅಂಕೆ * ಗಾತ್ರ(blas_arg_t));
  
  ಆದರೆ(ವಾದಗಳು_array == ಖಾಲಿ){
    openblas_warning(0, "memory alloc failed!\n");
    ಹಿಂತಿರುಗಿ;
  }

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMP
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ಇಲ್ಲದಿದ್ದರೆ
  ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_REAL;
#ಆದರೆ_ಅಂತ್ಯ
#ಇಲ್ಲದಿದ್ದರೆ
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ XDOUBLE
  ಬಹುಲಕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ಅಥವಾಆದರೆ defined(DOUBLE)
  ಬಹುಲಕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ಇಲ್ಲದಿದ್ದರೆ
  ಬಹುಲಕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ
  
  ಗಾಗಿ(i=0; i<ಗುಂಪು_ಎಣಿಕೆ; ಮಾತೃಕೆ_idx+=ಗುಂಪು_ಗಾತ್ರ[i], i++){
    ಗುಂಪು_alpha = (ಶೂನ್ಯ *)&alpha_array[i * COMPSIZE];
    ಗುಂಪು_beta = (ಶೂನ್ಯ *)&beta_array[i * COMPSIZE];

    ಗುಂಪು_m = ಗುಂಪು_n = ಗುಂಪು_k = 0;
    ಗುಂಪು_lda = ಗುಂಪು_ldb = ಗುಂಪು_ldc = 0;
    ಗುಂಪು_transa = -1;
    ಗುಂಪು_transb = -1;
    info   = 0;
   
#ಆದರೆ defined(CBLAS) 
    ಆದರೆ (ಕ್ರಮ == Cblas_ಲಂಬಸಾಲು_Major) {
      ಗುಂಪು_m = m_array[i];
      ಗುಂಪು_n = n_array[i];
      ಗುಂಪು_k = k_array[i];

      ಗುಂಪು_lda = lda_array[i];
      ಗುಂಪು_ldb = ldb_array[i];
      ಗುಂಪು_ldc = ldc_array[i];
      
      ಆದರೆ (transa_array[i] == CblasNoTrans)     ಗುಂಪು_transa = 0;
      ಆದರೆ (transa_array[i] == CblasTrans)       ಗುಂಪು_transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transa_array[i] == CblasConjNoTrans) ಗುಂಪು_transa = 0;
      ಆದರೆ (transa_array[i] == CblasConjTrans)   ಗುಂಪು_transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transa_array[i] == CblasConjNoTrans) ಗುಂಪು_transa = 2;
      ಆದರೆ (transa_array[i] == CblasConjTrans)   ಗುಂಪು_transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
      ಆದರೆ (transb_array[i] == CblasNoTrans)     ಗುಂಪು_transb = 0;
      ಆದರೆ (transb_array[i] == CblasTrans)       ಗುಂಪು_transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transb_array[i] == CblasConjNoTrans) ಗುಂಪು_transb = 0;
      ಆದರೆ (transb_array[i] == CblasConjTrans)   ಗುಂಪು_transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transb_array[i] == CblasConjNoTrans) ಗುಂಪು_transb = 2;
      ಆದರೆ (transb_array[i] == CblasConjTrans)   ಗುಂಪು_transb = 3;
#ಆದರೆ_ಅಂತ್ಯ
      ಗುಂಪು_nrowa = ಗುಂಪು_m;
      ಆದರೆ (ಗುಂಪು_transa & 1) ಗುಂಪು_nrowa = ಗುಂಪು_k;
      ಗುಂಪು_nrowb = ಗುಂಪು_k;
      ಆದರೆ (ಗುಂಪು_transb & 1) ಗುಂಪು_nrowb = ಗುಂಪು_n;

      info=-1;
	
      ಆದರೆ (ಗುಂಪು_ldc < ಗುಂಪು_m) info = 13;
      ಆದರೆ (ಗುಂಪು_ldb < ಗುಂಪು_nrowb)  info = 10;
      ಆದರೆ (ಗುಂಪು_lda < ಗುಂಪು_nrowa)  info =  8;
      ಆದರೆ (ಗುಂಪು_k < 0)        info =  5;
      ಆದರೆ (ಗುಂಪು_n < 0)        info =  4;
      ಆದರೆ (ಗುಂಪು_m < 0)        info =  3;
      ಆದರೆ (ಗುಂಪು_transb < 0)        info =  2;
      ಆದರೆ (ಗುಂಪು_transa < 0)        info =  1;

    }ಇಲ್ಲದಿದ್ದರೆ ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
      
      ಗುಂಪು_m = n_array[i];
      ಗುಂಪು_n = m_array[i];
      ಗುಂಪು_k = k_array[i];

      ಗುಂಪು_lda = ldb_array[i];
      ಗುಂಪು_ldb = lda_array[i];
      ಗುಂಪು_ldc = ldc_array[i];

      ಆದರೆ (transb_array[i] == CblasNoTrans)     ಗುಂಪು_transa = 0;
      ಆದರೆ (transb_array[i] == CblasTrans)       ಗುಂಪು_transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transb_array[i] == CblasConjNoTrans) ಗುಂಪು_transa = 0;
      ಆದರೆ (transb_array[i] == CblasConjTrans)   ಗುಂಪು_transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transb_array[i] == CblasConjNoTrans) ಗುಂಪು_transa = 2;
      ಆದರೆ (transb_array[i] == CblasConjTrans)   ಗುಂಪು_transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
      ಆದರೆ (transa_array[i] == CblasNoTrans)     ಗುಂಪು_transb = 0;
      ಆದರೆ (transa_array[i] == CblasTrans)       ಗುಂಪು_transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transa_array[i] == CblasConjNoTrans) ಗುಂಪು_transb = 0;
      ಆದರೆ (transa_array[i] == CblasConjTrans)   ಗುಂಪು_transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transa_array[i] == CblasConjNoTrans) ಗುಂಪು_transb = 2;
      ಆದರೆ (transa_array[i] == CblasConjTrans)   ಗುಂಪು_transb = 3;
#ಆದರೆ_ಅಂತ್ಯ

#ಇಲ್ಲದಿದ್ದರೆ
      ಗುಂಪು_m = m_array[i];
      ಗುಂಪು_n = n_array[i];
      ಗುಂಪು_k = k_array[i];

      ಗುಂಪು_lda = lda_array[i];
      ಗುಂಪು_ldb = ldb_array[i];
      ಗುಂಪು_ldc = ldc_array[i];
      
      ಆದರೆ (transa_array[i] == 'N')   ಗುಂಪು_transa = 0;
      ಆದರೆ (transa_array[i] == 'T')   ಗುಂಪು_transa = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transa_array[i] == 'C')   ಗುಂಪು_transa = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transa_array[i] == 'C')   ಗುಂಪು_transa = 3;
#ಆದರೆ_ಅಂತ್ಯ
      ಆದರೆ (transb_array[i] == 'N')   ಗುಂಪು_transb = 0;
      ಆದರೆ (transb_array[i] == 'T')   ಗುಂಪು_transb = 1;
#ವ್ಯಾಖ್ಯಾನಿತವಾಗದಿದ್ದರೆ COMPLEX
      ಆದರೆ (transb_array[i] == 'C')   ಗುಂಪು_transb = 1;
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ (transb_array[i] == 'C')   ಗುಂಪು_transb = 3;
#ಆದರೆ_ಅಂತ್ಯ
#ಆದರೆ_ಅಂತ್ಯ

      ಗುಂಪು_nrowa = ಗುಂಪು_m;
      ಆದರೆ (ಗುಂಪು_transa & 1) ಗುಂಪು_nrowa = ಗುಂಪು_k;
      ಗುಂಪು_nrowb = ಗುಂಪು_k;
      ಆದರೆ (ಗುಂಪು_transb & 1) ಗುಂಪು_nrowb = ಗುಂಪು_n;

      info=-1;
	
      ಆದರೆ (ಗುಂಪು_ldc < ಗುಂಪು_m) info = 13;
      ಆದರೆ (ಗುಂಪು_ldb < ಗುಂಪು_nrowb)  info = 10;
      ಆದರೆ (ಗುಂಪು_lda < ಗುಂಪು_nrowa)  info =  8;
      ಆದರೆ (ಗುಂಪು_k < 0)        info =  5;
      ಆದರೆ (ಗುಂಪು_n < 0)        info =  4;
      ಆದರೆ (ಗುಂಪು_m < 0)        info =  3;
      ಆದರೆ (ಗುಂಪು_transb < 0)        info =  2;
      ಆದರೆ (ಗುಂಪು_transa < 0)        info =  1;      
#ಆದರೆ defined(CBLAS)
    }
#ಆದರೆ_ಅಂತ್ಯ

    ಆದರೆ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, ಗಾತ್ರ(ERROR_NAME));
      ಬಿಡುಗಡೆ(ವಾದಗಳು_array);
      ಹಿಂತಿರುಗಿ;
    }

    ಆದರೆ (ಗುಂಪು_m == 0 || ಗುಂಪು_n == 0) ಮುಂದುವರಿ;

    ಗುಂಪು_ಬಹುಲಕ=ಬಹುಲಕ;

#ಆದರೆ defined(SMP) || defined(SMALL_ಮಾತೃಕೆ_OPT)
    MNK = (ದ್ವಿಗುಣ) ಗುಂಪು_m * (ದ್ವಿಗುಣ) ಗುಂಪು_n * (ದ್ವಿಗುಣ) ಗುಂಪು_k;
#ಆದರೆ_ಅಂತ್ಯ

#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMALL_ಮಾತೃಕೆ_OPT
#ಆದರೆ !defined(COMPLEX)
    ಆದರೆ (*(FLOAT *)(ಗುಂಪು_alpha) != ZERO &&
        GEMM_SMALL_ಮಾತೃಕೆ_PERMIT(ಗುಂಪು_transa, ಗುಂಪು_transb, ಗುಂಪು_m, ಗುಂಪು_n, ಗುಂಪು_k,
                                 *(FLOAT *)(ಗುಂಪು_alpha), *(FLOAT *)(ಗುಂಪು_beta))){
#ಇಲ್ಲದಿದ್ದರೆ
    ಆದರೆ (((FLOAT *)(ಗುಂಪು_alpha))[0] != ZERO &&
        GEMM_SMALL_ಮಾತೃಕೆ_PERMIT(ಗುಂಪು_transa, ಗುಂಪು_transb, ಗುಂಪು_m, ಗುಂಪು_n, ಗುಂಪು_k,
                                 ((FLOAT *)(ಗುಂಪು_alpha))[0], ((FLOAT *)(ಗುಂಪು_alpha))[1],
                                 ((FLOAT *)(ಗುಂಪು_beta))[0], ((FLOAT *)(ಗುಂಪು_beta))[1])){
#ಆದರೆ_ಅಂತ್ಯ
      ಗುಂಪು_routine=ಖಾಲಿ;
#ಆದರೆ !defined(COMPLEX)
      ಆದರೆ(*(FLOAT *)(ಗುಂಪು_beta) == 0.0){
	ಗುಂಪು_ಬಹುಲಕ=ಬಹುಲಕ | BLAS_SMALL_B0_OPT;
	ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((ಗುಂಪು_transb<<2)|ಗುಂಪು_transa));
      }ಇಲ್ಲದಿದ್ದರೆ{
	ಗುಂಪು_ಬಹುಲಕ=ಬಹುಲಕ | BLAS_SMALL_OPT;
	ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((ಗುಂಪು_transb<<2)|ಗುಂಪು_transa));
      }
#ಇಲ್ಲದಿದ್ದರೆ
      ಆದರೆ(((FLOAT *)(ಗುಂಪು_beta))[0] == 0.0 && ((FLOAT *)(ಗುಂಪು_beta))[1] == 0.0){
	ಗುಂಪು_ಬಹುಲಕ=ಬಹುಲಕ | BLAS_SMALL_B0_OPT;
	ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((ಗುಂಪು_transb<<2)|ಗುಂಪು_transa));
      }ಇಲ್ಲದಿದ್ದರೆ{
	ಗುಂಪು_ಬಹುಲಕ=ಬಹುಲಕ | BLAS_SMALL_OPT;
	ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((ಗುಂಪು_transb<<2)|ಗುಂಪು_transa));
      }

#ಆದರೆ_ಅಂತ್ಯ
      
    }ಇಲ್ಲದಿದ್ದರೆ{
#ಆದರೆ_ಅಂತ್ಯ
      ಗುಂಪು_routine=(ಶೂನ್ಯ*)(gemm[(ಗುಂಪು_transb<<2)|ಗುಂಪು_transa]);
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMALL_ಮಾತೃಕೆ_OPT
    }
#ಆದರೆ_ಅಂತ್ಯ    

    
    ಗಾಗಿ(j=0; j<ಗುಂಪು_ಗಾತ್ರ[i]; j++){
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].m=ಗುಂಪು_m;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].n=ಗುಂಪು_n;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].k=ಗುಂಪು_k;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].lda=ಗುಂಪು_lda;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].ldb=ಗುಂಪು_ldb;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].ldc=ಗುಂಪು_ldc;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].alpha=ಗುಂಪು_alpha;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].beta=ಗುಂಪು_beta;

#ಆದರೆ defined(CBLAS)
      ಆದರೆ (ಕ್ರಮ == Cblas_ಅಡ್ಡಸಾಲು_Major) {
	ವಾದಗಳು_array[ಎಣಿಕೆ_2].a=(b_array[ಮಾತೃಕೆ_idx+j]);
	ವಾದಗಳು_array[ಎಣಿಕೆ_2].b=(a_array[ಮಾತೃಕೆ_idx+j]);
      } ಇಲ್ಲದಿದ್ದರೆ
#ಆದರೆ_ಅಂತ್ಯ
      {
	ವಾದಗಳು_array[ಎಣಿಕೆ_2].a=(a_array[ಮಾತೃಕೆ_idx+j]);
	ವಾದಗಳು_array[ಎಣಿಕೆ_2].b=(b_array[ಮಾತೃಕೆ_idx+j]);
      }
      
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].c=(c_array[ಮಾತೃಕೆ_idx+j]);
      
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].routine_ಬಹುಲಕ=ಗುಂಪು_ಬಹುಲಕ;
      ವಾದಗಳು_array[ಎಣಿಕೆ_2].routine=ಗುಂಪು_routine;
#ವ್ಯಾಖ್ಯಾನಿತವಾದರೆ SMALL_ಮಾತೃಕೆ_OPT
      ಆದರೆ (!ಗುಂಪು_routine)
        ವಾದಗಳು_array[ಎಣಿಕೆ_2].routine=ಗುಂಪು_small_ಮಾತೃಕೆ_opt_routine;
#ಆದರೆ_ಅಂತ್ಯ      
      ಎಣಿಕೆ_2++;
    }
  }

  ಆದರೆ(ಎಣಿಕೆ_2>0){
    GEMM_BATCH_THREAD(ವಾದಗಳು_array,ಎಣಿಕೆ_2);
  }

  ಬಿಡುಗಡೆ(ವಾದಗಳು_array);
}
