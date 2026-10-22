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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव "common.h"

शून्य openblas_warning(पूर्णाङ्क verbose, ध्रुव वर्ण * सन्देशः);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#अञ्ञथा
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#अञ्ञथा
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#अञ्ञथा
#define ERROR_NAME "CGEMM_BATCH "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
स्थिर पूर्णाङ्क (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#सचे defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#अञ्ञथा
#define USE_SMALL_MATRIX_OPT 0
#endif

#सचे USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(table[idx]))
#अञ्ञथा
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(*(uintptr_t *)((वर्ण *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
स्थिर size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अञ्ञथा

स्थिर size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

स्थिर size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
शून्य NAME(वर्ण *transa_array, वर्ण *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_गणनम् = *gcount;
#अञ्ञथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_TRANSPOSE *  transa_array, गणना CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_गणनम्, blasint * group_size) {
#अञ्ञथा
	   शून्य * valpha_array,
	   शून्य ** va_array, blasint * lda_array,
	   शून्य ** vb_array, blasint * ldb_array,
	   शून्य * vbeta_array,
	   शून्य ** vc_array, blasint * ldc_array, blasint group_गणनम्, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=रिक्त;

  पूर्णाङ्क mode=0, group_mode=0;
  blasint total_अङ्कः=0;

  blasint i=0, j=0, matrix_idx=0, गणनम्=0;

  पूर्णाङ्क group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  शून्य * group_alpha, * group_beta;
  शून्य * group_routine=रिक्त;
#ifdef SMALL_MATRIX_OPT
  शून्य * group_small_matrix_opt_routine=रिक्त;
#endif
  
#सचे defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुण MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  अत्थाय(i=0; i<group_गणनम्; i++){
    total_अङ्कः+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_अङ्कः * परिमाण(blas_arg_t));
  
  सचे(args_array == रिक्त){
    openblas_warning(0, "memory alloc failed!\n");
    पटिनिवत्त;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अञ्ञथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अञ्ञथा
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अञ्ञथा
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  अत्थाय(i=0; i<group_गणनम्; matrix_idx+=group_size[i], i++){
    group_alpha = (शून्य *)&alpha_array[i * COMPSIZE];
    group_beta = (शून्य *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#सचे defined(CBLAS) 
    सचे (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      सचे (transa_array[i] == CblasNoTrans)     group_transa = 0;
      सचे (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      सचे (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      सचे (transa_array[i] == CblasConjTrans)   group_transa = 1;
#अञ्ञथा
      सचे (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      सचे (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      सचे (transb_array[i] == CblasNoTrans)     group_transb = 0;
      सचे (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      सचे (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      सचे (transb_array[i] == CblasConjTrans)   group_transb = 1;
#अञ्ञथा
      सचे (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      सचे (transb_array[i] == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      सचे (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      सचे (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      सचे (group_ldc < group_m) info = 13;
      सचे (group_ldb < group_nrowb)  info = 10;
      सचे (group_lda < group_nrowa)  info =  8;
      सचे (group_k < 0)        info =  5;
      सचे (group_n < 0)        info =  4;
      सचे (group_m < 0)        info =  3;
      सचे (group_transb < 0)        info =  2;
      सचे (group_transa < 0)        info =  1;

    }अञ्ञथा सचे (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      सचे (transb_array[i] == CblasNoTrans)     group_transa = 0;
      सचे (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      सचे (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      सचे (transb_array[i] == CblasConjTrans)   group_transa = 1;
#अञ्ञथा
      सचे (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      सचे (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      सचे (transa_array[i] == CblasNoTrans)     group_transb = 0;
      सचे (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      सचे (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      सचे (transa_array[i] == CblasConjTrans)   group_transb = 1;
#अञ्ञथा
      सचे (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      सचे (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#अञ्ञथा
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      सचे (transa_array[i] == 'N')   group_transa = 0;
      सचे (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      सचे (transa_array[i] == 'C')   group_transa = 1;
#अञ्ञथा
      सचे (transa_array[i] == 'C')   group_transa = 3;
#endif
      सचे (transb_array[i] == 'N')   group_transb = 0;
      सचे (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      सचे (transb_array[i] == 'C')   group_transb = 1;
#अञ्ञथा
      सचे (transb_array[i] == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      सचे (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      सचे (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      सचे (group_ldc < group_m) info = 13;
      सचे (group_ldb < group_nrowb)  info = 10;
      सचे (group_lda < group_nrowa)  info =  8;
      सचे (group_k < 0)        info =  5;
      सचे (group_n < 0)        info =  4;
      सचे (group_m < 0)        info =  3;
      सचे (group_transb < 0)        info =  2;
      सचे (group_transa < 0)        info =  1;      
#सचे defined(CBLAS)
    }
#endif

    सचे (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
      मुक्त(args_array);
      पटिनिवत्त;
    }

    सचे (group_m == 0 || group_n == 0) अनुवर्त;

    group_mode=mode;

#सचे defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुण) group_m * (द्विगुण) group_n * (द्विगुण) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#सचे !defined(COMPLEX)
    सचे (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#अञ्ञथा
    सचे (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=रिक्त;
#सचे !defined(COMPLEX)
      सचे(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अञ्ञथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#अञ्ञथा
      सचे(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अञ्ञथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }अञ्ञथा{
#endif
      group_routine=(शून्य*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
    अत्थाय(j=0; j<group_size[i]; j++){
      args_array[गणनम्].m=group_m;
      args_array[गणनम्].n=group_n;
      args_array[गणनम्].k=group_k;
      args_array[गणनम्].lda=group_lda;
      args_array[गणनम्].ldb=group_ldb;
      args_array[गणनम्].ldc=group_ldc;
      args_array[गणनम्].alpha=group_alpha;
      args_array[गणनम्].beta=group_beta;

#सचे defined(CBLAS)
      सचे (order == CblasRowMajor) {
	args_array[गणनम्].a=(b_array[matrix_idx+j]);
	args_array[गणनम्].b=(a_array[matrix_idx+j]);
      } अञ्ञथा
#endif
      {
	args_array[गणनम्].a=(a_array[matrix_idx+j]);
	args_array[गणनम्].b=(b_array[matrix_idx+j]);
      }
      
      args_array[गणनम्].c=(c_array[matrix_idx+j]);
      
      args_array[गणनम्].routine_mode=group_mode;
      args_array[गणनम्].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      सचे (!group_routine)
        args_array[गणनम्].routine=group_small_matrix_opt_routine;
#endif      
      गणनम्++;
    }
  }

  सचे(गणनम्>0){
    GEMM_BATCH_THREAD(args_array,गणनम्);
  }

  मुक्त(args_array);
}
