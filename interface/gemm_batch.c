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

#অন্তর্ভাৱ <stdio.h>
#অন্তর্ভাৱ <stdlib.h>
#অন্তর্ভাৱ "common.h"

শূন্য openblas_warning(পূর্ণাঙ্ক verbose, ধ্রুৱ ৱর্ণ * বাৰ্তা);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#অন্যথা
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#অন্যথা
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#নোচেদ্যদি defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#অন্যথা
#define ERROR_NAME "CGEMM_BATCH "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
স্থির পূর্ণাঙ্ক (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#যদি defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#অন্যথা
#define USE_SMALL_MATRIX_OPT 0
#endif

#যদি USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((শূন্য *)(table[idx]))
#অন্যথা
#define SMALL_KERNEL_ADDR(table, idx) ((শূন্য *)(*(uintptr_t *)((ৱর্ণ *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
স্থির size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


স্থির size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#অন্যথা

স্থির size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

স্থির size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (পূর্ণাঙ্ক (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
শূন্য NAME(ৱর্ণ *transa_array, ৱর্ণ *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_গুনতি = *gcount;
#অন্যথা

শূন্য CNAME(গণনা CBLAS_ORDER order, গণনা CBLAS_TRANSPOSE *  transa_array, গণনা CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_গুনতি, blasint * group_size) {
#অন্যথা
	   শূন্য * valpha_array,
	   শূন্য ** va_array, blasint * lda_array,
	   শূন্য ** vb_array, blasint * ldb_array,
	   শূন্য * vbeta_array,
	   শূন্য ** vc_array, blasint * ldc_array, blasint group_গুনতি, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=রিক্ত;

  পূর্ণাঙ্ক mode=0, group_mode=0;
  blasint total_অঙ্ক=0;

  blasint i=0, j=0, matrix_idx=0, গুনতি=0;

  পূর্ণাঙ্ক group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  শূন্য * group_alpha, * group_beta;
  শূন্য * group_routine=রিক্ত;
#ifdef SMALL_MATRIX_OPT
  শূন্য * group_small_matrix_opt_routine=রিক্ত;
#endif
  
#যদি defined (SMP) || defined(SMALL_MATRIX_OPT)
  দ্ৱিগুণ MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  বাবে(i=0; i<group_গুনতি; i++){
    total_অঙ্ক+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_অঙ্ক * পরিমাণ(blas_arg_t));
  
  যদি(args_array == রিক্ত){
    openblas_warning(0, "memory alloc failed!\n");
    ঘূৰাই;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#নোচেদ্যদি defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#অন্যথা
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#অন্যথা
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#নোচেদ্যদি defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#অন্যথা
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  বাবে(i=0; i<group_গুনতি; matrix_idx+=group_size[i], i++){
    group_alpha = (শূন্য *)&alpha_array[i * COMPSIZE];
    group_beta = (শূন্য *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#যদি defined(CBLAS) 
    যদি (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      যদি (transa_array[i] == CblasNoTrans)     group_transa = 0;
      যদি (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      যদি (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      যদি (transa_array[i] == CblasConjTrans)   group_transa = 1;
#অন্যথা
      যদি (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      যদি (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      যদি (transb_array[i] == CblasNoTrans)     group_transb = 0;
      যদি (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      যদি (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      যদি (transb_array[i] == CblasConjTrans)   group_transb = 1;
#অন্যথা
      যদি (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      যদি (transb_array[i] == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      যদি (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      যদি (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      যদি (group_ldc < group_m) info = 13;
      যদি (group_ldb < group_nrowb)  info = 10;
      যদি (group_lda < group_nrowa)  info =  8;
      যদি (group_k < 0)        info =  5;
      যদি (group_n < 0)        info =  4;
      যদি (group_m < 0)        info =  3;
      যদি (group_transb < 0)        info =  2;
      যদি (group_transa < 0)        info =  1;

    }অন্যথা যদি (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      যদি (transb_array[i] == CblasNoTrans)     group_transa = 0;
      যদি (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      যদি (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      যদি (transb_array[i] == CblasConjTrans)   group_transa = 1;
#অন্যথা
      যদি (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      যদি (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      যদি (transa_array[i] == CblasNoTrans)     group_transb = 0;
      যদি (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      যদি (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      যদি (transa_array[i] == CblasConjTrans)   group_transb = 1;
#অন্যথা
      যদি (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      যদি (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#অন্যথা
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      যদি (transa_array[i] == 'N')   group_transa = 0;
      যদি (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      যদি (transa_array[i] == 'C')   group_transa = 1;
#অন্যথা
      যদি (transa_array[i] == 'C')   group_transa = 3;
#endif
      যদি (transb_array[i] == 'N')   group_transb = 0;
      যদি (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      যদি (transb_array[i] == 'C')   group_transb = 1;
#অন্যথা
      যদি (transb_array[i] == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      যদি (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      যদি (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      যদি (group_ldc < group_m) info = 13;
      যদি (group_ldb < group_nrowb)  info = 10;
      যদি (group_lda < group_nrowa)  info =  8;
      যদি (group_k < 0)        info =  5;
      যদি (group_n < 0)        info =  4;
      যদি (group_m < 0)        info =  3;
      যদি (group_transb < 0)        info =  2;
      যদি (group_transa < 0)        info =  1;      
#যদি defined(CBLAS)
    }
#endif

    যদি (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, পরিমাণ(ERROR_NAME));
      মুক্ত(args_array);
      ঘূৰাই;
    }

    যদি (group_m == 0 || group_n == 0) অনুৱর্ত;

    group_mode=mode;

#যদি defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (দ্ৱিগুণ) group_m * (দ্ৱিগুণ) group_n * (দ্ৱিগুণ) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#যদি !defined(COMPLEX)
    যদি (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#অন্যথা
    যদি (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=রিক্ত;
#যদি !defined(COMPLEX)
      যদি(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }অন্যথা{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#অন্যথা
      যদি(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }অন্যথা{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }অন্যথা{
#endif
      group_routine=(শূন্য*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
    বাবে(j=0; j<group_size[i]; j++){
      args_array[গুনতি].m=group_m;
      args_array[গুনতি].n=group_n;
      args_array[গুনতি].k=group_k;
      args_array[গুনতি].lda=group_lda;
      args_array[গুনতি].ldb=group_ldb;
      args_array[গুনতি].ldc=group_ldc;
      args_array[গুনতি].alpha=group_alpha;
      args_array[গুনতি].beta=group_beta;

#যদি defined(CBLAS)
      যদি (order == CblasRowMajor) {
	args_array[গুনতি].a=(b_array[matrix_idx+j]);
	args_array[গুনতি].b=(a_array[matrix_idx+j]);
      } অন্যথা
#endif
      {
	args_array[গুনতি].a=(a_array[matrix_idx+j]);
	args_array[গুনতি].b=(b_array[matrix_idx+j]);
      }
      
      args_array[গুনতি].c=(c_array[matrix_idx+j]);
      
      args_array[গুনতি].routine_mode=group_mode;
      args_array[গুনতি].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      যদি (!group_routine)
        args_array[গুনতি].routine=group_small_matrix_opt_routine;
#endif      
      গুনতি++;
    }
  }

  যদি(গুনতি>0){
    GEMM_BATCH_THREAD(args_array,গুনতি);
  }

  মুক্ত(args_array);
}
