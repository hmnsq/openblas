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

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن "common.h"

فراغ openblas_warning(صحيح verbose, ثابت حرف * msg);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#אחרת
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#אחרת
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#אחרת
#define ERROR_NAME "CGEMM_BATCH "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#אם defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#אחרת
#define USE_SMALL_MATRIX_OPT 0
#endif

#אם USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((فراغ *)(table[idx]))
#אחרת
#define SMALL_KERNEL_ADDR(table, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#אחרת

ساكن size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ساكن size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
فراغ NAME(حرف *transa_array, حرف *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_count = *gcount;
#אחרת

فراغ CNAME(تعداد CBLAS_ORDER order, تعداد CBLAS_TRANSPOSE *  transa_array, تعداد CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_count, blasint * group_size) {
#אחרת
	   فراغ * valpha_array,
	   فراغ ** va_array, blasint * lda_array,
	   فراغ ** vb_array, blasint * ldb_array,
	   فراغ * vbeta_array,
	   فراغ ** vc_array, blasint * ldc_array, blasint group_count, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=عدم;

  صحيح mode=0, group_mode=0;
  blasint total_num=0;

  blasint i=0, j=0, matrix_idx=0, count=0;

  صحيح group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  فراغ * group_alpha, * group_beta;
  فراغ * group_routine=عدم;
#ifdef SMALL_MATRIX_OPT
  فراغ * group_small_matrix_opt_routine=عدم;
#endif
  
#אם defined (SMP) || defined(SMALL_MATRIX_OPT)
  مزدوج MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  עבור(i=0; i<group_count; i++){
    total_num+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_num * حجم(blas_arg_t));
  
  אם(args_array == عدم){
    openblas_warning(0, "memory alloc failed!\n");
    החזר;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#אחרת
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#אחרת
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#אחרת
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  עבור(i=0; i<group_count; matrix_idx+=group_size[i], i++){
    group_alpha = (فراغ *)&alpha_array[i * COMPSIZE];
    group_beta = (فراغ *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#אם defined(CBLAS) 
    אם (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      אם (transa_array[i] == CblasNoTrans)     group_transa = 0;
      אם (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      אם (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      אם (transa_array[i] == CblasConjTrans)   group_transa = 1;
#אחרת
      אם (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      אם (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      אם (transb_array[i] == CblasNoTrans)     group_transb = 0;
      אם (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      אם (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      אם (transb_array[i] == CblasConjTrans)   group_transb = 1;
#אחרת
      אם (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      אם (transb_array[i] == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      אם (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      אם (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      אם (group_ldc < group_m) info = 13;
      אם (group_ldb < group_nrowb)  info = 10;
      אם (group_lda < group_nrowa)  info =  8;
      אם (group_k < 0)        info =  5;
      אם (group_n < 0)        info =  4;
      אם (group_m < 0)        info =  3;
      אם (group_transb < 0)        info =  2;
      אם (group_transa < 0)        info =  1;

    }אחרת אם (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      אם (transb_array[i] == CblasNoTrans)     group_transa = 0;
      אם (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      אם (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      אם (transb_array[i] == CblasConjTrans)   group_transa = 1;
#אחרת
      אם (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      אם (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      אם (transa_array[i] == CblasNoTrans)     group_transb = 0;
      אם (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      אם (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      אם (transa_array[i] == CblasConjTrans)   group_transb = 1;
#אחרת
      אם (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      אם (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#אחרת
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      אם (transa_array[i] == 'N')   group_transa = 0;
      אם (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      אם (transa_array[i] == 'C')   group_transa = 1;
#אחרת
      אם (transa_array[i] == 'C')   group_transa = 3;
#endif
      אם (transb_array[i] == 'N')   group_transb = 0;
      אם (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      אם (transb_array[i] == 'C')   group_transb = 1;
#אחרת
      אם (transb_array[i] == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      אם (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      אם (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      אם (group_ldc < group_m) info = 13;
      אם (group_ldb < group_nrowb)  info = 10;
      אם (group_lda < group_nrowa)  info =  8;
      אם (group_k < 0)        info =  5;
      אם (group_n < 0)        info =  4;
      אם (group_m < 0)        info =  3;
      אם (group_transb < 0)        info =  2;
      אם (group_transa < 0)        info =  1;      
#אם defined(CBLAS)
    }
#endif

    אם (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      حرر(args_array);
      החזר;
    }

    אם (group_m == 0 || group_n == 0) استمر;

    group_mode=mode;

#אם defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (مزدوج) group_m * (مزدوج) group_n * (مزدوج) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#אם !defined(COMPLEX)
    אם (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#אחרת
    אם (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=عدم;
#אם !defined(COMPLEX)
      אם(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }אחרת{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#אחרת
      אם(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }אחרת{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }אחרת{
#endif
      group_routine=(فراغ*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
    עבור(j=0; j<group_size[i]; j++){
      args_array[count].m=group_m;
      args_array[count].n=group_n;
      args_array[count].k=group_k;
      args_array[count].lda=group_lda;
      args_array[count].ldb=group_ldb;
      args_array[count].ldc=group_ldc;
      args_array[count].alpha=group_alpha;
      args_array[count].beta=group_beta;

#אם defined(CBLAS)
      אם (order == CblasRowMajor) {
	args_array[count].a=(b_array[matrix_idx+j]);
	args_array[count].b=(a_array[matrix_idx+j]);
      } אחרת
#endif
      {
	args_array[count].a=(a_array[matrix_idx+j]);
	args_array[count].b=(b_array[matrix_idx+j]);
      }
      
      args_array[count].c=(c_array[matrix_idx+j]);
      
      args_array[count].routine_mode=group_mode;
      args_array[count].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      אם (!group_routine)
        args_array[count].routine=group_small_matrix_opt_routine;
#endif      
      count++;
    }
  }

  אם(count>0){
    GEMM_BATCH_THREAD(args_array,count);
  }

  حرر(args_array);
}
