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

#include <stdio.h>
#include <stdlib.h>
#include "common.h"

void openblas_warning(ကိန်းပြည့် verbose, const var_na * msg);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#elif defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#မဟုတ်ပါက
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#elif defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#မဟုတ်ပါက
#define ERROR_NAME "CGEMM_BATCH "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
static ကိန်းပြည့် (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#အကယ်၍ defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#မဟုတ်ပါက
#define USE_SMALL_MATRIX_OPT 0
#endif

#အကယ်၍ USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(table[idx]))
#မဟုတ်ပါက
#define SMALL_KERNEL_ADDR(table, idx) ((void *)(*(uintptr_t *)((var_na *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
static size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


static size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (ကိန်းပြည့် (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (ကိန်းပြည့် (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#မဟုတ်ပါက

static size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

static size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (ကိန်းပြည့် (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (ကိန်းပြည့် (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(var_na *transa_array, var_na *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_count = *gcount;
#မဟုတ်ပါက

void CNAME(ga_nanaa CBLAS_ORDER order, ga_nanaa CBLAS_TRANSPOSE *  transa_array, ga_nanaa CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_count, blasint * group_size) {
#မဟုတ်ပါက
	   void * valpha_array,
	   void ** va_array, blasint * lda_array,
	   void ** vb_array, blasint * ldb_array,
	   void * vbeta_array,
	   void ** vc_array, blasint * ldc_array, blasint group_count, blasint * group_size) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=NULL;

  ကိန်းပြည့် mode=0, group_mode=0;
  blasint total_num=0;

  blasint i=0, j=0, matrix_idx=0, count=0;

  ကိန်းပြည့် group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  void * group_alpha, * group_beta;
  void * group_routine=NULL;
#ifdef SMALL_MATRIX_OPT
  void * group_small_matrix_opt_routine=NULL;
#endif
  
#အကယ်၍ defined (SMP) || defined(SMALL_MATRIX_OPT)
  dviguu_na MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  အတွက်(i=0; i<group_count; i++){
    total_num+=group_size[i];
  }

  args_array=(blas_arg_t *)malloc(total_num * parimaa_na(blas_arg_t));
  
  အကယ်၍(args_array == NULL){
    openblas_warning(0, "memory alloc failed!\n");
    ပြန်ပေး;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#မဟုတ်ပါက
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#မဟုတ်ပါက
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#elif defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#မဟုတ်ပါက
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  အတွက်(i=0; i<group_count; matrix_idx+=group_size[i], i++){
    group_alpha = (void *)&alpha_array[i * COMPSIZE];
    group_beta = (void *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#အကယ်၍ defined(CBLAS) 
    အကယ်၍ (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      အကယ်၍ (transa_array[i] == CblasNoTrans)     group_transa = 0;
      အကယ်၍ (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      အကယ်၍ (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      အကယ်၍ (transa_array[i] == CblasConjTrans)   group_transa = 1;
#မဟုတ်ပါက
      အကယ်၍ (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      အကယ်၍ (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      အကယ်၍ (transb_array[i] == CblasNoTrans)     group_transb = 0;
      အကယ်၍ (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      အကယ်၍ (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      အကယ်၍ (transb_array[i] == CblasConjTrans)   group_transb = 1;
#မဟုတ်ပါက
      အကယ်၍ (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      အကယ်၍ (transb_array[i] == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      အကယ်၍ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      အကယ်၍ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      အကယ်၍ (group_ldc < group_m) info = 13;
      အကယ်၍ (group_ldb < group_nrowb)  info = 10;
      အကယ်၍ (group_lda < group_nrowa)  info =  8;
      အကယ်၍ (group_k < 0)        info =  5;
      အကယ်၍ (group_n < 0)        info =  4;
      အကယ်၍ (group_m < 0)        info =  3;
      အကယ်၍ (group_transb < 0)        info =  2;
      အကယ်၍ (group_transa < 0)        info =  1;

    }မဟုတ်ပါက အကယ်၍ (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      အကယ်၍ (transb_array[i] == CblasNoTrans)     group_transa = 0;
      အကယ်၍ (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      အကယ်၍ (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      အကယ်၍ (transb_array[i] == CblasConjTrans)   group_transa = 1;
#မဟုတ်ပါက
      အကယ်၍ (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      အကယ်၍ (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      အကယ်၍ (transa_array[i] == CblasNoTrans)     group_transb = 0;
      အကယ်၍ (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      အကယ်၍ (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      အကယ်၍ (transa_array[i] == CblasConjTrans)   group_transb = 1;
#မဟုတ်ပါက
      အကယ်၍ (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      အကယ်၍ (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#မဟုတ်ပါက
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      အကယ်၍ (transa_array[i] == 'N')   group_transa = 0;
      အကယ်၍ (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      အကယ်၍ (transa_array[i] == 'C')   group_transa = 1;
#မဟုတ်ပါက
      အကယ်၍ (transa_array[i] == 'C')   group_transa = 3;
#endif
      အကယ်၍ (transb_array[i] == 'N')   group_transb = 0;
      အကယ်၍ (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      အကယ်၍ (transb_array[i] == 'C')   group_transb = 1;
#မဟုတ်ပါက
      အကယ်၍ (transb_array[i] == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      အကယ်၍ (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      အကယ်၍ (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      အကယ်၍ (group_ldc < group_m) info = 13;
      အကယ်၍ (group_ldb < group_nrowb)  info = 10;
      အကယ်၍ (group_lda < group_nrowa)  info =  8;
      အကယ်၍ (group_k < 0)        info =  5;
      အကယ်၍ (group_n < 0)        info =  4;
      အကယ်၍ (group_m < 0)        info =  3;
      အကယ်၍ (group_transb < 0)        info =  2;
      အကယ်၍ (group_transa < 0)        info =  1;      
#အကယ်၍ defined(CBLAS)
    }
#endif

    အကယ်၍ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, parimaa_na(ERROR_NAME));
      free(args_array);
      ပြန်ပေး;
    }

    အကယ်၍ (group_m == 0 || group_n == 0) continue;

    group_mode=mode;

#အကယ်၍ defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (dviguu_na) group_m * (dviguu_na) group_n * (dviguu_na) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#အကယ်၍ !defined(COMPLEX)
    အကယ်၍ (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#မဟုတ်ပါက
    အကယ်၍ (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=NULL;
#အကယ်၍ !defined(COMPLEX)
      အကယ်၍(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }မဟုတ်ပါက{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#မဟုတ်ပါက
      အကယ်၍(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }မဟုတ်ပါက{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }မဟုတ်ပါက{
#endif
      group_routine=(void*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
    အတွက်(j=0; j<group_size[i]; j++){
      args_array[count].m=group_m;
      args_array[count].n=group_n;
      args_array[count].k=group_k;
      args_array[count].lda=group_lda;
      args_array[count].ldb=group_ldb;
      args_array[count].ldc=group_ldc;
      args_array[count].alpha=group_alpha;
      args_array[count].beta=group_beta;

#အကယ်၍ defined(CBLAS)
      အကယ်၍ (order == CblasRowMajor) {
	args_array[count].a=(b_array[matrix_idx+j]);
	args_array[count].b=(a_array[matrix_idx+j]);
      } မဟုတ်ပါက
#endif
      {
	args_array[count].a=(a_array[matrix_idx+j]);
	args_array[count].b=(b_array[matrix_idx+j]);
      }
      
      args_array[count].c=(c_array[matrix_idx+j]);
      
      args_array[count].routine_mode=group_mode;
      args_array[count].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      အကယ်၍ (!group_routine)
        args_array[count].routine=group_small_matrix_opt_routine;
#endif      
      count++;
    }
  }

  အကယ်၍(count>0){
    GEMM_BATCH_THREAD(args_array,count);
  }

  free(args_array);
}
