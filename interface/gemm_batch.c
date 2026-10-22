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

शून्य openblas_warning(पूर्णांक verbose, स्थिरांक अक्षर * संदेश);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#define ERROR_NAME "SGEMM_BATCH "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH "
#किंवाजर defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#define ERROR_NAME "CGEMM_BATCH "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
स्थिर पूर्णांक (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#जर defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#अन्यथा
#define USE_SMALL_MATRIX_OPT 0
#endif

#जर USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(table[idx]))
#अन्यथा
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(*(uintptr_t *)((अक्षर *)gotoblas + (माप_प्रकार)(table[idx]))))
#endif


#ifndef COMPLEX
स्थिर माप_प्रकार gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर माप_प्रकार gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

स्थिर माप_प्रकार zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

स्थिर माप_प्रकार zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
शून्य NAME(अक्षर *transa_array, अक्षर *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_size) {
       blasint group_मोजणी = *gcount;
#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_TRANSPOSE *  transa_array, गणना CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#ifndef COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_मोजणी, blasint * group_size) {
#अन्यथा
	   शून्य * valpha_array,
	   शून्य ** va_array, blasint * lda_array,
	   शून्य ** vb_array, blasint * ldb_array,
	   शून्य * vbeta_array,
	   शून्य ** vc_array, blasint * ldc_array, blasint group_मोजणी, blasint * group_size) {

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

  पूर्णांक mode=0, group_mode=0;
  blasint total_अंक=0;

  blasint i=0, j=0, matrix_idx=0, मोजणी=0;

  पूर्णांक group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  शून्य * group_alpha, * group_beta;
  शून्य * group_routine=रिक्त;
#ifdef SMALL_MATRIX_OPT
  शून्य * group_small_matrix_opt_routine=रिक्त;
#endif
  
#जर defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुण MNK;
#endif

  PRINT_DEBUG_CNAME;
    
  साठी(i=0; i<group_मोजणी; i++){
    total_अंक+=group_size[i];
  }

  args_array=(blas_arg_t *)दो_स्मृति(total_अंक * आकार(blas_arg_t));
  
  जर(args_array == रिक्त){
    openblas_warning(0, "memory alloc failed!\n");
    परत;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#किंवाजर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#किंवाजर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  साठी(i=0; i<group_मोजणी; matrix_idx+=group_size[i], i++){
    group_alpha = (शून्य *)&alpha_array[i * COMPSIZE];
    group_beta = (शून्य *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#जर defined(CBLAS) 
    जर (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      जर (transa_array[i] == CblasNoTrans)     group_transa = 0;
      जर (transa_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      जर (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      जर (transa_array[i] == CblasConjTrans)   group_transa = 1;
#अन्यथा
      जर (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      जर (transa_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      जर (transb_array[i] == CblasNoTrans)     group_transb = 0;
      जर (transb_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      जर (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      जर (transb_array[i] == CblasConjTrans)   group_transb = 1;
#अन्यथा
      जर (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      जर (transb_array[i] == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      जर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      जर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      जर (group_ldc < group_m) info = 13;
      जर (group_ldb < group_nrowb)  info = 10;
      जर (group_lda < group_nrowa)  info =  8;
      जर (group_k < 0)        info =  5;
      जर (group_n < 0)        info =  4;
      जर (group_m < 0)        info =  3;
      जर (group_transb < 0)        info =  2;
      जर (group_transa < 0)        info =  1;

    }अन्यथा जर (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      जर (transb_array[i] == CblasNoTrans)     group_transa = 0;
      जर (transb_array[i] == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      जर (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      जर (transb_array[i] == CblasConjTrans)   group_transa = 1;
#अन्यथा
      जर (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      जर (transb_array[i] == CblasConjTrans)   group_transa = 3;
#endif
      जर (transa_array[i] == CblasNoTrans)     group_transb = 0;
      जर (transa_array[i] == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      जर (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      जर (transa_array[i] == CblasConjTrans)   group_transb = 1;
#अन्यथा
      जर (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      जर (transa_array[i] == CblasConjTrans)   group_transb = 3;
#endif

#अन्यथा
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      जर (transa_array[i] == 'N')   group_transa = 0;
      जर (transa_array[i] == 'T')   group_transa = 1;
#ifndef COMPLEX
      जर (transa_array[i] == 'C')   group_transa = 1;
#अन्यथा
      जर (transa_array[i] == 'C')   group_transa = 3;
#endif
      जर (transb_array[i] == 'N')   group_transb = 0;
      जर (transb_array[i] == 'T')   group_transb = 1;
#ifndef COMPLEX
      जर (transb_array[i] == 'C')   group_transb = 1;
#अन्यथा
      जर (transb_array[i] == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      जर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      जर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      जर (group_ldc < group_m) info = 13;
      जर (group_ldb < group_nrowb)  info = 10;
      जर (group_lda < group_nrowa)  info =  8;
      जर (group_k < 0)        info =  5;
      जर (group_n < 0)        info =  4;
      जर (group_m < 0)        info =  3;
      जर (group_transb < 0)        info =  2;
      जर (group_transa < 0)        info =  1;      
#जर defined(CBLAS)
    }
#endif

    जर (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, आकार(ERROR_NAME));
      मुक्त(args_array);
      परत;
    }

    जर (group_m == 0 || group_n == 0) सुरूठेवा;

    group_mode=mode;

#जर defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुण) group_m * (द्विगुण) group_n * (द्विगुण) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#जर !defined(COMPLEX)
    जर (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#अन्यथा
    जर (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#endif
      group_routine=रिक्त;
#जर !defined(COMPLEX)
      जर(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#अन्यथा
      जर(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }अन्यथा{
#endif
      group_routine=(शून्य*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
    साठी(j=0; j<group_size[i]; j++){
      args_array[मोजणी].m=group_m;
      args_array[मोजणी].n=group_n;
      args_array[मोजणी].k=group_k;
      args_array[मोजणी].lda=group_lda;
      args_array[मोजणी].ldb=group_ldb;
      args_array[मोजणी].ldc=group_ldc;
      args_array[मोजणी].alpha=group_alpha;
      args_array[मोजणी].beta=group_beta;

#जर defined(CBLAS)
      जर (order == CblasRowMajor) {
	args_array[मोजणी].a=(b_array[matrix_idx+j]);
	args_array[मोजणी].b=(a_array[matrix_idx+j]);
      } अन्यथा
#endif
      {
	args_array[मोजणी].a=(a_array[matrix_idx+j]);
	args_array[मोजणी].b=(b_array[matrix_idx+j]);
      }
      
      args_array[मोजणी].c=(c_array[matrix_idx+j]);
      
      args_array[मोजणी].routine_mode=group_mode;
      args_array[मोजणी].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      जर (!group_routine)
        args_array[मोजणी].routine=group_small_matrix_opt_routine;
#endif      
      मोजणी++;
    }
  }

  जर(मोजणी>0){
    GEMM_BATCH_THREAD(args_array,मोजणी);
  }

  मुक्त(args_array);
}
