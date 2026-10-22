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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश "common.h"

व्योम openblas_warning(पूर्णांक verbose, स्थिर अक्षर * संदेश);

#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMM_BATCH "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMM_BATCH "
#घोषणा GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#घोषणा ERROR_NAME "SGEMM_BATCH "
#घोषणा GEMM_BATCH_THREAD sgemm_batch_thread
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMM_BATCH "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMM_BATCH "
#घोषणा GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#घोषणा ERROR_NAME "CGEMM_BATCH "
#घोषणा GEMM_BATCH_THREAD cgemm_batch_thread
#पूर्ण_अगर
#पूर्ण_अगर
अटल पूर्णांक (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#अगर defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#घोषणा USE_SMALL_MATRIX_OPT 1
#अन्यथा
#घोषणा USE_SMALL_MATRIX_OPT 0
#पूर्ण_अगर

#अगर USE_SMALL_MATRIX_OPT
#अगर_अघोषित DYNAMIC_ARCH
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(table[idx]))
#अन्यथा
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(*(uintptr_t *)((अक्षर *)gotoblas + (माप_प्रकार)(table[idx]))))
#पूर्ण_अगर


#अगर_अघोषित COMPLEX
अटल माप_प्रकार gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


अटल माप_प्रकार gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#घोषणा GEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#घोषणा GEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

अटल माप_प्रकार zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

अटल माप_प्रकार zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#घोषणा ZGEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#घोषणा ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#पूर्ण_अगर
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS
व्योम NAME(अक्षर *transa_array, अक्षर *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * group_आकार) {
       blasint group_गिनती = *gcount;
#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_TRANSPOSE *  transa_array, क्रमागत CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#अगर_अघोषित COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint group_गिनती, blasint * group_आकार) {
#अन्यथा
	   व्योम * valpha_array,
	   व्योम ** va_array, blasint * lda_array,
	   व्योम ** vb_array, blasint * ldb_array,
	   व्योम * vbeta_array,
	   व्योम ** vc_array, blasint * ldc_array, blasint group_गिनती, blasint * group_आकार) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#पूर्ण_अगर
#पूर्ण_अगर
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=शून्य;

  पूर्णांक mode=0, group_mode=0;
  blasint योग_अंक=0;

  blasint i=0, j=0, matrix_idx=0, गिनती=0;

  पूर्णांक group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  व्योम * group_alpha, * group_beta;
  व्योम * group_routine=शून्य;
#अगर_घोषित SMALL_MATRIX_OPT
  व्योम * group_small_matrix_opt_routine=शून्य;
#पूर्ण_अगर
  
#अगर defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुन MNK;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;
    
  क्रम(i=0; i<group_गिनती; i++){
    योग_अंक+=group_आकार[i];
  }

  args_array=(blas_arg_t *)दो_स्मृति(योग_अंक * माप(blas_arg_t));
  
  अगर(args_array == शून्य){
    openblas_warning(0, "memory alloc failed!\n");
    वापस;
  }

#अगर_घोषित SMP
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
  
  क्रम(i=0; i<group_गिनती; matrix_idx+=group_आकार[i], i++){
    group_alpha = (व्योम *)&alpha_array[i * COMPSIZE];
    group_beta = (व्योम *)&beta_array[i * COMPSIZE];

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#अगर defined(CBLAS) 
    अगर (order == CblasColMajor) {
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      अगर (transa_array[i] == CblasNoTrans)     group_transa = 0;
      अगर (transa_array[i] == CblasTrans)       group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (transa_array[i] == CblasConjNoTrans) group_transa = 0;
      अगर (transa_array[i] == CblasConjTrans)   group_transa = 1;
#अन्यथा
      अगर (transa_array[i] == CblasConjNoTrans) group_transa = 2;
      अगर (transa_array[i] == CblasConjTrans)   group_transa = 3;
#पूर्ण_अगर
      अगर (transb_array[i] == CblasNoTrans)     group_transb = 0;
      अगर (transb_array[i] == CblasTrans)       group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (transb_array[i] == CblasConjNoTrans) group_transb = 0;
      अगर (transb_array[i] == CblasConjTrans)   group_transb = 1;
#अन्यथा
      अगर (transb_array[i] == CblasConjNoTrans) group_transb = 2;
      अगर (transb_array[i] == CblasConjTrans)   group_transb = 3;
#पूर्ण_अगर
      group_nrowa = group_m;
      अगर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      अगर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      अगर (group_ldc < group_m) info = 13;
      अगर (group_ldb < group_nrowb)  info = 10;
      अगर (group_lda < group_nrowa)  info =  8;
      अगर (group_k < 0)        info =  5;
      अगर (group_n < 0)        info =  4;
      अगर (group_m < 0)        info =  3;
      अगर (group_transb < 0)        info =  2;
      अगर (group_transa < 0)        info =  1;

    }अन्यथा अगर (order == CblasRowMajor) {
      
      group_m = n_array[i];
      group_n = m_array[i];
      group_k = k_array[i];

      group_lda = ldb_array[i];
      group_ldb = lda_array[i];
      group_ldc = ldc_array[i];

      अगर (transb_array[i] == CblasNoTrans)     group_transa = 0;
      अगर (transb_array[i] == CblasTrans)       group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (transb_array[i] == CblasConjNoTrans) group_transa = 0;
      अगर (transb_array[i] == CblasConjTrans)   group_transa = 1;
#अन्यथा
      अगर (transb_array[i] == CblasConjNoTrans) group_transa = 2;
      अगर (transb_array[i] == CblasConjTrans)   group_transa = 3;
#पूर्ण_अगर
      अगर (transa_array[i] == CblasNoTrans)     group_transb = 0;
      अगर (transa_array[i] == CblasTrans)       group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (transa_array[i] == CblasConjNoTrans) group_transb = 0;
      अगर (transa_array[i] == CblasConjTrans)   group_transb = 1;
#अन्यथा
      अगर (transa_array[i] == CblasConjNoTrans) group_transb = 2;
      अगर (transa_array[i] == CblasConjTrans)   group_transb = 3;
#पूर्ण_अगर

#अन्यथा
      group_m = m_array[i];
      group_n = n_array[i];
      group_k = k_array[i];

      group_lda = lda_array[i];
      group_ldb = ldb_array[i];
      group_ldc = ldc_array[i];
      
      अगर (transa_array[i] == 'N')   group_transa = 0;
      अगर (transa_array[i] == 'T')   group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (transa_array[i] == 'C')   group_transa = 1;
#अन्यथा
      अगर (transa_array[i] == 'C')   group_transa = 3;
#पूर्ण_अगर
      अगर (transb_array[i] == 'N')   group_transb = 0;
      अगर (transb_array[i] == 'T')   group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (transb_array[i] == 'C')   group_transb = 1;
#अन्यथा
      अगर (transb_array[i] == 'C')   group_transb = 3;
#पूर्ण_अगर
#पूर्ण_अगर

      group_nrowa = group_m;
      अगर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      अगर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      अगर (group_ldc < group_m) info = 13;
      अगर (group_ldb < group_nrowb)  info = 10;
      अगर (group_lda < group_nrowa)  info =  8;
      अगर (group_k < 0)        info =  5;
      अगर (group_n < 0)        info =  4;
      अगर (group_m < 0)        info =  3;
      अगर (group_transb < 0)        info =  2;
      अगर (group_transa < 0)        info =  1;      
#अगर defined(CBLAS)
    }
#पूर्ण_अगर

    अगर (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
      मुक्त(args_array);
      वापस;
    }

    अगर (group_m == 0 || group_n == 0) जारी;

    group_mode=mode;

#अगर defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुन) group_m * (द्विगुन) group_n * (द्विगुन) group_k;
#पूर्ण_अगर

#अगर_घोषित SMALL_MATRIX_OPT
#अगर !defined(COMPLEX)
    अगर (*(FLOAT *)(group_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(group_alpha), *(FLOAT *)(group_beta))){
#अन्यथा
    अगर (((FLOAT *)(group_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(group_alpha))[0], ((FLOAT *)(group_alpha))[1],
                                 ((FLOAT *)(group_beta))[0], ((FLOAT *)(group_beta))[1])){
#पूर्ण_अगर
      group_routine=शून्य;
#अगर !defined(COMPLEX)
      अगर(*(FLOAT *)(group_beta) == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#अन्यथा
      अगर(((FLOAT *)(group_beta))[0] == 0.0 && ((FLOAT *)(group_beta))[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#पूर्ण_अगर
      
    }अन्यथा{
#पूर्ण_अगर
      group_routine=(व्योम*)(gemm[(group_transb<<2)|group_transa]);
#अगर_घोषित SMALL_MATRIX_OPT
    }
#पूर्ण_अगर    

    
    क्रम(j=0; j<group_आकार[i]; j++){
      args_array[गिनती].m=group_m;
      args_array[गिनती].n=group_n;
      args_array[गिनती].k=group_k;
      args_array[गिनती].lda=group_lda;
      args_array[गिनती].ldb=group_ldb;
      args_array[गिनती].ldc=group_ldc;
      args_array[गिनती].alpha=group_alpha;
      args_array[गिनती].beta=group_beta;

#अगर defined(CBLAS)
      अगर (order == CblasRowMajor) {
	args_array[गिनती].a=(b_array[matrix_idx+j]);
	args_array[गिनती].b=(a_array[matrix_idx+j]);
      } अन्यथा
#पूर्ण_अगर
      {
	args_array[गिनती].a=(a_array[matrix_idx+j]);
	args_array[गिनती].b=(b_array[matrix_idx+j]);
      }
      
      args_array[गिनती].c=(c_array[matrix_idx+j]);
      
      args_array[गिनती].routine_mode=group_mode;
      args_array[गिनती].routine=group_routine;
#अगर_घोषित SMALL_MATRIX_OPT
      अगर (!group_routine)
        args_array[गिनती].routine=group_small_matrix_opt_routine;
#पूर्ण_अगर      
      गिनती++;
    }
  }

  अगर(गिनती>0){
    GEMM_BATCH_THREAD(args_array,गिनती);
  }

  मुक्त(args_array);
}
