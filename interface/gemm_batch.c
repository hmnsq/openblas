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

#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGEMM_BATCH "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEMM_BATCH "
#परिभाषय GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#परिभाषय ERROR_NAME "SGEMM_BATCH "
#परिभाषय GEMM_BATCH_THREAD sgemm_batch_thread
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMM_BATCH "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMM_BATCH "
#परिभाषय GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#परिभाषय ERROR_NAME "CGEMM_BATCH "
#परिभाषय GEMM_BATCH_THREAD cgemm_batch_thread
#यदि_अन्तः
#यदि_अन्तः
स्थिर पूर्णाङ्क (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#यदि defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#परिभाषय USE_SMALL_MATRIX_OPT 1
#अन्यथा
#परिभाषय USE_SMALL_MATRIX_OPT 0
#यदि_अन्तः

#यदि USE_SMALL_MATRIX_OPT
#यदि_अपरिभाषितम् DYNAMIC_ARCH
#परिभाषय SMALL_KERNEL_ADDR(सारणी, idx) ((शून्य *)(सारणी[idx]))
#अन्यथा
#परिभाषय SMALL_KERNEL_ADDR(सारणी, idx) ((शून्य *)(*(uintptr_t *)((वर्ण *)gotoblas + (size_t)(सारणी[idx]))))
#यदि_अन्तः


#यदि_अपरिभाषितम् COMPLEX
स्थिर size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#परिभाषय GEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#परिभाषय GEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

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

#परिभाषय ZGEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#परिभाषय ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#यदि_अन्तः
#यदि_अन्तः

OPENBLAS_EXPORT

#यदि_अपरिभाषितम् CBLAS
शून्य NAME(वर्ण *transa_array, वर्ण *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * गणः_आकार) {
       blasint गणः_गणनम् = *gcount;
#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_TRANSPOSE *  transa_array, गणना CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint गणः_गणनम्, blasint * गणः_आकार) {
#अन्यथा
	   शून्य * valpha_array,
	   शून्य ** va_array, blasint * lda_array,
	   शून्य ** vb_array, blasint * ldb_array,
	   शून्य * vbeta_array,
	   शून्य ** vc_array, blasint * ldc_array, blasint गणः_गणनम्, blasint * गणः_आकार) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#यदि_अन्तः
#यदि_अन्तः
  BLASLONG गणः_m, गणः_n, गणः_k;
  BLASLONG गणः_lda, गणः_ldb, गणः_ldc;

  blas_arg_t * तर्काः_array=रिक्त;

  पूर्णाङ्क बहुलकम्=0, गणः_बहुलकम्=0;
  blasint योग_अङ्कः=0;

  blasint i=0, j=0, matrix_idx=0, गणनम्=0;

  पूर्णाङ्क गणः_transa, गणः_transb;
  BLASLONG गणः_nrowa, गणः_nrowb;
  blasint info;

  शून्य * गणः_alpha, * गणः_beta;
  शून्य * गणः_routine=रिक्त;
#यदि_परिभाषितम् SMALL_MATRIX_OPT
  शून्य * गणः_small_matrix_opt_routine=रिक्त;
#यदि_अन्तः
  
#यदि defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुण MNK;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;
    
  प्रति(i=0; i<गणः_गणनम्; i++){
    योग_अङ्कः+=गणः_आकार[i];
  }

  तर्काः_array=(blas_arg_t *)स्मृतिं_देहि(योग_अङ्कः * परिमाण(blas_arg_t));
  
  यदि(तर्काः_array == रिक्त){
    openblas_warning(0, "memory alloc failed!\n");
    प्रत्यागम;
  }

#यदि_परिभाषितम् SMP
#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
  बहुलकम्  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  बहुलकम्  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  बहुलकम्  =  BLAS_SINGLE  | BLAS_REAL;
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
  बहुलकम्  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  बहुलकम्  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  बहुलकम्  =  BLAS_SINGLE  | BLAS_COMPLEX;
#यदि_अन्तः
#यदि_अन्तः
#यदि_अन्तः
  
  प्रति(i=0; i<गणः_गणनम्; matrix_idx+=गणः_आकार[i], i++){
    गणः_alpha = (शून्य *)&alpha_array[i * COMPSIZE];
    गणः_beta = (शून्य *)&beta_array[i * COMPSIZE];

    गणः_m = गणः_n = गणः_k = 0;
    गणः_lda = गणः_ldb = गणः_ldc = 0;
    गणः_transa = -1;
    गणः_transb = -1;
    info   = 0;
   
#यदि defined(CBLAS) 
    यदि (क्रमः == CblasColMajor) {
      गणः_m = m_array[i];
      गणः_n = n_array[i];
      गणः_k = k_array[i];

      गणः_lda = lda_array[i];
      गणः_ldb = ldb_array[i];
      गणः_ldc = ldc_array[i];
      
      यदि (transa_array[i] == CblasNoTrans)     गणः_transa = 0;
      यदि (transa_array[i] == CblasTrans)       गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa_array[i] == CblasConjNoTrans) गणः_transa = 0;
      यदि (transa_array[i] == CblasConjTrans)   गणः_transa = 1;
#अन्यथा
      यदि (transa_array[i] == CblasConjNoTrans) गणः_transa = 2;
      यदि (transa_array[i] == CblasConjTrans)   गणः_transa = 3;
#यदि_अन्तः
      यदि (transb_array[i] == CblasNoTrans)     गणः_transb = 0;
      यदि (transb_array[i] == CblasTrans)       गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transb_array[i] == CblasConjNoTrans) गणः_transb = 0;
      यदि (transb_array[i] == CblasConjTrans)   गणः_transb = 1;
#अन्यथा
      यदि (transb_array[i] == CblasConjNoTrans) गणः_transb = 2;
      यदि (transb_array[i] == CblasConjTrans)   गणः_transb = 3;
#यदि_अन्तः
      गणः_nrowa = गणः_m;
      यदि (गणः_transa & 1) गणः_nrowa = गणः_k;
      गणः_nrowb = गणः_k;
      यदि (गणः_transb & 1) गणः_nrowb = गणः_n;

      info=-1;
	
      यदि (गणः_ldc < गणः_m) info = 13;
      यदि (गणः_ldb < गणः_nrowb)  info = 10;
      यदि (गणः_lda < गणः_nrowa)  info =  8;
      यदि (गणः_k < 0)        info =  5;
      यदि (गणः_n < 0)        info =  4;
      यदि (गणः_m < 0)        info =  3;
      यदि (गणः_transb < 0)        info =  2;
      यदि (गणः_transa < 0)        info =  1;

    }अन्यथा यदि (क्रमः == CblasRowMajor) {
      
      गणः_m = n_array[i];
      गणः_n = m_array[i];
      गणः_k = k_array[i];

      गणः_lda = ldb_array[i];
      गणः_ldb = lda_array[i];
      गणः_ldc = ldc_array[i];

      यदि (transb_array[i] == CblasNoTrans)     गणः_transa = 0;
      यदि (transb_array[i] == CblasTrans)       गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transb_array[i] == CblasConjNoTrans) गणः_transa = 0;
      यदि (transb_array[i] == CblasConjTrans)   गणः_transa = 1;
#अन्यथा
      यदि (transb_array[i] == CblasConjNoTrans) गणः_transa = 2;
      यदि (transb_array[i] == CblasConjTrans)   गणः_transa = 3;
#यदि_अन्तः
      यदि (transa_array[i] == CblasNoTrans)     गणः_transb = 0;
      यदि (transa_array[i] == CblasTrans)       गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa_array[i] == CblasConjNoTrans) गणः_transb = 0;
      यदि (transa_array[i] == CblasConjTrans)   गणः_transb = 1;
#अन्यथा
      यदि (transa_array[i] == CblasConjNoTrans) गणः_transb = 2;
      यदि (transa_array[i] == CblasConjTrans)   गणः_transb = 3;
#यदि_अन्तः

#अन्यथा
      गणः_m = m_array[i];
      गणः_n = n_array[i];
      गणः_k = k_array[i];

      गणः_lda = lda_array[i];
      गणः_ldb = ldb_array[i];
      गणः_ldc = ldc_array[i];
      
      यदि (transa_array[i] == 'N')   गणः_transa = 0;
      यदि (transa_array[i] == 'T')   गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa_array[i] == 'C')   गणः_transa = 1;
#अन्यथा
      यदि (transa_array[i] == 'C')   गणः_transa = 3;
#यदि_अन्तः
      यदि (transb_array[i] == 'N')   गणः_transb = 0;
      यदि (transb_array[i] == 'T')   गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transb_array[i] == 'C')   गणः_transb = 1;
#अन्यथा
      यदि (transb_array[i] == 'C')   गणः_transb = 3;
#यदि_अन्तः
#यदि_अन्तः

      गणः_nrowa = गणः_m;
      यदि (गणः_transa & 1) गणः_nrowa = गणः_k;
      गणः_nrowb = गणः_k;
      यदि (गणः_transb & 1) गणः_nrowb = गणः_n;

      info=-1;
	
      यदि (गणः_ldc < गणः_m) info = 13;
      यदि (गणः_ldb < गणः_nrowb)  info = 10;
      यदि (गणः_lda < गणः_nrowa)  info =  8;
      यदि (गणः_k < 0)        info =  5;
      यदि (गणः_n < 0)        info =  4;
      यदि (गणः_m < 0)        info =  3;
      यदि (गणः_transb < 0)        info =  2;
      यदि (गणः_transa < 0)        info =  1;      
#यदि defined(CBLAS)
    }
#यदि_अन्तः

    यदि (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
      मुक्त(तर्काः_array);
      प्रत्यागम;
    }

    यदि (गणः_m == 0 || गणः_n == 0) अनुवर्त;

    गणः_बहुलकम्=बहुलकम्;

#यदि defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुण) गणः_m * (द्विगुण) गणः_n * (द्विगुण) गणः_k;
#यदि_अन्तः

#यदि_परिभाषितम् SMALL_MATRIX_OPT
#यदि !defined(COMPLEX)
    यदि (*(FLOAT *)(गणः_alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(गणः_transa, गणः_transb, गणः_m, गणः_n, गणः_k,
                                 *(FLOAT *)(गणः_alpha), *(FLOAT *)(गणः_beta))){
#अन्यथा
    यदि (((FLOAT *)(गणः_alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(गणः_transa, गणः_transb, गणः_m, गणः_n, गणः_k,
                                 ((FLOAT *)(गणः_alpha))[0], ((FLOAT *)(गणः_alpha))[1],
                                 ((FLOAT *)(गणः_beta))[0], ((FLOAT *)(गणः_beta))[1])){
#यदि_अन्तः
      गणः_routine=रिक्त;
#यदि !defined(COMPLEX)
      यदि(*(FLOAT *)(गणः_beta) == 0.0){
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_B0_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((गणः_transb<<2)|गणः_transa));
      }अन्यथा{
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((गणः_transb<<2)|गणः_transa));
      }
#अन्यथा
      यदि(((FLOAT *)(गणः_beta))[0] == 0.0 && ((FLOAT *)(गणः_beta))[1] == 0.0){
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_B0_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((गणः_transb<<2)|गणः_transa));
      }अन्यथा{
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((गणः_transb<<2)|गणः_transa));
      }

#यदि_अन्तः
      
    }अन्यथा{
#यदि_अन्तः
      गणः_routine=(शून्य*)(gemm[(गणः_transb<<2)|गणः_transa]);
#यदि_परिभाषितम् SMALL_MATRIX_OPT
    }
#यदि_अन्तः    

    
    प्रति(j=0; j<गणः_आकार[i]; j++){
      तर्काः_array[गणनम्].m=गणः_m;
      तर्काः_array[गणनम्].n=गणः_n;
      तर्काः_array[गणनम्].k=गणः_k;
      तर्काः_array[गणनम्].lda=गणः_lda;
      तर्काः_array[गणनम्].ldb=गणः_ldb;
      तर्काः_array[गणनम्].ldc=गणः_ldc;
      तर्काः_array[गणनम्].alpha=गणः_alpha;
      तर्काः_array[गणनम्].beta=गणः_beta;

#यदि defined(CBLAS)
      यदि (क्रमः == CblasRowMajor) {
	तर्काः_array[गणनम्].a=(b_array[matrix_idx+j]);
	तर्काः_array[गणनम्].b=(a_array[matrix_idx+j]);
      } अन्यथा
#यदि_अन्तः
      {
	तर्काः_array[गणनम्].a=(a_array[matrix_idx+j]);
	तर्काः_array[गणनम्].b=(b_array[matrix_idx+j]);
      }
      
      तर्काः_array[गणनम्].c=(c_array[matrix_idx+j]);
      
      तर्काः_array[गणनम्].routine_बहुलकम्=गणः_बहुलकम्;
      तर्काः_array[गणनम्].routine=गणः_routine;
#यदि_परिभाषितम् SMALL_MATRIX_OPT
      यदि (!गणः_routine)
        तर्काः_array[गणनम्].routine=गणः_small_matrix_opt_routine;
#यदि_अन्तः      
      गणनम्++;
    }
  }

  यदि(गणनम्>0){
    GEMM_BATCH_THREAD(तर्काः_array,गणनम्);
  }

  मुक्त(तर्काः_array);
}
