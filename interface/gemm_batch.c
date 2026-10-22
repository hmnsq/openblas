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

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل "common.h"

خالی openblas_warning(صحیح verbose, ثابت حرف * پیغام);

#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMM_BATCH "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD dgemm_batch_thread
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD sgemm_batch_thread
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM_BATCH "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD zgemm_batch_thread
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD cgemm_batch_thread
#اگر_ختم
#اگر_ختم
ساکن صحیح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#اگر defined(SMALL_میٹرکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#تعریف_کرو USE_SMALL_میٹرکس_OPT 1
#ورنہ
#تعریف_کرو USE_SMALL_میٹرکس_OPT 0
#اگر_ختم

#اگر USE_SMALL_میٹرکس_OPT
#اگر_غیر_تعریف_شدہ DYNAMIC_ARCH
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((خالی *)(جدول[idx]))
#ورنہ
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((خالی *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#اگر_ختم


#اگر_غیر_تعریف_شدہ COMPLEX
ساکن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساکن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کرو GEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کرو GEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ورنہ

ساکن size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ساکن size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#تعریف_کرو ZGEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کرو ZGEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#اگر_ختم
#اگر_ختم

OPENBLAS_EXPORT

#اگر_غیر_تعریف_شدہ CBLAS
خالی NAME(حرف *transa_array, حرف *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * گروہ_سائز) {
       blasint گروہ_گنتی = *gcount;
#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_TRANSPOSE *  transa_array, شمار CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint گروہ_گنتی, blasint * گروہ_سائز) {
#ورنہ
	   خالی * valpha_array,
	   خالی ** va_array, blasint * lda_array,
	   خالی ** vb_array, blasint * ldb_array,
	   خالی * vbeta_array,
	   خالی ** vc_array, blasint * ldc_array, blasint گروہ_گنتی, blasint * گروہ_سائز) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#اگر_ختم
#اگر_ختم
  BLASLONG گروہ_m, گروہ_n, گروہ_k;
  BLASLONG گروہ_lda, گروہ_ldb, گروہ_ldc;

  blas_arg_t * دلائل_array=تہی;

  صحیح منوال=0, گروہ_منوال=0;
  blasint کل_عدد=0;

  blasint i=0, j=0, میٹرکس_idx=0, گنتی=0;

  صحیح گروہ_transa, گروہ_transb;
  BLASLONG گروہ_nrowa, گروہ_nrowb;
  blasint info;

  خالی * گروہ_alpha, * گروہ_beta;
  خالی * گروہ_routine=تہی;
#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
  خالی * گروہ_small_میٹرکس_opt_routine=تہی;
#اگر_ختم
  
#اگر defined (SMP) || defined(SMALL_میٹرکس_OPT)
  دوہرا MNK;
#اگر_ختم

  PRINT_DEBUG_CNAME;
    
  برائے(i=0; i<گروہ_گنتی; i++){
    کل_عدد+=گروہ_سائز[i];
  }

  دلائل_array=(blas_arg_t *)حافظہ_دو(کل_عدد * حجم(blas_arg_t));
  
  اگر(دلائل_array == تہی){
    openblas_warning(0, "memory alloc failed!\n");
    واپس;
  }

#اگر_تعریف_شدہ SMP
#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
  منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنہ_اگر defined(DOUBLE)
  منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#ورنہ
  منوال  =  BLAS_SINGLE  | BLAS_REAL;
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
  منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنہ_اگر defined(DOUBLE)
  منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ورنہ
  منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#اگر_ختم
#اگر_ختم
#اگر_ختم
  
  برائے(i=0; i<گروہ_گنتی; میٹرکس_idx+=گروہ_سائز[i], i++){
    گروہ_alpha = (خالی *)&alpha_array[i * COMPSIZE];
    گروہ_beta = (خالی *)&beta_array[i * COMPSIZE];

    گروہ_m = گروہ_n = گروہ_k = 0;
    گروہ_lda = گروہ_ldb = گروہ_ldc = 0;
    گروہ_transa = -1;
    گروہ_transb = -1;
    info   = 0;
   
#اگر defined(CBLAS) 
    اگر (ترتیب == Cblas_ستون_Major) {
      گروہ_m = m_array[i];
      گروہ_n = n_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = lda_array[i];
      گروہ_ldb = ldb_array[i];
      گروہ_ldc = ldc_array[i];
      
      اگر (transa_array[i] == CblasNoTrans)     گروہ_transa = 0;
      اگر (transa_array[i] == CblasTrans)       گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa_array[i] == CblasConjNoTrans) گروہ_transa = 0;
      اگر (transa_array[i] == CblasConjTrans)   گروہ_transa = 1;
#ورنہ
      اگر (transa_array[i] == CblasConjNoTrans) گروہ_transa = 2;
      اگر (transa_array[i] == CblasConjTrans)   گروہ_transa = 3;
#اگر_ختم
      اگر (transb_array[i] == CblasNoTrans)     گروہ_transb = 0;
      اگر (transb_array[i] == CblasTrans)       گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transb_array[i] == CblasConjNoTrans) گروہ_transb = 0;
      اگر (transb_array[i] == CblasConjTrans)   گروہ_transb = 1;
#ورنہ
      اگر (transb_array[i] == CblasConjNoTrans) گروہ_transb = 2;
      اگر (transb_array[i] == CblasConjTrans)   گروہ_transb = 3;
#اگر_ختم
      گروہ_nrowa = گروہ_m;
      اگر (گروہ_transa & 1) گروہ_nrowa = گروہ_k;
      گروہ_nrowb = گروہ_k;
      اگر (گروہ_transb & 1) گروہ_nrowb = گروہ_n;

      info=-1;
	
      اگر (گروہ_ldc < گروہ_m) info = 13;
      اگر (گروہ_ldb < گروہ_nrowb)  info = 10;
      اگر (گروہ_lda < گروہ_nrowa)  info =  8;
      اگر (گروہ_k < 0)        info =  5;
      اگر (گروہ_n < 0)        info =  4;
      اگر (گروہ_m < 0)        info =  3;
      اگر (گروہ_transb < 0)        info =  2;
      اگر (گروہ_transa < 0)        info =  1;

    }ورنہ اگر (ترتیب == Cblas_صف_Major) {
      
      گروہ_m = n_array[i];
      گروہ_n = m_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = ldb_array[i];
      گروہ_ldb = lda_array[i];
      گروہ_ldc = ldc_array[i];

      اگر (transb_array[i] == CblasNoTrans)     گروہ_transa = 0;
      اگر (transb_array[i] == CblasTrans)       گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transb_array[i] == CblasConjNoTrans) گروہ_transa = 0;
      اگر (transb_array[i] == CblasConjTrans)   گروہ_transa = 1;
#ورنہ
      اگر (transb_array[i] == CblasConjNoTrans) گروہ_transa = 2;
      اگر (transb_array[i] == CblasConjTrans)   گروہ_transa = 3;
#اگر_ختم
      اگر (transa_array[i] == CblasNoTrans)     گروہ_transb = 0;
      اگر (transa_array[i] == CblasTrans)       گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa_array[i] == CblasConjNoTrans) گروہ_transb = 0;
      اگر (transa_array[i] == CblasConjTrans)   گروہ_transb = 1;
#ورنہ
      اگر (transa_array[i] == CblasConjNoTrans) گروہ_transb = 2;
      اگر (transa_array[i] == CblasConjTrans)   گروہ_transb = 3;
#اگر_ختم

#ورنہ
      گروہ_m = m_array[i];
      گروہ_n = n_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = lda_array[i];
      گروہ_ldb = ldb_array[i];
      گروہ_ldc = ldc_array[i];
      
      اگر (transa_array[i] == 'N')   گروہ_transa = 0;
      اگر (transa_array[i] == 'T')   گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa_array[i] == 'C')   گروہ_transa = 1;
#ورنہ
      اگر (transa_array[i] == 'C')   گروہ_transa = 3;
#اگر_ختم
      اگر (transb_array[i] == 'N')   گروہ_transb = 0;
      اگر (transb_array[i] == 'T')   گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transb_array[i] == 'C')   گروہ_transb = 1;
#ورنہ
      اگر (transb_array[i] == 'C')   گروہ_transb = 3;
#اگر_ختم
#اگر_ختم

      گروہ_nrowa = گروہ_m;
      اگر (گروہ_transa & 1) گروہ_nrowa = گروہ_k;
      گروہ_nrowb = گروہ_k;
      اگر (گروہ_transb & 1) گروہ_nrowb = گروہ_n;

      info=-1;
	
      اگر (گروہ_ldc < گروہ_m) info = 13;
      اگر (گروہ_ldb < گروہ_nrowb)  info = 10;
      اگر (گروہ_lda < گروہ_nrowa)  info =  8;
      اگر (گروہ_k < 0)        info =  5;
      اگر (گروہ_n < 0)        info =  4;
      اگر (گروہ_m < 0)        info =  3;
      اگر (گروہ_transb < 0)        info =  2;
      اگر (گروہ_transa < 0)        info =  1;      
#اگر defined(CBLAS)
    }
#اگر_ختم

    اگر (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      آزاد(دلائل_array);
      واپس;
    }

    اگر (گروہ_m == 0 || گروہ_n == 0) جاری;

    گروہ_منوال=منوال;

#اگر defined(SMP) || defined(SMALL_میٹرکس_OPT)
    MNK = (دوہرا) گروہ_m * (دوہرا) گروہ_n * (دوہرا) گروہ_k;
#اگر_ختم

#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
#اگر !defined(COMPLEX)
    اگر (*(FLOAT *)(گروہ_alpha) != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 *(FLOAT *)(گروہ_alpha), *(FLOAT *)(گروہ_beta))){
#ورنہ
    اگر (((FLOAT *)(گروہ_alpha))[0] != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 ((FLOAT *)(گروہ_alpha))[0], ((FLOAT *)(گروہ_alpha))[1],
                                 ((FLOAT *)(گروہ_beta))[0], ((FLOAT *)(گروہ_beta))[1])){
#اگر_ختم
      گروہ_routine=تہی;
#اگر !defined(COMPLEX)
      اگر(*(FLOAT *)(گروہ_beta) == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }ورنہ{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }
#ورنہ
      اگر(((FLOAT *)(گروہ_beta))[0] == 0.0 && ((FLOAT *)(گروہ_beta))[1] == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }ورنہ{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }

#اگر_ختم
      
    }ورنہ{
#اگر_ختم
      گروہ_routine=(خالی*)(gemm[(گروہ_transb<<2)|گروہ_transa]);
#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
    }
#اگر_ختم    

    
    برائے(j=0; j<گروہ_سائز[i]; j++){
      دلائل_array[گنتی].m=گروہ_m;
      دلائل_array[گنتی].n=گروہ_n;
      دلائل_array[گنتی].k=گروہ_k;
      دلائل_array[گنتی].lda=گروہ_lda;
      دلائل_array[گنتی].ldb=گروہ_ldb;
      دلائل_array[گنتی].ldc=گروہ_ldc;
      دلائل_array[گنتی].alpha=گروہ_alpha;
      دلائل_array[گنتی].beta=گروہ_beta;

#اگر defined(CBLAS)
      اگر (ترتیب == Cblas_صف_Major) {
	دلائل_array[گنتی].a=(b_array[میٹرکس_idx+j]);
	دلائل_array[گنتی].b=(a_array[میٹرکس_idx+j]);
      } ورنہ
#اگر_ختم
      {
	دلائل_array[گنتی].a=(a_array[میٹرکس_idx+j]);
	دلائل_array[گنتی].b=(b_array[میٹرکس_idx+j]);
      }
      
      دلائل_array[گنتی].c=(c_array[میٹرکس_idx+j]);
      
      دلائل_array[گنتی].routine_منوال=گروہ_منوال;
      دلائل_array[گنتی].routine=گروہ_routine;
#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
      اگر (!گروہ_routine)
        دلائل_array[گنتی].routine=گروہ_small_میٹرکس_opt_routine;
#اگر_ختم      
      گنتی++;
    }
  }

  اگر(گنتی>0){
    GEMM_BATCH_THREAD(دلائل_array,گنتی);
  }

  آزاد(دلائل_array);
}
