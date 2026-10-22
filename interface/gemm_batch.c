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

فراغ openblas_warning(صحيح verbose, ثابت حرف * پیغام);

#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD dgemm_batch_thread
#وإلا
#تعریف_کرو ERROR_NAME "SGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD sgemm_batch_thread
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD zgemm_batch_thread
#وإلا
#تعریف_کرو ERROR_NAME "CGEMM_BATCH "
#تعریف_کرو GEMM_BATCH_THREAD cgemm_batch_thread
#جے_انت
#جے_انت
ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#إذا defined(SMALL_میٹرکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#تعریف_کرو USE_SMALL_میٹرکس_OPT 1
#وإلا
#تعریف_کرو USE_SMALL_میٹرکس_OPT 0
#جے_انت

#إذا USE_SMALL_میٹرکس_OPT
#جے_بنا_تعریف DYNAMIC_ARCH
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#تعریف_کرو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#جے_انت


#جے_بنا_تعریف COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کرو GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کرو GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#وإلا

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

#تعریف_کرو ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کرو ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#جے_انت
#جے_انت

OPENBLAS_EXPORT

#جے_بنا_تعریف CBLAS
فراغ NAME(حرف *transa_array, حرف *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * گروہ_ناپ) {
       blasint گروہ_گنتی = *gcount;
#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE *  transa_array, تعداد CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#جے_بنا_تعریف COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint گروہ_گنتی, blasint * گروہ_ناپ) {
#وإلا
	   فراغ * valpha_array,
	   فراغ ** va_array, blasint * lda_array,
	   فراغ ** vb_array, blasint * ldb_array,
	   فراغ * vbeta_array,
	   فراغ ** vc_array, blasint * ldc_array, blasint گروہ_گنتی, blasint * گروہ_ناپ) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#جے_انت
#جے_انت
  BLASLONG گروہ_m, گروہ_n, گروہ_k;
  BLASLONG گروہ_lda, گروہ_ldb, گروہ_ldc;

  blas_arg_t * دلیلاں_array=عدم;

  صحيح منوال=0, گروہ_منوال=0;
  blasint کل_عدد=0;

  blasint i=0, j=0, میٹرکس_idx=0, گنتی=0;

  صحيح گروہ_transa, گروہ_transb;
  BLASLONG گروہ_nrowa, گروہ_nrowb;
  blasint info;

  فراغ * گروہ_alpha, * گروہ_beta;
  فراغ * گروہ_routine=عدم;
#جے_تعریف_شدہ SMALL_میٹرکس_OPT
  فراغ * گروہ_small_میٹرکس_opt_routine=عدم;
#جے_انت
  
#إذا defined (SMP) || defined(SMALL_میٹرکس_OPT)
  مزدوج MNK;
#جے_انت

  PRINT_DEBUG_CNAME;
    
  لكل(i=0; i<گروہ_گنتی; i++){
    کل_عدد+=گروہ_ناپ[i];
  }

  دلیلاں_array=(blas_arg_t *)حافظہ_دیو(کل_عدد * حجم(blas_arg_t));
  
  إذا(دلیلاں_array == عدم){
    openblas_warning(0, "memory alloc failed!\n");
    أرجع;
  }

#جے_تعریف_شدہ SMP
#جے_بنا_تعریف COMPLEX
#جے_تعریف_شدہ XDOUBLE
  منوال  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  منوال  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  منوال  =  BLAS_SINGLE  | BLAS_REAL;
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
  منوال  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  منوال  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  منوال  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جے_انت
#جے_انت
#جے_انت
  
  لكل(i=0; i<گروہ_گنتی; میٹرکس_idx+=گروہ_ناپ[i], i++){
    گروہ_alpha = (فراغ *)&alpha_array[i * COMPSIZE];
    گروہ_beta = (فراغ *)&beta_array[i * COMPSIZE];

    گروہ_m = گروہ_n = گروہ_k = 0;
    گروہ_lda = گروہ_ldb = گروہ_ldc = 0;
    گروہ_transa = -1;
    گروہ_transb = -1;
    info   = 0;
   
#إذا defined(CBLAS) 
    إذا (ترتیب == Cblas_خانہ_Major) {
      گروہ_m = m_array[i];
      گروہ_n = n_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = lda_array[i];
      گروہ_ldb = ldb_array[i];
      گروہ_ldc = ldc_array[i];
      
      إذا (transa_array[i] == CblasNoTrans)     گروہ_transa = 0;
      إذا (transa_array[i] == CblasTrans)       گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) گروہ_transa = 0;
      إذا (transa_array[i] == CblasConjTrans)   گروہ_transa = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) گروہ_transa = 2;
      إذا (transa_array[i] == CblasConjTrans)   گروہ_transa = 3;
#جے_انت
      إذا (transb_array[i] == CblasNoTrans)     گروہ_transb = 0;
      إذا (transb_array[i] == CblasTrans)       گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) گروہ_transb = 0;
      إذا (transb_array[i] == CblasConjTrans)   گروہ_transb = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) گروہ_transb = 2;
      إذا (transb_array[i] == CblasConjTrans)   گروہ_transb = 3;
#جے_انت
      گروہ_nrowa = گروہ_m;
      إذا (گروہ_transa & 1) گروہ_nrowa = گروہ_k;
      گروہ_nrowb = گروہ_k;
      إذا (گروہ_transb & 1) گروہ_nrowb = گروہ_n;

      info=-1;
	
      إذا (گروہ_ldc < گروہ_m) info = 13;
      إذا (گروہ_ldb < گروہ_nrowb)  info = 10;
      إذا (گروہ_lda < گروہ_nrowa)  info =  8;
      إذا (گروہ_k < 0)        info =  5;
      إذا (گروہ_n < 0)        info =  4;
      إذا (گروہ_m < 0)        info =  3;
      إذا (گروہ_transb < 0)        info =  2;
      إذا (گروہ_transa < 0)        info =  1;

    }وإلا إذا (ترتیب == Cblas_صف_Major) {
      
      گروہ_m = n_array[i];
      گروہ_n = m_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = ldb_array[i];
      گروہ_ldb = lda_array[i];
      گروہ_ldc = ldc_array[i];

      إذا (transb_array[i] == CblasNoTrans)     گروہ_transa = 0;
      إذا (transb_array[i] == CblasTrans)       گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) گروہ_transa = 0;
      إذا (transb_array[i] == CblasConjTrans)   گروہ_transa = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) گروہ_transa = 2;
      إذا (transb_array[i] == CblasConjTrans)   گروہ_transa = 3;
#جے_انت
      إذا (transa_array[i] == CblasNoTrans)     گروہ_transb = 0;
      إذا (transa_array[i] == CblasTrans)       گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) گروہ_transb = 0;
      إذا (transa_array[i] == CblasConjTrans)   گروہ_transb = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) گروہ_transb = 2;
      إذا (transa_array[i] == CblasConjTrans)   گروہ_transb = 3;
#جے_انت

#وإلا
      گروہ_m = m_array[i];
      گروہ_n = n_array[i];
      گروہ_k = k_array[i];

      گروہ_lda = lda_array[i];
      گروہ_ldb = ldb_array[i];
      گروہ_ldc = ldc_array[i];
      
      إذا (transa_array[i] == 'N')   گروہ_transa = 0;
      إذا (transa_array[i] == 'T')   گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transa_array[i] == 'C')   گروہ_transa = 1;
#وإلا
      إذا (transa_array[i] == 'C')   گروہ_transa = 3;
#جے_انت
      إذا (transb_array[i] == 'N')   گروہ_transb = 0;
      إذا (transb_array[i] == 'T')   گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transb_array[i] == 'C')   گروہ_transb = 1;
#وإلا
      إذا (transb_array[i] == 'C')   گروہ_transb = 3;
#جے_انت
#جے_انت

      گروہ_nrowa = گروہ_m;
      إذا (گروہ_transa & 1) گروہ_nrowa = گروہ_k;
      گروہ_nrowb = گروہ_k;
      إذا (گروہ_transb & 1) گروہ_nrowb = گروہ_n;

      info=-1;
	
      إذا (گروہ_ldc < گروہ_m) info = 13;
      إذا (گروہ_ldb < گروہ_nrowb)  info = 10;
      إذا (گروہ_lda < گروہ_nrowa)  info =  8;
      إذا (گروہ_k < 0)        info =  5;
      إذا (گروہ_n < 0)        info =  4;
      إذا (گروہ_m < 0)        info =  3;
      إذا (گروہ_transb < 0)        info =  2;
      إذا (گروہ_transa < 0)        info =  1;      
#إذا defined(CBLAS)
    }
#جے_انت

    إذا (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      حرر(دلیلاں_array);
      أرجع;
    }

    إذا (گروہ_m == 0 || گروہ_n == 0) استمر;

    گروہ_منوال=منوال;

#إذا defined(SMP) || defined(SMALL_میٹرکس_OPT)
    MNK = (مزدوج) گروہ_m * (مزدوج) گروہ_n * (مزدوج) گروہ_k;
#جے_انت

#جے_تعریف_شدہ SMALL_میٹرکس_OPT
#إذا !defined(COMPLEX)
    إذا (*(FLOAT *)(گروہ_alpha) != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 *(FLOAT *)(گروہ_alpha), *(FLOAT *)(گروہ_beta))){
#وإلا
    إذا (((FLOAT *)(گروہ_alpha))[0] != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 ((FLOAT *)(گروہ_alpha))[0], ((FLOAT *)(گروہ_alpha))[1],
                                 ((FLOAT *)(گروہ_beta))[0], ((FLOAT *)(گروہ_beta))[1])){
#جے_انت
      گروہ_routine=عدم;
#إذا !defined(COMPLEX)
      إذا(*(FLOAT *)(گروہ_beta) == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }وإلا{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }
#وإلا
      إذا(((FLOAT *)(گروہ_beta))[0] == 0.0 && ((FLOAT *)(گروہ_beta))[1] == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }وإلا{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }

#جے_انت
      
    }وإلا{
#جے_انت
      گروہ_routine=(فراغ*)(gemm[(گروہ_transb<<2)|گروہ_transa]);
#جے_تعریف_شدہ SMALL_میٹرکس_OPT
    }
#جے_انت    

    
    لكل(j=0; j<گروہ_ناپ[i]; j++){
      دلیلاں_array[گنتی].m=گروہ_m;
      دلیلاں_array[گنتی].n=گروہ_n;
      دلیلاں_array[گنتی].k=گروہ_k;
      دلیلاں_array[گنتی].lda=گروہ_lda;
      دلیلاں_array[گنتی].ldb=گروہ_ldb;
      دلیلاں_array[گنتی].ldc=گروہ_ldc;
      دلیلاں_array[گنتی].alpha=گروہ_alpha;
      دلیلاں_array[گنتی].beta=گروہ_beta;

#إذا defined(CBLAS)
      إذا (ترتیب == Cblas_صف_Major) {
	دلیلاں_array[گنتی].a=(b_array[میٹرکس_idx+j]);
	دلیلاں_array[گنتی].b=(a_array[میٹرکس_idx+j]);
      } وإلا
#جے_انت
      {
	دلیلاں_array[گنتی].a=(a_array[میٹرکس_idx+j]);
	دلیلاں_array[گنتی].b=(b_array[میٹرکس_idx+j]);
      }
      
      دلیلاں_array[گنتی].c=(c_array[میٹرکس_idx+j]);
      
      دلیلاں_array[گنتی].routine_منوال=گروہ_منوال;
      دلیلاں_array[گنتی].routine=گروہ_routine;
#جے_تعریف_شدہ SMALL_میٹرکس_OPT
      إذا (!گروہ_routine)
        دلیلاں_array[گنتی].routine=گروہ_small_میٹرکس_opt_routine;
#جے_انت      
      گنتی++;
    }
  }

  إذا(گنتی>0){
    GEMM_BATCH_THREAD(دلیلاں_array,گنتی);
  }

  حرر(دلیلاں_array);
}
