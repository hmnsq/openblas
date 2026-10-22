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

فراغ openblas_warning(صحيح verbose, ثابت حرف * پیام);

#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "QGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "DGEMM_BATCH "
#تعریف_کن GEMM_BATCH_THREAD dgemm_batch_thread
#وإلا
#تعریف_کن ERROR_NAME "SGEMM_BATCH "
#تعریف_کن GEMM_BATCH_THREAD sgemm_batch_thread
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
#تعریف_کن ERROR_NAME "XGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کن ERROR_NAME "ZGEMM_BATCH "
#تعریف_کن GEMM_BATCH_THREAD zgemm_batch_thread
#وإلا
#تعریف_کن ERROR_NAME "CGEMM_BATCH "
#تعریف_کن GEMM_BATCH_THREAD cgemm_batch_thread
#پایان_اگر
#پایان_اگر
ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#إذا defined(SMALL_ماتریس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#تعریف_کن USE_SMALL_ماتریس_OPT 1
#وإلا
#تعریف_کن USE_SMALL_ماتریس_OPT 0
#پایان_اگر

#إذا USE_SMALL_ماتریس_OPT
#اگر_تعریف_نشده DYNAMIC_ARCH
#تعریف_کن SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#تعریف_کن SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#پایان_اگر


#اگر_تعریف_نشده COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کن GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کن GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
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

#تعریف_کن ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کن ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#پایان_اگر
#پایان_اگر

OPENBLAS_EXPORT

#اگر_تعریف_نشده CBLAS
فراغ NAME(حرف *transa_array, حرف *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * گروه_اندازه) {
       blasint گروه_شمار = *gcount;
#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE *  transa_array, تعداد CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#اگر_تعریف_نشده COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint گروه_شمار, blasint * گروه_اندازه) {
#وإلا
	   فراغ * valpha_array,
	   فراغ ** va_array, blasint * lda_array,
	   فراغ ** vb_array, blasint * ldb_array,
	   فراغ * vbeta_array,
	   فراغ ** vc_array, blasint * ldc_array, blasint گروه_شمار, blasint * گروه_اندازه) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#پایان_اگر
#پایان_اگر
  BLASLONG گروه_m, گروه_n, گروه_k;
  BLASLONG گروه_lda, گروه_ldb, گروه_ldc;

  blas_arg_t * آرگومانها_array=عدم;

  صحيح نما=0, گروه_نما=0;
  blasint مجموع_رقم=0;

  blasint i=0, j=0, ماتریس_idx=0, شمار=0;

  صحيح گروه_transa, گروه_transb;
  BLASLONG گروه_nrowa, گروه_nrowb;
  blasint info;

  فراغ * گروه_alpha, * گروه_beta;
  فراغ * گروه_routine=عدم;
#اگر_تعریف_شده SMALL_ماتریس_OPT
  فراغ * گروه_small_ماتریس_opt_routine=عدم;
#پایان_اگر
  
#إذا defined (SMP) || defined(SMALL_ماتریس_OPT)
  مزدوج MNK;
#پایان_اگر

  PRINT_DEBUG_CNAME;
    
  لكل(i=0; i<گروه_شمار; i++){
    مجموع_رقم+=گروه_اندازه[i];
  }

  آرگومانها_array=(blas_arg_t *)حافظه_بده(مجموع_رقم * حجم(blas_arg_t));
  
  إذا(آرگومانها_array == عدم){
    openblas_warning(0, "memory alloc failed!\n");
    أرجع;
  }

#اگر_تعریف_شده SMP
#اگر_تعریف_نشده COMPLEX
#اگر_تعریف_شده XDOUBLE
  نما  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  نما  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  نما  =  BLAS_SINGLE  | BLAS_REAL;
#پایان_اگر
#وإلا
#اگر_تعریف_شده XDOUBLE
  نما  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  نما  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  نما  =  BLAS_SINGLE  | BLAS_COMPLEX;
#پایان_اگر
#پایان_اگر
#پایان_اگر
  
  لكل(i=0; i<گروه_شمار; ماتریس_idx+=گروه_اندازه[i], i++){
    گروه_alpha = (فراغ *)&alpha_array[i * COMPSIZE];
    گروه_beta = (فراغ *)&beta_array[i * COMPSIZE];

    گروه_m = گروه_n = گروه_k = 0;
    گروه_lda = گروه_ldb = گروه_ldc = 0;
    گروه_transa = -1;
    گروه_transb = -1;
    info   = 0;
   
#إذا defined(CBLAS) 
    إذا (ترتیب == Cblas_ستون_Major) {
      گروه_m = m_array[i];
      گروه_n = n_array[i];
      گروه_k = k_array[i];

      گروه_lda = lda_array[i];
      گروه_ldb = ldb_array[i];
      گروه_ldc = ldc_array[i];
      
      إذا (transa_array[i] == CblasNoTrans)     گروه_transa = 0;
      إذا (transa_array[i] == CblasTrans)       گروه_transa = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) گروه_transa = 0;
      إذا (transa_array[i] == CblasConjTrans)   گروه_transa = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) گروه_transa = 2;
      إذا (transa_array[i] == CblasConjTrans)   گروه_transa = 3;
#پایان_اگر
      إذا (transb_array[i] == CblasNoTrans)     گروه_transb = 0;
      إذا (transb_array[i] == CblasTrans)       گروه_transb = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) گروه_transb = 0;
      إذا (transb_array[i] == CblasConjTrans)   گروه_transb = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) گروه_transb = 2;
      إذا (transb_array[i] == CblasConjTrans)   گروه_transb = 3;
#پایان_اگر
      گروه_nrowa = گروه_m;
      إذا (گروه_transa & 1) گروه_nrowa = گروه_k;
      گروه_nrowb = گروه_k;
      إذا (گروه_transb & 1) گروه_nrowb = گروه_n;

      info=-1;
	
      إذا (گروه_ldc < گروه_m) info = 13;
      إذا (گروه_ldb < گروه_nrowb)  info = 10;
      إذا (گروه_lda < گروه_nrowa)  info =  8;
      إذا (گروه_k < 0)        info =  5;
      إذا (گروه_n < 0)        info =  4;
      إذا (گروه_m < 0)        info =  3;
      إذا (گروه_transb < 0)        info =  2;
      إذا (گروه_transa < 0)        info =  1;

    }وإلا إذا (ترتیب == Cblas_ردیف_Major) {
      
      گروه_m = n_array[i];
      گروه_n = m_array[i];
      گروه_k = k_array[i];

      گروه_lda = ldb_array[i];
      گروه_ldb = lda_array[i];
      گروه_ldc = ldc_array[i];

      إذا (transb_array[i] == CblasNoTrans)     گروه_transa = 0;
      إذا (transb_array[i] == CblasTrans)       گروه_transa = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) گروه_transa = 0;
      إذا (transb_array[i] == CblasConjTrans)   گروه_transa = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) گروه_transa = 2;
      إذا (transb_array[i] == CblasConjTrans)   گروه_transa = 3;
#پایان_اگر
      إذا (transa_array[i] == CblasNoTrans)     گروه_transb = 0;
      إذا (transa_array[i] == CblasTrans)       گروه_transb = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) گروه_transb = 0;
      إذا (transa_array[i] == CblasConjTrans)   گروه_transb = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) گروه_transb = 2;
      إذا (transa_array[i] == CblasConjTrans)   گروه_transb = 3;
#پایان_اگر

#وإلا
      گروه_m = m_array[i];
      گروه_n = n_array[i];
      گروه_k = k_array[i];

      گروه_lda = lda_array[i];
      گروه_ldb = ldb_array[i];
      گروه_ldc = ldc_array[i];
      
      إذا (transa_array[i] == 'N')   گروه_transa = 0;
      إذا (transa_array[i] == 'T')   گروه_transa = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transa_array[i] == 'C')   گروه_transa = 1;
#وإلا
      إذا (transa_array[i] == 'C')   گروه_transa = 3;
#پایان_اگر
      إذا (transb_array[i] == 'N')   گروه_transb = 0;
      إذا (transb_array[i] == 'T')   گروه_transb = 1;
#اگر_تعریف_نشده COMPLEX
      إذا (transb_array[i] == 'C')   گروه_transb = 1;
#وإلا
      إذا (transb_array[i] == 'C')   گروه_transb = 3;
#پایان_اگر
#پایان_اگر

      گروه_nrowa = گروه_m;
      إذا (گروه_transa & 1) گروه_nrowa = گروه_k;
      گروه_nrowb = گروه_k;
      إذا (گروه_transb & 1) گروه_nrowb = گروه_n;

      info=-1;
	
      إذا (گروه_ldc < گروه_m) info = 13;
      إذا (گروه_ldb < گروه_nrowb)  info = 10;
      إذا (گروه_lda < گروه_nrowa)  info =  8;
      إذا (گروه_k < 0)        info =  5;
      إذا (گروه_n < 0)        info =  4;
      إذا (گروه_m < 0)        info =  3;
      إذا (گروه_transb < 0)        info =  2;
      إذا (گروه_transa < 0)        info =  1;      
#إذا defined(CBLAS)
    }
#پایان_اگر

    إذا (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      حرر(آرگومانها_array);
      أرجع;
    }

    إذا (گروه_m == 0 || گروه_n == 0) استمر;

    گروه_نما=نما;

#إذا defined(SMP) || defined(SMALL_ماتریس_OPT)
    MNK = (مزدوج) گروه_m * (مزدوج) گروه_n * (مزدوج) گروه_k;
#پایان_اگر

#اگر_تعریف_شده SMALL_ماتریس_OPT
#إذا !defined(COMPLEX)
    إذا (*(FLOAT *)(گروه_alpha) != ZERO &&
        GEMM_SMALL_ماتریس_PERMIT(گروه_transa, گروه_transb, گروه_m, گروه_n, گروه_k,
                                 *(FLOAT *)(گروه_alpha), *(FLOAT *)(گروه_beta))){
#وإلا
    إذا (((FLOAT *)(گروه_alpha))[0] != ZERO &&
        GEMM_SMALL_ماتریس_PERMIT(گروه_transa, گروه_transb, گروه_m, گروه_n, گروه_k,
                                 ((FLOAT *)(گروه_alpha))[0], ((FLOAT *)(گروه_alpha))[1],
                                 ((FLOAT *)(گروه_beta))[0], ((FLOAT *)(گروه_beta))[1])){
#پایان_اگر
      گروه_routine=عدم;
#إذا !defined(COMPLEX)
      إذا(*(FLOAT *)(گروه_beta) == 0.0){
	گروه_نما=نما | BLAS_SMALL_B0_OPT;
	گروه_small_ماتریس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروه_transb<<2)|گروه_transa));
      }وإلا{
	گروه_نما=نما | BLAS_SMALL_OPT;
	گروه_small_ماتریس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروه_transb<<2)|گروه_transa));
      }
#وإلا
      إذا(((FLOAT *)(گروه_beta))[0] == 0.0 && ((FLOAT *)(گروه_beta))[1] == 0.0){
	گروه_نما=نما | BLAS_SMALL_B0_OPT;
	گروه_small_ماتریس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((گروه_transb<<2)|گروه_transa));
      }وإلا{
	گروه_نما=نما | BLAS_SMALL_OPT;
	گروه_small_ماتریس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((گروه_transb<<2)|گروه_transa));
      }

#پایان_اگر
      
    }وإلا{
#پایان_اگر
      گروه_routine=(فراغ*)(gemm[(گروه_transb<<2)|گروه_transa]);
#اگر_تعریف_شده SMALL_ماتریس_OPT
    }
#پایان_اگر    

    
    لكل(j=0; j<گروه_اندازه[i]; j++){
      آرگومانها_array[شمار].m=گروه_m;
      آرگومانها_array[شمار].n=گروه_n;
      آرگومانها_array[شمار].k=گروه_k;
      آرگومانها_array[شمار].lda=گروه_lda;
      آرگومانها_array[شمار].ldb=گروه_ldb;
      آرگومانها_array[شمار].ldc=گروه_ldc;
      آرگومانها_array[شمار].alpha=گروه_alpha;
      آرگومانها_array[شمار].beta=گروه_beta;

#إذا defined(CBLAS)
      إذا (ترتیب == Cblas_ردیف_Major) {
	آرگومانها_array[شمار].a=(b_array[ماتریس_idx+j]);
	آرگومانها_array[شمار].b=(a_array[ماتریس_idx+j]);
      } وإلا
#پایان_اگر
      {
	آرگومانها_array[شمار].a=(a_array[ماتریس_idx+j]);
	آرگومانها_array[شمار].b=(b_array[ماتریس_idx+j]);
      }
      
      آرگومانها_array[شمار].c=(c_array[ماتریس_idx+j]);
      
      آرگومانها_array[شمار].routine_نما=گروه_نما;
      آرگومانها_array[شمار].routine=گروه_routine;
#اگر_تعریف_شده SMALL_ماتریس_OPT
      إذا (!گروه_routine)
        آرگومانها_array[شمار].routine=گروه_small_ماتریس_opt_routine;
#پایان_اگر      
      شمار++;
    }
  }

  إذا(شمار>0){
    GEMM_BATCH_THREAD(آرگومانها_array,شمار);
  }

  حرر(آرگومانها_array);
}
