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

#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "QGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "DGEMM_BATCH "
#تعریف_کړه GEMM_BATCH_THREAD dgemm_batch_thread
#وإلا
#تعریف_کړه ERROR_NAME "SGEMM_BATCH "
#تعریف_کړه GEMM_BATCH_THREAD sgemm_batch_thread
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
#تعریف_کړه ERROR_NAME "XGEMM_BATCH "
#وإلاإذا defined(DOUBLE)
#تعریف_کړه ERROR_NAME "ZGEMM_BATCH "
#تعریف_کړه GEMM_BATCH_THREAD zgemm_batch_thread
#وإلا
#تعریف_کړه ERROR_NAME "CGEMM_BATCH "
#تعریف_کړه GEMM_BATCH_THREAD cgemm_batch_thread
#د_که_پای
#د_که_پای
ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#إذا defined(SMALL_ماتریکس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#تعریف_کړه USE_SMALL_ماتریکس_OPT 1
#وإلا
#تعریف_کړه USE_SMALL_ماتریکس_OPT 0
#د_که_پای

#إذا USE_SMALL_ماتریکس_OPT
#که_نه_تعریف_شوی DYNAMIC_ARCH
#تعریف_کړه SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#تعریف_کړه SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#د_که_پای


#که_نه_تعریف_شوی COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#تعریف_کړه GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#تعریف_کړه GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
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

#تعریف_کړه ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#تعریف_کړه ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#د_که_پای
#د_که_پای

OPENBLAS_EXPORT

#که_نه_تعریف_شوی CBLAS
فراغ NAME(حرف *transa_array, حرف *transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint * gcount, blasint * ډله_اندازه) {
       blasint ډله_شمېر = *gcount;
#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE *  transa_array, تعداد CBLAS_TRANSPOSE * transb_array,
	   blasint * m_array, blasint * n_array, blasint * k_array,
#که_نه_تعریف_شوی COMPLEX
	   FLOAT * alpha_array,
	   IFLOAT ** a_array, blasint * lda_array,
	   IFLOAT ** b_array, blasint * ldb_array,
	   FLOAT * beta_array,
	   FLOAT ** c_array, blasint * ldc_array, blasint ډله_شمېر, blasint * ډله_اندازه) {
#وإلا
	   فراغ * valpha_array,
	   فراغ ** va_array, blasint * lda_array,
	   فراغ ** vb_array, blasint * ldb_array,
	   فراغ * vbeta_array,
	   فراغ ** vc_array, blasint * ldc_array, blasint ډله_شمېر, blasint * ډله_اندازه) {

  FLOAT * alpha_array=(FLOAT *)valpha_array;
  FLOAT * beta_array=(FLOAT *)vbeta_array;
  FLOAT ** a_array=(FLOAT**)va_array;
  FLOAT ** b_array=(FLOAT**)vb_array;
  FLOAT ** c_array=(FLOAT**)vc_array;
#د_که_پای
#د_که_پای
  BLASLONG ډله_m, ډله_n, ډله_k;
  BLASLONG ډله_lda, ډله_ldb, ډله_ldc;

  blas_arg_t * دلیلونه_array=عدم;

  صحيح ډېرپېښ=0, ډله_ډېرپېښ=0;
  blasint ټول_عدد=0;

  blasint i=0, j=0, ماتریکس_idx=0, شمېر=0;

  صحيح ډله_transa, ډله_transb;
  BLASLONG ډله_nrowa, ډله_nrowb;
  blasint info;

  فراغ * ډله_alpha, * ډله_beta;
  فراغ * ډله_routine=عدم;
#که_تعریف_شوی SMALL_ماتریکس_OPT
  فراغ * ډله_small_ماتریکس_opt_routine=عدم;
#د_که_پای
  
#إذا defined (SMP) || defined(SMALL_ماتریکس_OPT)
  مزدوج MNK;
#د_که_پای

  PRINT_DEBUG_CNAME;
    
  لكل(i=0; i<ډله_شمېر; i++){
    ټول_عدد+=ډله_اندازه[i];
  }

  دلیلونه_array=(blas_arg_t *)حافظه_ورکړه(ټول_عدد * حجم(blas_arg_t));
  
  إذا(دلیلونه_array == عدم){
    openblas_warning(0, "memory alloc failed!\n");
    أرجع;
  }

#که_تعریف_شوی SMP
#که_نه_تعریف_شوی COMPLEX
#که_تعریف_شوی XDOUBLE
  ډېرپېښ  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  ډېرپېښ  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  ډېرپېښ  =  BLAS_SINGLE  | BLAS_REAL;
#د_که_پای
#وإلا
#که_تعریف_شوی XDOUBLE
  ډېرپېښ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  ډېرپېښ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  ډېرپېښ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#د_که_پای
#د_که_پای
#د_که_پای
  
  لكل(i=0; i<ډله_شمېر; ماتریکس_idx+=ډله_اندازه[i], i++){
    ډله_alpha = (فراغ *)&alpha_array[i * COMPSIZE];
    ډله_beta = (فراغ *)&beta_array[i * COMPSIZE];

    ډله_m = ډله_n = ډله_k = 0;
    ډله_lda = ډله_ldb = ډله_ldc = 0;
    ډله_transa = -1;
    ډله_transb = -1;
    info   = 0;
   
#إذا defined(CBLAS) 
    إذا (ترتیب == Cblas_ستون_Major) {
      ډله_m = m_array[i];
      ډله_n = n_array[i];
      ډله_k = k_array[i];

      ډله_lda = lda_array[i];
      ډله_ldb = ldb_array[i];
      ډله_ldc = ldc_array[i];
      
      إذا (transa_array[i] == CblasNoTrans)     ډله_transa = 0;
      إذا (transa_array[i] == CblasTrans)       ډله_transa = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) ډله_transa = 0;
      إذا (transa_array[i] == CblasConjTrans)   ډله_transa = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) ډله_transa = 2;
      إذا (transa_array[i] == CblasConjTrans)   ډله_transa = 3;
#د_که_پای
      إذا (transb_array[i] == CblasNoTrans)     ډله_transb = 0;
      إذا (transb_array[i] == CblasTrans)       ډله_transb = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) ډله_transb = 0;
      إذا (transb_array[i] == CblasConjTrans)   ډله_transb = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) ډله_transb = 2;
      إذا (transb_array[i] == CblasConjTrans)   ډله_transb = 3;
#د_که_پای
      ډله_nrowa = ډله_m;
      إذا (ډله_transa & 1) ډله_nrowa = ډله_k;
      ډله_nrowb = ډله_k;
      إذا (ډله_transb & 1) ډله_nrowb = ډله_n;

      info=-1;
	
      إذا (ډله_ldc < ډله_m) info = 13;
      إذا (ډله_ldb < ډله_nrowb)  info = 10;
      إذا (ډله_lda < ډله_nrowa)  info =  8;
      إذا (ډله_k < 0)        info =  5;
      إذا (ډله_n < 0)        info =  4;
      إذا (ډله_m < 0)        info =  3;
      إذا (ډله_transb < 0)        info =  2;
      إذا (ډله_transa < 0)        info =  1;

    }وإلا إذا (ترتیب == Cblas_لیکه_Major) {
      
      ډله_m = n_array[i];
      ډله_n = m_array[i];
      ډله_k = k_array[i];

      ډله_lda = ldb_array[i];
      ډله_ldb = lda_array[i];
      ډله_ldc = ldc_array[i];

      إذا (transb_array[i] == CblasNoTrans)     ډله_transa = 0;
      إذا (transb_array[i] == CblasTrans)       ډله_transa = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transb_array[i] == CblasConjNoTrans) ډله_transa = 0;
      إذا (transb_array[i] == CblasConjTrans)   ډله_transa = 1;
#وإلا
      إذا (transb_array[i] == CblasConjNoTrans) ډله_transa = 2;
      إذا (transb_array[i] == CblasConjTrans)   ډله_transa = 3;
#د_که_پای
      إذا (transa_array[i] == CblasNoTrans)     ډله_transb = 0;
      إذا (transa_array[i] == CblasTrans)       ډله_transb = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transa_array[i] == CblasConjNoTrans) ډله_transb = 0;
      إذا (transa_array[i] == CblasConjTrans)   ډله_transb = 1;
#وإلا
      إذا (transa_array[i] == CblasConjNoTrans) ډله_transb = 2;
      إذا (transa_array[i] == CblasConjTrans)   ډله_transb = 3;
#د_که_پای

#وإلا
      ډله_m = m_array[i];
      ډله_n = n_array[i];
      ډله_k = k_array[i];

      ډله_lda = lda_array[i];
      ډله_ldb = ldb_array[i];
      ډله_ldc = ldc_array[i];
      
      إذا (transa_array[i] == 'N')   ډله_transa = 0;
      إذا (transa_array[i] == 'T')   ډله_transa = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transa_array[i] == 'C')   ډله_transa = 1;
#وإلا
      إذا (transa_array[i] == 'C')   ډله_transa = 3;
#د_که_پای
      إذا (transb_array[i] == 'N')   ډله_transb = 0;
      إذا (transb_array[i] == 'T')   ډله_transb = 1;
#که_نه_تعریف_شوی COMPLEX
      إذا (transb_array[i] == 'C')   ډله_transb = 1;
#وإلا
      إذا (transb_array[i] == 'C')   ډله_transb = 3;
#د_که_پای
#د_که_پای

      ډله_nrowa = ډله_m;
      إذا (ډله_transa & 1) ډله_nrowa = ډله_k;
      ډله_nrowb = ډله_k;
      إذا (ډله_transb & 1) ډله_nrowb = ډله_n;

      info=-1;
	
      إذا (ډله_ldc < ډله_m) info = 13;
      إذا (ډله_ldb < ډله_nrowb)  info = 10;
      إذا (ډله_lda < ډله_nrowa)  info =  8;
      إذا (ډله_k < 0)        info =  5;
      إذا (ډله_n < 0)        info =  4;
      إذا (ډله_m < 0)        info =  3;
      إذا (ډله_transb < 0)        info =  2;
      إذا (ډله_transa < 0)        info =  1;      
#إذا defined(CBLAS)
    }
#د_که_پای

    إذا (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      حرر(دلیلونه_array);
      أرجع;
    }

    إذا (ډله_m == 0 || ډله_n == 0) استمر;

    ډله_ډېرپېښ=ډېرپېښ;

#إذا defined(SMP) || defined(SMALL_ماتریکس_OPT)
    MNK = (مزدوج) ډله_m * (مزدوج) ډله_n * (مزدوج) ډله_k;
#د_که_پای

#که_تعریف_شوی SMALL_ماتریکس_OPT
#إذا !defined(COMPLEX)
    إذا (*(FLOAT *)(ډله_alpha) != ZERO &&
        GEMM_SMALL_ماتریکس_PERMIT(ډله_transa, ډله_transb, ډله_m, ډله_n, ډله_k,
                                 *(FLOAT *)(ډله_alpha), *(FLOAT *)(ډله_beta))){
#وإلا
    إذا (((FLOAT *)(ډله_alpha))[0] != ZERO &&
        GEMM_SMALL_ماتریکس_PERMIT(ډله_transa, ډله_transb, ډله_m, ډله_n, ډله_k,
                                 ((FLOAT *)(ډله_alpha))[0], ((FLOAT *)(ډله_alpha))[1],
                                 ((FLOAT *)(ډله_beta))[0], ((FLOAT *)(ډله_beta))[1])){
#د_که_پای
      ډله_routine=عدم;
#إذا !defined(COMPLEX)
      إذا(*(FLOAT *)(ډله_beta) == 0.0){
	ډله_ډېرپېښ=ډېرپېښ | BLAS_SMALL_B0_OPT;
	ډله_small_ماتریکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((ډله_transb<<2)|ډله_transa));
      }وإلا{
	ډله_ډېرپېښ=ډېرپېښ | BLAS_SMALL_OPT;
	ډله_small_ماتریکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((ډله_transb<<2)|ډله_transa));
      }
#وإلا
      إذا(((FLOAT *)(ډله_beta))[0] == 0.0 && ((FLOAT *)(ډله_beta))[1] == 0.0){
	ډله_ډېرپېښ=ډېرپېښ | BLAS_SMALL_B0_OPT;
	ډله_small_ماتریکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((ډله_transb<<2)|ډله_transa));
      }وإلا{
	ډله_ډېرپېښ=ډېرپېښ | BLAS_SMALL_OPT;
	ډله_small_ماتریکس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((ډله_transb<<2)|ډله_transa));
      }

#د_که_پای
      
    }وإلا{
#د_که_پای
      ډله_routine=(فراغ*)(gemm[(ډله_transb<<2)|ډله_transa]);
#که_تعریف_شوی SMALL_ماتریکس_OPT
    }
#د_که_پای    

    
    لكل(j=0; j<ډله_اندازه[i]; j++){
      دلیلونه_array[شمېر].m=ډله_m;
      دلیلونه_array[شمېر].n=ډله_n;
      دلیلونه_array[شمېر].k=ډله_k;
      دلیلونه_array[شمېر].lda=ډله_lda;
      دلیلونه_array[شمېر].ldb=ډله_ldb;
      دلیلونه_array[شمېر].ldc=ډله_ldc;
      دلیلونه_array[شمېر].alpha=ډله_alpha;
      دلیلونه_array[شمېر].beta=ډله_beta;

#إذا defined(CBLAS)
      إذا (ترتیب == Cblas_لیکه_Major) {
	دلیلونه_array[شمېر].a=(b_array[ماتریکس_idx+j]);
	دلیلونه_array[شمېر].b=(a_array[ماتریکس_idx+j]);
      } وإلا
#د_که_پای
      {
	دلیلونه_array[شمېر].a=(a_array[ماتریکس_idx+j]);
	دلیلونه_array[شمېر].b=(b_array[ماتریکس_idx+j]);
      }
      
      دلیلونه_array[شمېر].c=(c_array[ماتریکس_idx+j]);
      
      دلیلونه_array[شمېر].routine_ډېرپېښ=ډله_ډېرپېښ;
      دلیلونه_array[شمېر].routine=ډله_routine;
#که_تعریف_شوی SMALL_ماتریکس_OPT
      إذا (!ډله_routine)
        دلیلونه_array[شمېر].routine=ډله_small_ماتریکس_opt_routine;
#د_که_پای      
      شمېر++;
    }
  }

  إذا(شمېر>0){
    GEMM_BATCH_THREAD(دلیلونه_array,شمېر);
  }

  حرر(دلیلونه_array);
}
