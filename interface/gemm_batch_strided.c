/*****************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
#تعریف_کرو ERROR_NAME "QGEMM_BATCH_STRIDED "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD dgemm_batch_thread
#وإلا
#تعریف_کرو ERROR_NAME "SGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD sgemm_batch_thread
#جے_انت
#وإلا
#جے_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM_BATCH_STRIDED "
#وإلاإذا defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD zgemm_batch_thread
#وإلا
#تعریف_کرو ERROR_NAME "CGEMM_BATCH_STRIDED "
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
فراغ NAME(حرف *transa, حرف *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           حرف ta = *transa;
   	   حرف tb = *transb;	   
	   blasint گنتی = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#إذا !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#وإلا
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#جے_انت
#وإلا

فراغ CNAME(تعداد CBLAS_ترتیب ترتیب, تعداد CBLAS_TRANSPOSE transa, تعداد CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#جے_بنا_تعریف COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint گنتی) {
#وإلا
	   فراغ * valpha,
	   فراغ * va, blasint lda, blasint stridea,
	   فراغ * vb, blasint ldb, blasint strideb,
	   فراغ * vbeta,
	   فراغ * vc, blasint ldc, blasint stridec, blasint گنتی) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#جے_انت
#جے_انت
  BLASLONG گروہ_m, گروہ_n, گروہ_k;
  BLASLONG گروہ_lda, گروہ_ldb, گروہ_ldc;

  blas_arg_t * دلیلاں_array=عدم;

  صحيح منوال=0, گروہ_منوال=0;

  blasint i=0;

  صحيح گروہ_transa, گروہ_transb;
  BLASLONG گروہ_nrowa, گروہ_nrowb;
  blasint info;

  فراغ * گروہ_routine=عدم;
#جے_تعریف_شدہ SMALL_میٹرکس_OPT
  فراغ * گروہ_small_میٹرکس_opt_routine=عدم;
#جے_انت
  
#إذا defined (SMP) || defined(SMALL_میٹرکس_OPT)
  مزدوج MNK;
#جے_انت

  PRINT_DEBUG_CNAME;

  دلیلاں_array=(blas_arg_t *)حافظہ_دیو(گنتی * حجم(blas_arg_t));
  
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
  
  لكل(i=0; i<گنتی; i++) {

    گروہ_m = گروہ_n = گروہ_k = 0;
    گروہ_lda = گروہ_ldb = گروہ_ldc = 0;
    گروہ_transa = -1;
    گروہ_transb = -1;
    info   = 0;
   
#إذا defined(CBLAS)
    إذا (ترتیب == Cblas_خانہ_Major) {
      گروہ_m = m;
      گروہ_n = n;
      گروہ_k = k;

      گروہ_lda = lda;
      گروہ_ldb = ldb;
      گروہ_ldc = ldc;

      إذا (transa == CblasNoTrans)     گروہ_transa = 0;
      إذا (transa == CblasTrans)       گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transa == CblasConjNoTrans) گروہ_transa = 0;
      إذا (transa == CblasConjTrans)   گروہ_transa = 1;
#وإلا
      إذا (transa == CblasConjNoTrans) گروہ_transa = 2;
      إذا (transa == CblasConjTrans)   گروہ_transa = 3;
#جے_انت
      إذا (transb == CblasNoTrans)     گروہ_transb = 0;
      إذا (transb == CblasTrans)       گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transb == CblasConjNoTrans) گروہ_transb = 0;
      إذا (transb == CblasConjTrans)   گروہ_transb = 1;
#وإلا
      إذا (transb == CblasConjNoTrans) گروہ_transb = 2;
      إذا (transb == CblasConjTrans)   گروہ_transb = 3;
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
      
      گروہ_m = n;
      گروہ_n = m;
      گروہ_k = k;

      گروہ_lda = ldb;
      گروہ_ldb = lda;
      گروہ_ldc = ldc;

      إذا (transb == CblasNoTrans)     گروہ_transa = 0;
      إذا (transb == CblasTrans)       گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transb == CblasConjNoTrans) گروہ_transa = 0;
      إذا (transb == CblasConjTrans)   گروہ_transa = 1;
#وإلا
      إذا (transb == CblasConjNoTrans) گروہ_transa = 2;
      إذا (transb == CblasConjTrans)   گروہ_transa = 3;
#جے_انت
      إذا (transa == CblasNoTrans)     گروہ_transb = 0;
      إذا (transa == CblasTrans)       گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (transa == CblasConjNoTrans) گروہ_transb = 0;
      إذا (transa == CblasConjTrans)   گروہ_transb = 1;
#وإلا
      إذا (transa == CblasConjNoTrans) گروہ_transb = 2;
      إذا (transa == CblasConjTrans)   گروہ_transb = 3;
#جے_انت

#وإلا
      گروہ_m = m;
      گروہ_n = n;
      گروہ_k = k;

      گروہ_lda = lda;
      گروہ_ldb = ldb;
      گروہ_ldc = ldc;
      
      إذا (ta == 'N')   گروہ_transa = 0;
      إذا (ta == 'T')   گروہ_transa = 1;
#جے_بنا_تعریف COMPLEX
      إذا (ta == 'C')   گروہ_transa = 1;
#وإلا
      إذا (ta == 'C')   گروہ_transa = 3;
#جے_انت
      إذا (tb == 'N')   گروہ_transb = 0;
      إذا (tb == 'T')   گروہ_transb = 1;
#جے_بنا_تعریف COMPLEX
      إذا (tb == 'C')   گروہ_transb = 1;
#وإلا
      إذا (tb == 'C')   گروہ_transb = 3;
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

    /* Every member shares M and N: an empty member means an empty batch. */
    إذا (گروہ_m == 0 || گروہ_n == 0) { حرر(دلیلاں_array); أرجع; }

    گروہ_منوال=منوال;

#إذا defined(SMP) || defined(SMALL_میٹرکس_OPT)
    MNK = (مزدوج) گروہ_m * (مزدوج) گروہ_n * (مزدوج) گروہ_k;
#جے_انت

#جے_تعریف_شدہ SMALL_میٹرکس_OPT
#إذا !defined(COMPLEX)
    إذا (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#وإلا
    إذا (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#جے_انت
      گروہ_routine=عدم;
#إذا !defined(COMPLEX)
      إذا(beta == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }وإلا{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }
#وإلا
      إذا(beta[0] == 0.0 && beta[1] == 0.0){
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

    
      دلیلاں_array[i].m=گروہ_m;
      دلیلاں_array[i].n=گروہ_n;
      دلیلاں_array[i].k=گروہ_k;
      دلیلاں_array[i].lda=گروہ_lda;
      دلیلاں_array[i].ldb=گروہ_ldb;
      دلیلاں_array[i].ldc=گروہ_ldc;
#جے_تعریف_شدہ COMPLEX
      دلیلاں_array[i].alpha=alpha;
      دلیلاں_array[i].beta=beta;
#وإلا
      دلیلاں_array[i].alpha=&alpha;
      دلیلاں_array[i].beta=&beta;
#جے_انت

#إذا defined(CBLAS)
      إذا (ترتیب == Cblas_خانہ_Major) {      
	دلیلاں_array[i].a=&(a[i*stridea]);
	دلیلاں_array[i].b=&(b[i*strideb]);
      }وإلا إذا(ترتیب == Cblas_صف_Major){
	دلیلاں_array[i].a=&(b[i*strideb]);
	دلیلاں_array[i].b=&(a[i*stridea]);
      }
#وإلا
	دلیلاں_array[i].a=&(a[i*stridea]);
	دلیلاں_array[i].b=&(b[i*strideb]);
#جے_انت
      
      دلیلاں_array[i].c= &c[i*stridec];
      
      دلیلاں_array[i].routine_منوال=گروہ_منوال;
      دلیلاں_array[i].routine=گروہ_routine;
#جے_تعریف_شدہ SMALL_میٹرکس_OPT
      إذا (!گروہ_routine)
        دلیلاں_array[i].routine=گروہ_small_میٹرکس_opt_routine;
#جے_انت      
  }

  إذا(گنتی>0) {
    GEMM_BATCH_THREAD(دلیلاں_array,گنتی);
  }

  حرر(دلیلاں_array);
}
