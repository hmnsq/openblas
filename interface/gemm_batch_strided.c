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

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل "common.h"

خالی openblas_warning(صحیح verbose, ثابت حرف * پیغام);

#اگر_غیر_تعریف_شدہ COMPLEX
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "QGEMM_BATCH_STRIDED "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "DGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD dgemm_batch_thread
#ورنہ
#تعریف_کرو ERROR_NAME "SGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD sgemm_batch_thread
#اگر_ختم
#ورنہ
#اگر_تعریف_شدہ XDOUBLE
#تعریف_کرو ERROR_NAME "XGEMM_BATCH_STRIDED "
#وگرنہ_اگر defined(DOUBLE)
#تعریف_کرو ERROR_NAME "ZGEMM_BATCH_STRIDED "
#تعریف_کرو GEMM_BATCH_THREAD zgemm_batch_thread
#ورنہ
#تعریف_کرو ERROR_NAME "CGEMM_BATCH_STRIDED "
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
خالی NAME(حرف *transa, حرف *transb,
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
#اگر !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#ورنہ
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#اگر_ختم
#ورنہ

خالی CNAME(شمار CBLAS_ترتیب ترتیب, شمار CBLAS_TRANSPOSE transa, شمار CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#اگر_غیر_تعریف_شدہ COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint گنتی) {
#ورنہ
	   خالی * valpha,
	   خالی * va, blasint lda, blasint stridea,
	   خالی * vb, blasint ldb, blasint strideb,
	   خالی * vbeta,
	   خالی * vc, blasint ldc, blasint stridec, blasint گنتی) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#اگر_ختم
#اگر_ختم
  BLASLONG گروہ_m, گروہ_n, گروہ_k;
  BLASLONG گروہ_lda, گروہ_ldb, گروہ_ldc;

  blas_arg_t * دلائل_array=تہی;

  صحیح منوال=0, گروہ_منوال=0;

  blasint i=0;

  صحیح گروہ_transa, گروہ_transb;
  BLASLONG گروہ_nrowa, گروہ_nrowb;
  blasint info;

  خالی * گروہ_routine=تہی;
#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
  خالی * گروہ_small_میٹرکس_opt_routine=تہی;
#اگر_ختم
  
#اگر defined (SMP) || defined(SMALL_میٹرکس_OPT)
  دوہرا MNK;
#اگر_ختم

  PRINT_DEBUG_CNAME;

  دلائل_array=(blas_arg_t *)حافظہ_دو(گنتی * حجم(blas_arg_t));
  
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
  
  برائے(i=0; i<گنتی; i++) {

    گروہ_m = گروہ_n = گروہ_k = 0;
    گروہ_lda = گروہ_ldb = گروہ_ldc = 0;
    گروہ_transa = -1;
    گروہ_transb = -1;
    info   = 0;
   
#اگر defined(CBLAS)
    اگر (ترتیب == Cblas_ستون_Major) {
      گروہ_m = m;
      گروہ_n = n;
      گروہ_k = k;

      گروہ_lda = lda;
      گروہ_ldb = ldb;
      گروہ_ldc = ldc;

      اگر (transa == CblasNoTrans)     گروہ_transa = 0;
      اگر (transa == CblasTrans)       گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa == CblasConjNoTrans) گروہ_transa = 0;
      اگر (transa == CblasConjTrans)   گروہ_transa = 1;
#ورنہ
      اگر (transa == CblasConjNoTrans) گروہ_transa = 2;
      اگر (transa == CblasConjTrans)   گروہ_transa = 3;
#اگر_ختم
      اگر (transb == CblasNoTrans)     گروہ_transb = 0;
      اگر (transb == CblasTrans)       گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transb == CblasConjNoTrans) گروہ_transb = 0;
      اگر (transb == CblasConjTrans)   گروہ_transb = 1;
#ورنہ
      اگر (transb == CblasConjNoTrans) گروہ_transb = 2;
      اگر (transb == CblasConjTrans)   گروہ_transb = 3;
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
      
      گروہ_m = n;
      گروہ_n = m;
      گروہ_k = k;

      گروہ_lda = ldb;
      گروہ_ldb = lda;
      گروہ_ldc = ldc;

      اگر (transb == CblasNoTrans)     گروہ_transa = 0;
      اگر (transb == CblasTrans)       گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transb == CblasConjNoTrans) گروہ_transa = 0;
      اگر (transb == CblasConjTrans)   گروہ_transa = 1;
#ورنہ
      اگر (transb == CblasConjNoTrans) گروہ_transa = 2;
      اگر (transb == CblasConjTrans)   گروہ_transa = 3;
#اگر_ختم
      اگر (transa == CblasNoTrans)     گروہ_transb = 0;
      اگر (transa == CblasTrans)       گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (transa == CblasConjNoTrans) گروہ_transb = 0;
      اگر (transa == CblasConjTrans)   گروہ_transb = 1;
#ورنہ
      اگر (transa == CblasConjNoTrans) گروہ_transb = 2;
      اگر (transa == CblasConjTrans)   گروہ_transb = 3;
#اگر_ختم

#ورنہ
      گروہ_m = m;
      گروہ_n = n;
      گروہ_k = k;

      گروہ_lda = lda;
      گروہ_ldb = ldb;
      گروہ_ldc = ldc;
      
      اگر (ta == 'N')   گروہ_transa = 0;
      اگر (ta == 'T')   گروہ_transa = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (ta == 'C')   گروہ_transa = 1;
#ورنہ
      اگر (ta == 'C')   گروہ_transa = 3;
#اگر_ختم
      اگر (tb == 'N')   گروہ_transb = 0;
      اگر (tb == 'T')   گروہ_transb = 1;
#اگر_غیر_تعریف_شدہ COMPLEX
      اگر (tb == 'C')   گروہ_transb = 1;
#ورنہ
      اگر (tb == 'C')   گروہ_transb = 3;
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

    /* Every member shares M and N: an empty member means an empty batch. */
    اگر (گروہ_m == 0 || گروہ_n == 0) { آزاد(دلائل_array); واپس; }

    گروہ_منوال=منوال;

#اگر defined(SMP) || defined(SMALL_میٹرکس_OPT)
    MNK = (دوہرا) گروہ_m * (دوہرا) گروہ_n * (دوہرا) گروہ_k;
#اگر_ختم

#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
#اگر !defined(COMPLEX)
    اگر (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#ورنہ
    اگر (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_میٹرکس_PERMIT(گروہ_transa, گروہ_transb, گروہ_m, گروہ_n, گروہ_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#اگر_ختم
      گروہ_routine=تہی;
#اگر !defined(COMPLEX)
      اگر(beta == 0.0){
	گروہ_منوال=منوال | BLAS_SMALL_B0_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروہ_transb<<2)|گروہ_transa));
      }ورنہ{
	گروہ_منوال=منوال | BLAS_SMALL_OPT;
	گروہ_small_میٹرکس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروہ_transb<<2)|گروہ_transa));
      }
#ورنہ
      اگر(beta[0] == 0.0 && beta[1] == 0.0){
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

    
      دلائل_array[i].m=گروہ_m;
      دلائل_array[i].n=گروہ_n;
      دلائل_array[i].k=گروہ_k;
      دلائل_array[i].lda=گروہ_lda;
      دلائل_array[i].ldb=گروہ_ldb;
      دلائل_array[i].ldc=گروہ_ldc;
#اگر_تعریف_شدہ COMPLEX
      دلائل_array[i].alpha=alpha;
      دلائل_array[i].beta=beta;
#ورنہ
      دلائل_array[i].alpha=&alpha;
      دلائل_array[i].beta=&beta;
#اگر_ختم

#اگر defined(CBLAS)
      اگر (ترتیب == Cblas_ستون_Major) {      
	دلائل_array[i].a=&(a[i*stridea]);
	دلائل_array[i].b=&(b[i*strideb]);
      }ورنہ اگر(ترتیب == Cblas_صف_Major){
	دلائل_array[i].a=&(b[i*strideb]);
	دلائل_array[i].b=&(a[i*stridea]);
      }
#ورنہ
	دلائل_array[i].a=&(a[i*stridea]);
	دلائل_array[i].b=&(b[i*strideb]);
#اگر_ختم
      
      دلائل_array[i].c= &c[i*stridec];
      
      دلائل_array[i].routine_منوال=گروہ_منوال;
      دلائل_array[i].routine=گروہ_routine;
#اگر_تعریف_شدہ SMALL_میٹرکس_OPT
      اگر (!گروہ_routine)
        دلائل_array[i].routine=گروہ_small_میٹرکس_opt_routine;
#اگر_ختم      
  }

  اگر(گنتی>0) {
    GEMM_BATCH_THREAD(دلائل_array,گنتی);
  }

  آزاد(دلائل_array);
}
