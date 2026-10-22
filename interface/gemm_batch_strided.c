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

فراغ openblas_warning(صحيح verbose, ثابت حرف * نياپو);

#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "QGEMM_BATCH_STRIDED "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "DGEMM_BATCH_STRIDED "
#وصف_ڪريو GEMM_BATCH_THREAD dgemm_batch_thread
#وإلا
#وصف_ڪريو ERROR_NAME "SGEMM_BATCH_STRIDED "
#وصف_ڪريو GEMM_BATCH_THREAD sgemm_batch_thread
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
#وصف_ڪريو ERROR_NAME "XGEMM_BATCH_STRIDED "
#وإلاإذا defined(DOUBLE)
#وصف_ڪريو ERROR_NAME "ZGEMM_BATCH_STRIDED "
#وصف_ڪريو GEMM_BATCH_THREAD zgemm_batch_thread
#وإلا
#وصف_ڪريو ERROR_NAME "CGEMM_BATCH_STRIDED "
#وصف_ڪريو GEMM_BATCH_THREAD cgemm_batch_thread
#جيڪڏهن_پورو
#جيڪڏهن_پورو
ساكن صحيح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#إذا defined(SMALL_ميٽرڪس_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#وصف_ڪريو USE_SMALL_ميٽرڪس_OPT 1
#وإلا
#وصف_ڪريو USE_SMALL_ميٽرڪس_OPT 0
#جيڪڏهن_پورو

#إذا USE_SMALL_ميٽرڪس_OPT
#جيڪڏهن_وصف_نه_ٿيل DYNAMIC_ARCH
#وصف_ڪريو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(جدول[idx]))
#وإلا
#وصف_ڪريو SMALL_KERNEL_ADDR(جدول, idx) ((فراغ *)(*(uintptr_t *)((حرف *)gotoblas + (size_t)(جدول[idx]))))
#جيڪڏهن_پورو


#جيڪڏهن_وصف_نه_ٿيل COMPLEX
ساكن size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ساكن size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#وصف_ڪريو GEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#وصف_ڪريو GEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
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

#وصف_ڪريو ZGEMM_SMALL_KERNEL(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#وصف_ڪريو ZGEMM_SMALL_KERNEL_B0(idx) (صحيح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#جيڪڏهن_پورو
#جيڪڏهن_پورو

OPENBLAS_EXPORT

#جيڪڏهن_وصف_نه_ٿيل CBLAS
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
	   blasint ڳڻپ = *matcount;
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
#جيڪڏهن_پورو
#وإلا

فراغ CNAME(تعداد CBLAS_ترتيب ترتيب, تعداد CBLAS_TRANSPOSE transa, تعداد CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint ڳڻپ) {
#وإلا
	   فراغ * valpha,
	   فراغ * va, blasint lda, blasint stridea,
	   فراغ * vb, blasint ldb, blasint strideb,
	   فراغ * vbeta,
	   فراغ * vc, blasint ldc, blasint stridec, blasint ڳڻپ) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
  BLASLONG گروهه_m, گروهه_n, گروهه_k;
  BLASLONG گروهه_lda, گروهه_ldb, گروهه_ldc;

  blas_arg_t * دليلون_array=عدم;

  صحيح گهڻائي=0, گروهه_گهڻائي=0;

  blasint i=0;

  صحيح گروهه_transa, گروهه_transb;
  BLASLONG گروهه_nrowa, گروهه_nrowb;
  blasint info;

  فراغ * گروهه_routine=عدم;
#جيڪڏهن_وصف_ٿيل SMALL_ميٽرڪس_OPT
  فراغ * گروهه_small_ميٽرڪس_opt_routine=عدم;
#جيڪڏهن_پورو
  
#إذا defined (SMP) || defined(SMALL_ميٽرڪس_OPT)
  مزدوج MNK;
#جيڪڏهن_پورو

  PRINT_DEBUG_CNAME;

  دليلون_array=(blas_arg_t *)يادگيري_ڏيو(ڳڻپ * حجم(blas_arg_t));
  
  إذا(دليلون_array == عدم){
    openblas_warning(0, "memory alloc failed!\n");
    أرجع;
  }

#جيڪڏهن_وصف_ٿيل SMP
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
#جيڪڏهن_وصف_ٿيل XDOUBLE
  گهڻائي  =  BLAS_XDOUBLE | BLAS_REAL;
#وإلاإذا defined(DOUBLE)
  گهڻائي  =  BLAS_DOUBLE  | BLAS_REAL;
#وإلا
  گهڻائي  =  BLAS_SINGLE  | BLAS_REAL;
#جيڪڏهن_پورو
#وإلا
#جيڪڏهن_وصف_ٿيل XDOUBLE
  گهڻائي  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وإلاإذا defined(DOUBLE)
  گهڻائي  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وإلا
  گهڻائي  =  BLAS_SINGLE  | BLAS_COMPLEX;
#جيڪڏهن_پورو
#جيڪڏهن_پورو
#جيڪڏهن_پورو
  
  لكل(i=0; i<ڳڻپ; i++) {

    گروهه_m = گروهه_n = گروهه_k = 0;
    گروهه_lda = گروهه_ldb = گروهه_ldc = 0;
    گروهه_transa = -1;
    گروهه_transb = -1;
    info   = 0;
   
#إذا defined(CBLAS)
    إذا (ترتيب == Cblas_کانو_Major) {
      گروهه_m = m;
      گروهه_n = n;
      گروهه_k = k;

      گروهه_lda = lda;
      گروهه_ldb = ldb;
      گروهه_ldc = ldc;

      إذا (transa == CblasNoTrans)     گروهه_transa = 0;
      إذا (transa == CblasTrans)       گروهه_transa = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (transa == CblasConjNoTrans) گروهه_transa = 0;
      إذا (transa == CblasConjTrans)   گروهه_transa = 1;
#وإلا
      إذا (transa == CblasConjNoTrans) گروهه_transa = 2;
      إذا (transa == CblasConjTrans)   گروهه_transa = 3;
#جيڪڏهن_پورو
      إذا (transb == CblasNoTrans)     گروهه_transb = 0;
      إذا (transb == CblasTrans)       گروهه_transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (transb == CblasConjNoTrans) گروهه_transb = 0;
      إذا (transb == CblasConjTrans)   گروهه_transb = 1;
#وإلا
      إذا (transb == CblasConjNoTrans) گروهه_transb = 2;
      إذا (transb == CblasConjTrans)   گروهه_transb = 3;
#جيڪڏهن_پورو
      گروهه_nrowa = گروهه_m;
      إذا (گروهه_transa & 1) گروهه_nrowa = گروهه_k;
      گروهه_nrowb = گروهه_k;
      إذا (گروهه_transb & 1) گروهه_nrowb = گروهه_n;

      info=-1;
	
      إذا (گروهه_ldc < گروهه_m) info = 13;
      إذا (گروهه_ldb < گروهه_nrowb)  info = 10;
      إذا (گروهه_lda < گروهه_nrowa)  info =  8;
      إذا (گروهه_k < 0)        info =  5;
      إذا (گروهه_n < 0)        info =  4;
      إذا (گروهه_m < 0)        info =  3;
      إذا (گروهه_transb < 0)        info =  2;
      إذا (گروهه_transa < 0)        info =  1;

    }وإلا إذا (ترتيب == Cblas_صف_Major) {
      
      گروهه_m = n;
      گروهه_n = m;
      گروهه_k = k;

      گروهه_lda = ldb;
      گروهه_ldb = lda;
      گروهه_ldc = ldc;

      إذا (transb == CblasNoTrans)     گروهه_transa = 0;
      إذا (transb == CblasTrans)       گروهه_transa = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (transb == CblasConjNoTrans) گروهه_transa = 0;
      إذا (transb == CblasConjTrans)   گروهه_transa = 1;
#وإلا
      إذا (transb == CblasConjNoTrans) گروهه_transa = 2;
      إذا (transb == CblasConjTrans)   گروهه_transa = 3;
#جيڪڏهن_پورو
      إذا (transa == CblasNoTrans)     گروهه_transb = 0;
      إذا (transa == CblasTrans)       گروهه_transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (transa == CblasConjNoTrans) گروهه_transb = 0;
      إذا (transa == CblasConjTrans)   گروهه_transb = 1;
#وإلا
      إذا (transa == CblasConjNoTrans) گروهه_transb = 2;
      إذا (transa == CblasConjTrans)   گروهه_transb = 3;
#جيڪڏهن_پورو

#وإلا
      گروهه_m = m;
      گروهه_n = n;
      گروهه_k = k;

      گروهه_lda = lda;
      گروهه_ldb = ldb;
      گروهه_ldc = ldc;
      
      إذا (ta == 'N')   گروهه_transa = 0;
      إذا (ta == 'T')   گروهه_transa = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (ta == 'C')   گروهه_transa = 1;
#وإلا
      إذا (ta == 'C')   گروهه_transa = 3;
#جيڪڏهن_پورو
      إذا (tb == 'N')   گروهه_transb = 0;
      إذا (tb == 'T')   گروهه_transb = 1;
#جيڪڏهن_وصف_نه_ٿيل COMPLEX
      إذا (tb == 'C')   گروهه_transb = 1;
#وإلا
      إذا (tb == 'C')   گروهه_transb = 3;
#جيڪڏهن_پورو
#جيڪڏهن_پورو

      گروهه_nrowa = گروهه_m;
      إذا (گروهه_transa & 1) گروهه_nrowa = گروهه_k;
      گروهه_nrowb = گروهه_k;
      إذا (گروهه_transb & 1) گروهه_nrowb = گروهه_n;

      info=-1;
	
      إذا (گروهه_ldc < گروهه_m) info = 13;
      إذا (گروهه_ldb < گروهه_nrowb)  info = 10;
      إذا (گروهه_lda < گروهه_nrowa)  info =  8;
      إذا (گروهه_k < 0)        info =  5;
      إذا (گروهه_n < 0)        info =  4;
      إذا (گروهه_m < 0)        info =  3;
      إذا (گروهه_transb < 0)        info =  2;
      إذا (گروهه_transa < 0)        info =  1;      
#إذا defined(CBLAS)
    }
#جيڪڏهن_پورو

    إذا (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, حجم(ERROR_NAME));
      حرر(دليلون_array);
      أرجع;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    إذا (گروهه_m == 0 || گروهه_n == 0) { حرر(دليلون_array); أرجع; }

    گروهه_گهڻائي=گهڻائي;

#إذا defined(SMP) || defined(SMALL_ميٽرڪس_OPT)
    MNK = (مزدوج) گروهه_m * (مزدوج) گروهه_n * (مزدوج) گروهه_k;
#جيڪڏهن_پورو

#جيڪڏهن_وصف_ٿيل SMALL_ميٽرڪس_OPT
#إذا !defined(COMPLEX)
    إذا (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_ميٽرڪس_PERMIT(گروهه_transa, گروهه_transb, گروهه_m, گروهه_n, گروهه_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#وإلا
    إذا (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_ميٽرڪس_PERMIT(گروهه_transa, گروهه_transb, گروهه_m, گروهه_n, گروهه_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#جيڪڏهن_پورو
      گروهه_routine=عدم;
#إذا !defined(COMPLEX)
      إذا(beta == 0.0){
	گروهه_گهڻائي=گهڻائي | BLAS_SMALL_B0_OPT;
	گروهه_small_ميٽرڪس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((گروهه_transb<<2)|گروهه_transa));
      }وإلا{
	گروهه_گهڻائي=گهڻائي | BLAS_SMALL_OPT;
	گروهه_small_ميٽرڪس_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((گروهه_transb<<2)|گروهه_transa));
      }
#وإلا
      إذا(beta[0] == 0.0 && beta[1] == 0.0){
	گروهه_گهڻائي=گهڻائي | BLAS_SMALL_B0_OPT;
	گروهه_small_ميٽرڪس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((گروهه_transb<<2)|گروهه_transa));
      }وإلا{
	گروهه_گهڻائي=گهڻائي | BLAS_SMALL_OPT;
	گروهه_small_ميٽرڪس_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((گروهه_transb<<2)|گروهه_transa));
      }

#جيڪڏهن_پورو
      
    }وإلا{
#جيڪڏهن_پورو
      گروهه_routine=(فراغ*)(gemm[(گروهه_transb<<2)|گروهه_transa]);
#جيڪڏهن_وصف_ٿيل SMALL_ميٽرڪس_OPT
    }
#جيڪڏهن_پورو    

    
      دليلون_array[i].m=گروهه_m;
      دليلون_array[i].n=گروهه_n;
      دليلون_array[i].k=گروهه_k;
      دليلون_array[i].lda=گروهه_lda;
      دليلون_array[i].ldb=گروهه_ldb;
      دليلون_array[i].ldc=گروهه_ldc;
#جيڪڏهن_وصف_ٿيل COMPLEX
      دليلون_array[i].alpha=alpha;
      دليلون_array[i].beta=beta;
#وإلا
      دليلون_array[i].alpha=&alpha;
      دليلون_array[i].beta=&beta;
#جيڪڏهن_پورو

#إذا defined(CBLAS)
      إذا (ترتيب == Cblas_کانو_Major) {      
	دليلون_array[i].a=&(a[i*stridea]);
	دليلون_array[i].b=&(b[i*strideb]);
      }وإلا إذا(ترتيب == Cblas_صف_Major){
	دليلون_array[i].a=&(b[i*strideb]);
	دليلون_array[i].b=&(a[i*stridea]);
      }
#وإلا
	دليلون_array[i].a=&(a[i*stridea]);
	دليلون_array[i].b=&(b[i*strideb]);
#جيڪڏهن_پورو
      
      دليلون_array[i].c= &c[i*stridec];
      
      دليلون_array[i].routine_گهڻائي=گروهه_گهڻائي;
      دليلون_array[i].routine=گروهه_routine;
#جيڪڏهن_وصف_ٿيل SMALL_ميٽرڪس_OPT
      إذا (!گروهه_routine)
        دليلون_array[i].routine=گروهه_small_ميٽرڪس_opt_routine;
#جيڪڏهن_پورو      
  }

  إذا(ڳڻپ>0) {
    GEMM_BATCH_THREAD(دليلون_array,ڳڻپ);
  }

  حرر(دليلون_array);
}
