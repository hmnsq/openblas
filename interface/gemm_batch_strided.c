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

#ਅਨ੍ਤਰ੍ਭਾਵ <stdio.h>
#ਅਨ੍ਤਰ੍ਭਾਵ <stdlib.h>
#ਅਨ੍ਤਰ੍ਭਾਵ "common.h"

ਖਾਲੀ openblas_warning(ਪੂਰਨ verbose, ਅਚਲ ਅੱਖਰ * ਸੁਨੇਹਾ);

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "QGEMM_BATCH_STRIDED "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "DGEMM_BATCH_STRIDED "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_BATCH_THREAD dgemm_batch_thread
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "SGEMM_BATCH_STRIDED "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_BATCH_THREAD sgemm_batch_thread
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "XGEMM_BATCH_STRIDED "
#ਜਾਂਜੇ defined(DOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "ZGEMM_BATCH_STRIDED "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_BATCH_THREAD zgemm_batch_thread
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ERROR_NAME "CGEMM_BATCH_STRIDED "
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_BATCH_THREAD cgemm_batch_thread
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
ਸਥਿਰ ਪੂਰਨ (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ਜੇ defined(SMALL_ਮੈਟ੍ਰਿਕਸ_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT 1
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT 0
#ਜੇ_ਅੰਤ

#ਜੇ USE_SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ DYNAMIC_ARCH
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMALL_KERNEL_ADDR(ਸਾਰਣੀ, idx) ((ਖਾਲੀ *)(ਸਾਰਣੀ[idx]))
#ਨਹੀਂਤਾਂ
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ SMALL_KERNEL_ADDR(ਸਾਰਣੀ, idx) ((ਖਾਲੀ *)(*(uintptr_t *)((ਅੱਖਰ *)gotoblas + (size_t)(ਸਾਰਣੀ[idx]))))
#ਜੇ_ਅੰਤ


#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
ਸਥਿਰ size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ਸਥਿਰ size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_SMALL_KERNEL_B0(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ GEMM_SMALL_KERNEL(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#ਨਹੀਂਤਾਂ

ਸਥਿਰ size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ਸਥਿਰ size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ZGEMM_SMALL_KERNEL(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#ਪਰਿਭਾਸ਼ਿਤ_ਕਰੋ ZGEMM_SMALL_KERNEL_B0(idx) (ਪੂਰਨ (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

OPENBLAS_EXPORT

#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ CBLAS
ਖਾਲੀ NAME(ਅੱਖਰ *transa, ਅੱਖਰ *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           ਅੱਖਰ ta = *transa;
   	   ਅੱਖਰ tb = *transb;	   
	   blasint ਗਿਣਤੀ_2 = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#ਜੇ !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#ਨਹੀਂਤਾਂ
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ

ਖਾਲੀ CNAME(ਗਿਣਤੀ CBLAS_ਕ੍ਰਮ ਕ੍ਰਮ, ਗਿਣਤੀ CBLAS_TRANSPOSE transa, ਗਿਣਤੀ CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint ਗਿਣਤੀ_2) {
#ਨਹੀਂਤਾਂ
	   ਖਾਲੀ * valpha,
	   ਖਾਲੀ * va, blasint lda, blasint stridea,
	   ਖਾਲੀ * vb, blasint ldb, blasint strideb,
	   ਖਾਲੀ * vbeta,
	   ਖਾਲੀ * vc, blasint ldc, blasint stridec, blasint ਗਿਣਤੀ_2) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
  BLASLONG ਸਮੂਹ_m, ਸਮੂਹ_n, ਸਮੂਹ_k;
  BLASLONG ਸਮੂਹ_lda, ਸਮੂਹ_ldb, ਸਮੂਹ_ldc;

  blas_arg_t * ਦਲੀਲਾਂ_array=ਨਿਰਾ;

  ਪੂਰਨ ਬਹੁਲਕ=0, ਸਮੂਹ_ਬਹੁਲਕ=0;

  blasint i=0;

  ਪੂਰਨ ਸਮੂਹ_transa, ਸਮੂਹ_transb;
  BLASLONG ਸਮੂਹ_nrowa, ਸਮੂਹ_nrowb;
  blasint info;

  ਖਾਲੀ * ਸਮੂਹ_routine=ਨਿਰਾ;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
  ਖਾਲੀ * ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine=ਨਿਰਾ;
#ਜੇ_ਅੰਤ
  
#ਜੇ defined (SMP) || defined(SMALL_ਮੈਟ੍ਰਿਕਸ_OPT)
  ਦੁੱਗਣਾ MNK;
#ਜੇ_ਅੰਤ

  PRINT_DEBUG_CNAME;

  ਦਲੀਲਾਂ_array=(blas_arg_t *)ਸਮ੍ਰਿਤੀ_ਦਿਓ(ਗਿਣਤੀ_2 * ਆਕਾਰ(blas_arg_t));
  
  ਜੇ(ਦਲੀਲਾਂ_array == ਨਿਰਾ){
    openblas_warning(0, "memory alloc failed!\n");
    ਵਾਪਸ;
  }

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMP
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_REAL;
#ਜਾਂਜੇ defined(DOUBLE)
  ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_REAL;
#ਨਹੀਂਤਾਂ
  ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_REAL;
#ਜੇ_ਅੰਤ
#ਨਹੀਂਤਾਂ
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ XDOUBLE
  ਬਹੁਲਕ  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#ਜਾਂਜੇ defined(DOUBLE)
  ਬਹੁਲਕ  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#ਨਹੀਂਤਾਂ
  ਬਹੁਲਕ  =  BLAS_SINGLE  | BLAS_COMPLEX;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ
  
  ਲਈ(i=0; i<ਗਿਣਤੀ_2; i++) {

    ਸਮੂਹ_m = ਸਮੂਹ_n = ਸਮੂਹ_k = 0;
    ਸਮੂਹ_lda = ਸਮੂਹ_ldb = ਸਮੂਹ_ldc = 0;
    ਸਮੂਹ_transa = -1;
    ਸਮੂਹ_transb = -1;
    info   = 0;
   
#ਜੇ defined(CBLAS)
    ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {
      ਸਮੂਹ_m = m;
      ਸਮੂਹ_n = n;
      ਸਮੂਹ_k = k;

      ਸਮੂਹ_lda = lda;
      ਸਮੂਹ_ldb = ldb;
      ਸਮੂਹ_ldc = ldc;

      ਜੇ (transa == CblasNoTrans)     ਸਮੂਹ_transa = 0;
      ਜੇ (transa == CblasTrans)       ਸਮੂਹ_transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (transa == CblasConjNoTrans) ਸਮੂਹ_transa = 0;
      ਜੇ (transa == CblasConjTrans)   ਸਮੂਹ_transa = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (transa == CblasConjNoTrans) ਸਮੂਹ_transa = 2;
      ਜੇ (transa == CblasConjTrans)   ਸਮੂਹ_transa = 3;
#ਜੇ_ਅੰਤ
      ਜੇ (transb == CblasNoTrans)     ਸਮੂਹ_transb = 0;
      ਜੇ (transb == CblasTrans)       ਸਮੂਹ_transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (transb == CblasConjNoTrans) ਸਮੂਹ_transb = 0;
      ਜੇ (transb == CblasConjTrans)   ਸਮੂਹ_transb = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (transb == CblasConjNoTrans) ਸਮੂਹ_transb = 2;
      ਜੇ (transb == CblasConjTrans)   ਸਮੂਹ_transb = 3;
#ਜੇ_ਅੰਤ
      ਸਮੂਹ_nrowa = ਸਮੂਹ_m;
      ਜੇ (ਸਮੂਹ_transa & 1) ਸਮੂਹ_nrowa = ਸਮੂਹ_k;
      ਸਮੂਹ_nrowb = ਸਮੂਹ_k;
      ਜੇ (ਸਮੂਹ_transb & 1) ਸਮੂਹ_nrowb = ਸਮੂਹ_n;

      info=-1;
	
      ਜੇ (ਸਮੂਹ_ldc < ਸਮੂਹ_m) info = 13;
      ਜੇ (ਸਮੂਹ_ldb < ਸਮੂਹ_nrowb)  info = 10;
      ਜੇ (ਸਮੂਹ_lda < ਸਮੂਹ_nrowa)  info =  8;
      ਜੇ (ਸਮੂਹ_k < 0)        info =  5;
      ਜੇ (ਸਮੂਹ_n < 0)        info =  4;
      ਜੇ (ਸਮੂਹ_m < 0)        info =  3;
      ਜੇ (ਸਮੂਹ_transb < 0)        info =  2;
      ਜੇ (ਸਮੂਹ_transa < 0)        info =  1;

    }ਨਹੀਂਤਾਂ ਜੇ (ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major) {
      
      ਸਮੂਹ_m = n;
      ਸਮੂਹ_n = m;
      ਸਮੂਹ_k = k;

      ਸਮੂਹ_lda = ldb;
      ਸਮੂਹ_ldb = lda;
      ਸਮੂਹ_ldc = ldc;

      ਜੇ (transb == CblasNoTrans)     ਸਮੂਹ_transa = 0;
      ਜੇ (transb == CblasTrans)       ਸਮੂਹ_transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (transb == CblasConjNoTrans) ਸਮੂਹ_transa = 0;
      ਜੇ (transb == CblasConjTrans)   ਸਮੂਹ_transa = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (transb == CblasConjNoTrans) ਸਮੂਹ_transa = 2;
      ਜੇ (transb == CblasConjTrans)   ਸਮੂਹ_transa = 3;
#ਜੇ_ਅੰਤ
      ਜੇ (transa == CblasNoTrans)     ਸਮੂਹ_transb = 0;
      ਜੇ (transa == CblasTrans)       ਸਮੂਹ_transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (transa == CblasConjNoTrans) ਸਮੂਹ_transb = 0;
      ਜੇ (transa == CblasConjTrans)   ਸਮੂਹ_transb = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (transa == CblasConjNoTrans) ਸਮੂਹ_transb = 2;
      ਜੇ (transa == CblasConjTrans)   ਸਮੂਹ_transb = 3;
#ਜੇ_ਅੰਤ

#ਨਹੀਂਤਾਂ
      ਸਮੂਹ_m = m;
      ਸਮੂਹ_n = n;
      ਸਮੂਹ_k = k;

      ਸਮੂਹ_lda = lda;
      ਸਮੂਹ_ldb = ldb;
      ਸਮੂਹ_ldc = ldc;
      
      ਜੇ (ta == 'N')   ਸਮੂਹ_transa = 0;
      ਜੇ (ta == 'T')   ਸਮੂਹ_transa = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (ta == 'C')   ਸਮੂਹ_transa = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (ta == 'C')   ਸਮੂਹ_transa = 3;
#ਜੇ_ਅੰਤ
      ਜੇ (tb == 'N')   ਸਮੂਹ_transb = 0;
      ਜੇ (tb == 'T')   ਸਮੂਹ_transb = 1;
#ਜੇ_ਅਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਜੇ (tb == 'C')   ਸਮੂਹ_transb = 1;
#ਨਹੀਂਤਾਂ
      ਜੇ (tb == 'C')   ਸਮੂਹ_transb = 3;
#ਜੇ_ਅੰਤ
#ਜੇ_ਅੰਤ

      ਸਮੂਹ_nrowa = ਸਮੂਹ_m;
      ਜੇ (ਸਮੂਹ_transa & 1) ਸਮੂਹ_nrowa = ਸਮੂਹ_k;
      ਸਮੂਹ_nrowb = ਸਮੂਹ_k;
      ਜੇ (ਸਮੂਹ_transb & 1) ਸਮੂਹ_nrowb = ਸਮੂਹ_n;

      info=-1;
	
      ਜੇ (ਸਮੂਹ_ldc < ਸਮੂਹ_m) info = 13;
      ਜੇ (ਸਮੂਹ_ldb < ਸਮੂਹ_nrowb)  info = 10;
      ਜੇ (ਸਮੂਹ_lda < ਸਮੂਹ_nrowa)  info =  8;
      ਜੇ (ਸਮੂਹ_k < 0)        info =  5;
      ਜੇ (ਸਮੂਹ_n < 0)        info =  4;
      ਜੇ (ਸਮੂਹ_m < 0)        info =  3;
      ਜੇ (ਸਮੂਹ_transb < 0)        info =  2;
      ਜੇ (ਸਮੂਹ_transa < 0)        info =  1;      
#ਜੇ defined(CBLAS)
    }
#ਜੇ_ਅੰਤ

    ਜੇ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, ਆਕਾਰ(ERROR_NAME));
      ਮੁਕਤ(ਦਲੀਲਾਂ_array);
      ਵਾਪਸ;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    ਜੇ (ਸਮੂਹ_m == 0 || ਸਮੂਹ_n == 0) { ਮੁਕਤ(ਦਲੀਲਾਂ_array); ਵਾਪਸ; }

    ਸਮੂਹ_ਬਹੁਲਕ=ਬਹੁਲਕ;

#ਜੇ defined(SMP) || defined(SMALL_ਮੈਟ੍ਰਿਕਸ_OPT)
    MNK = (ਦੁੱਗਣਾ) ਸਮੂਹ_m * (ਦੁੱਗਣਾ) ਸਮੂਹ_n * (ਦੁੱਗਣਾ) ਸਮੂਹ_k;
#ਜੇ_ਅੰਤ

#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
#ਜੇ !defined(COMPLEX)
    ਜੇ (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_ਮੈਟ੍ਰਿਕਸ_PERMIT(ਸਮੂਹ_transa, ਸਮੂਹ_transb, ਸਮੂਹ_m, ਸਮੂਹ_n, ਸਮੂਹ_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#ਨਹੀਂਤਾਂ
    ਜੇ (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_ਮੈਟ੍ਰਿਕਸ_PERMIT(ਸਮੂਹ_transa, ਸਮੂਹ_transb, ਸਮੂਹ_m, ਸਮੂਹ_n, ਸਮੂਹ_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#ਜੇ_ਅੰਤ
      ਸਮੂਹ_routine=ਨਿਰਾ;
#ਜੇ !defined(COMPLEX)
      ਜੇ(beta == 0.0){
	ਸਮੂਹ_ਬਹੁਲਕ=ਬਹੁਲਕ | BLAS_SMALL_B0_OPT;
	ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((ਸਮੂਹ_transb<<2)|ਸਮੂਹ_transa));
      }ਨਹੀਂਤਾਂ{
	ਸਮੂਹ_ਬਹੁਲਕ=ਬਹੁਲਕ | BLAS_SMALL_OPT;
	ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((ਸਮੂਹ_transb<<2)|ਸਮੂਹ_transa));
      }
#ਨਹੀਂਤਾਂ
      ਜੇ(beta[0] == 0.0 && beta[1] == 0.0){
	ਸਮੂਹ_ਬਹੁਲਕ=ਬਹੁਲਕ | BLAS_SMALL_B0_OPT;
	ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((ਸਮੂਹ_transb<<2)|ਸਮੂਹ_transa));
      }ਨਹੀਂਤਾਂ{
	ਸਮੂਹ_ਬਹੁਲਕ=ਬਹੁਲਕ | BLAS_SMALL_OPT;
	ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((ਸਮੂਹ_transb<<2)|ਸਮੂਹ_transa));
      }

#ਜੇ_ਅੰਤ
      
    }ਨਹੀਂਤਾਂ{
#ਜੇ_ਅੰਤ
      ਸਮੂਹ_routine=(ਖਾਲੀ*)(gemm[(ਸਮੂਹ_transb<<2)|ਸਮੂਹ_transa]);
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
    }
#ਜੇ_ਅੰਤ    

    
      ਦਲੀਲਾਂ_array[i].m=ਸਮੂਹ_m;
      ਦਲੀਲਾਂ_array[i].n=ਸਮੂਹ_n;
      ਦਲੀਲਾਂ_array[i].k=ਸਮੂਹ_k;
      ਦਲੀਲਾਂ_array[i].lda=ਸਮੂਹ_lda;
      ਦਲੀਲਾਂ_array[i].ldb=ਸਮੂਹ_ldb;
      ਦਲੀਲਾਂ_array[i].ldc=ਸਮੂਹ_ldc;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ COMPLEX
      ਦਲੀਲਾਂ_array[i].alpha=alpha;
      ਦਲੀਲਾਂ_array[i].beta=beta;
#ਨਹੀਂਤਾਂ
      ਦਲੀਲਾਂ_array[i].alpha=&alpha;
      ਦਲੀਲਾਂ_array[i].beta=&beta;
#ਜੇ_ਅੰਤ

#ਜੇ defined(CBLAS)
      ਜੇ (ਕ੍ਰਮ == Cblas_ਖਾਨਾ_Major) {      
	ਦਲੀਲਾਂ_array[i].a=&(a[i*stridea]);
	ਦਲੀਲਾਂ_array[i].b=&(b[i*strideb]);
      }ਨਹੀਂਤਾਂ ਜੇ(ਕ੍ਰਮ == Cblas_ਪੰਕਤੀ_Major){
	ਦਲੀਲਾਂ_array[i].a=&(b[i*strideb]);
	ਦਲੀਲਾਂ_array[i].b=&(a[i*stridea]);
      }
#ਨਹੀਂਤਾਂ
	ਦਲੀਲਾਂ_array[i].a=&(a[i*stridea]);
	ਦਲੀਲਾਂ_array[i].b=&(b[i*strideb]);
#ਜੇ_ਅੰਤ
      
      ਦਲੀਲਾਂ_array[i].c= &c[i*stridec];
      
      ਦਲੀਲਾਂ_array[i].routine_ਬਹੁਲਕ=ਸਮੂਹ_ਬਹੁਲਕ;
      ਦਲੀਲਾਂ_array[i].routine=ਸਮੂਹ_routine;
#ਜੇ_ਪਰਿਭਾਸ਼ਿਤ SMALL_ਮੈਟ੍ਰਿਕਸ_OPT
      ਜੇ (!ਸਮੂਹ_routine)
        ਦਲੀਲਾਂ_array[i].routine=ਸਮੂਹ_small_ਮੈਟ੍ਰਿਕਸ_opt_routine;
#ਜੇ_ਅੰਤ      
  }

  ਜੇ(ਗਿਣਤੀ_2>0) {
    GEMM_BATCH_THREAD(ਦਲੀਲਾਂ_array,ਗਿਣਤੀ_2);
  }

  ਮੁਕਤ(ਦਲੀਲਾਂ_array);
}
