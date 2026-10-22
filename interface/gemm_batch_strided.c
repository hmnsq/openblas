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

#అన్తర్భావ <stdio.h>
#అన్తర్భావ <stdlib.h>
#అన్తర్భావ "common.h"

శూన్య openblas_warning(పూర్ణాంకం verbose, స్థిరరాశి అక్షరం * సందేశం);

#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "QGEMM_BATCH_STRIDED "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "DGEMM_BATCH_STRIDED "
#నిర్వచించు GEMM_BATCH_THREAD dgemm_batch_thread
#లేకపోతే
#నిర్వచించు ERROR_NAME "SGEMM_BATCH_STRIDED "
#నిర్వచించు GEMM_BATCH_THREAD sgemm_batch_thread
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
#నిర్వచించు ERROR_NAME "XGEMM_BATCH_STRIDED "
#లేదావేళ defined(DOUBLE)
#నిర్వచించు ERROR_NAME "ZGEMM_BATCH_STRIDED "
#నిర్వచించు GEMM_BATCH_THREAD zgemm_batch_thread
#లేకపోతే
#నిర్వచించు ERROR_NAME "CGEMM_BATCH_STRIDED "
#నిర్వచించు GEMM_BATCH_THREAD cgemm_batch_thread
#అయితే_ముగింపు
#అయితే_ముగింపు
స్థిర పూర్ణాంకం (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#ఒకవేళ defined(SMALL_మాతృక_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#నిర్వచించు USE_SMALL_మాతృక_OPT 1
#లేకపోతే
#నిర్వచించు USE_SMALL_మాతృక_OPT 0
#అయితే_ముగింపు

#ఒకవేళ USE_SMALL_మాతృక_OPT
#నిర్వచితం_కాకపోతే DYNAMIC_ARCH
#నిర్వచించు SMALL_KERNEL_ADDR(పట్టిక, idx) ((శూన్య *)(పట్టిక[idx]))
#లేకపోతే
#నిర్వచించు SMALL_KERNEL_ADDR(పట్టిక, idx) ((శూన్య *)(*(uintptr_t *)((అక్షరం *)gotoblas + (size_t)(పట్టిక[idx]))))
#అయితే_ముగింపు


#నిర్వచితం_కాకపోతే COMPLEX
స్థిర size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


స్థిర size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#నిర్వచించు GEMM_SMALL_KERNEL_B0(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#నిర్వచించు GEMM_SMALL_KERNEL(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#లేకపోతే

స్థిర size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

స్థిర size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#నిర్వచించు ZGEMM_SMALL_KERNEL(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#నిర్వచించు ZGEMM_SMALL_KERNEL_B0(idx) (పూర్ణాంకం (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#అయితే_ముగింపు
#అయితే_ముగింపు

OPENBLAS_EXPORT

#నిర్వచితం_కాకపోతే CBLAS
శూన్య NAME(అక్షరం *transa, అక్షరం *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           అక్షరం ta = *transa;
   	   అక్షరం tb = *transb;	   
	   blasint లెక్క = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#ఒకవేళ !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#లేకపోతే
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#అయితే_ముగింపు
#లేకపోతే

శూన్య CNAME(గణన CBLAS_క్రమం క్రమం, గణన CBLAS_TRANSPOSE transa, గణన CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#నిర్వచితం_కాకపోతే COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint లెక్క) {
#లేకపోతే
	   శూన్య * valpha,
	   శూన్య * va, blasint lda, blasint stridea,
	   శూన్య * vb, blasint ldb, blasint strideb,
	   శూన్య * vbeta,
	   శూన్య * vc, blasint ldc, blasint stridec, blasint లెక్క) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#అయితే_ముగింపు
#అయితే_ముగింపు
  BLASLONG సమూహం_m, సమూహం_n, సమూహం_k;
  BLASLONG సమూహం_lda, సమూహం_ldb, సమూహం_ldc;

  blas_arg_t * వాదనలు_array=ఖాళీ;

  పూర్ణాంకం బహుళకం=0, సమూహం_బహుళకం=0;

  blasint i=0;

  పూర్ణాంకం సమూహం_transa, సమూహం_transb;
  BLASLONG సమూహం_nrowa, సమూహం_nrowb;
  blasint info;

  శూన్య * సమూహం_routine=ఖాళీ;
#నిర్వచితమైతే SMALL_మాతృక_OPT
  శూన్య * సమూహం_small_మాతృక_opt_routine=ఖాళీ;
#అయితే_ముగింపు
  
#ఒకవేళ defined (SMP) || defined(SMALL_మాతృక_OPT)
  ద్విగుణ MNK;
#అయితే_ముగింపు

  PRINT_DEBUG_CNAME;

  వాదనలు_array=(blas_arg_t *)జ్ఞాపకం_ఇవ్వు(లెక్క * పరిమాణం(blas_arg_t));
  
  ఒకవేళ(వాదనలు_array == ఖాళీ){
    openblas_warning(0, "memory alloc failed!\n");
    తిరిగి;
  }

#నిర్వచితమైతే SMP
#నిర్వచితం_కాకపోతే COMPLEX
#నిర్వచితమైతే XDOUBLE
  బహుళకం  =  BLAS_XDOUBLE | BLAS_REAL;
#లేదావేళ defined(DOUBLE)
  బహుళకం  =  BLAS_DOUBLE  | BLAS_REAL;
#లేకపోతే
  బహుళకం  =  BLAS_SINGLE  | BLAS_REAL;
#అయితే_ముగింపు
#లేకపోతే
#నిర్వచితమైతే XDOUBLE
  బహుళకం  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#లేదావేళ defined(DOUBLE)
  బహుళకం  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#లేకపోతే
  బహుళకం  =  BLAS_SINGLE  | BLAS_COMPLEX;
#అయితే_ముగింపు
#అయితే_ముగింపు
#అయితే_ముగింపు
  
  కొరకు(i=0; i<లెక్క; i++) {

    సమూహం_m = సమూహం_n = సమూహం_k = 0;
    సమూహం_lda = సమూహం_ldb = సమూహం_ldc = 0;
    సమూహం_transa = -1;
    సమూహం_transb = -1;
    info   = 0;
   
#ఒకవేళ defined(CBLAS)
    ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {
      సమూహం_m = m;
      సమూహం_n = n;
      సమూహం_k = k;

      సమూహం_lda = lda;
      సమూహం_ldb = ldb;
      సమూహం_ldc = ldc;

      ఒకవేళ (transa == CblasNoTrans)     సమూహం_transa = 0;
      ఒకవేళ (transa == CblasTrans)       సమూహం_transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (transa == CblasConjNoTrans) సమూహం_transa = 0;
      ఒకవేళ (transa == CblasConjTrans)   సమూహం_transa = 1;
#లేకపోతే
      ఒకవేళ (transa == CblasConjNoTrans) సమూహం_transa = 2;
      ఒకవేళ (transa == CblasConjTrans)   సమూహం_transa = 3;
#అయితే_ముగింపు
      ఒకవేళ (transb == CblasNoTrans)     సమూహం_transb = 0;
      ఒకవేళ (transb == CblasTrans)       సమూహం_transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (transb == CblasConjNoTrans) సమూహం_transb = 0;
      ఒకవేళ (transb == CblasConjTrans)   సమూహం_transb = 1;
#లేకపోతే
      ఒకవేళ (transb == CblasConjNoTrans) సమూహం_transb = 2;
      ఒకవేళ (transb == CblasConjTrans)   సమూహం_transb = 3;
#అయితే_ముగింపు
      సమూహం_nrowa = సమూహం_m;
      ఒకవేళ (సమూహం_transa & 1) సమూహం_nrowa = సమూహం_k;
      సమూహం_nrowb = సమూహం_k;
      ఒకవేళ (సమూహం_transb & 1) సమూహం_nrowb = సమూహం_n;

      info=-1;
	
      ఒకవేళ (సమూహం_ldc < సమూహం_m) info = 13;
      ఒకవేళ (సమూహం_ldb < సమూహం_nrowb)  info = 10;
      ఒకవేళ (సమూహం_lda < సమూహం_nrowa)  info =  8;
      ఒకవేళ (సమూహం_k < 0)        info =  5;
      ఒకవేళ (సమూహం_n < 0)        info =  4;
      ఒకవేళ (సమూహం_m < 0)        info =  3;
      ఒకవేళ (సమూహం_transb < 0)        info =  2;
      ఒకవేళ (సమూహం_transa < 0)        info =  1;

    }లేకపోతే ఒకవేళ (క్రమం == Cblas_అడ్డువరుస_Major) {
      
      సమూహం_m = n;
      సమూహం_n = m;
      సమూహం_k = k;

      సమూహం_lda = ldb;
      సమూహం_ldb = lda;
      సమూహం_ldc = ldc;

      ఒకవేళ (transb == CblasNoTrans)     సమూహం_transa = 0;
      ఒకవేళ (transb == CblasTrans)       సమూహం_transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (transb == CblasConjNoTrans) సమూహం_transa = 0;
      ఒకవేళ (transb == CblasConjTrans)   సమూహం_transa = 1;
#లేకపోతే
      ఒకవేళ (transb == CblasConjNoTrans) సమూహం_transa = 2;
      ఒకవేళ (transb == CblasConjTrans)   సమూహం_transa = 3;
#అయితే_ముగింపు
      ఒకవేళ (transa == CblasNoTrans)     సమూహం_transb = 0;
      ఒకవేళ (transa == CblasTrans)       సమూహం_transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (transa == CblasConjNoTrans) సమూహం_transb = 0;
      ఒకవేళ (transa == CblasConjTrans)   సమూహం_transb = 1;
#లేకపోతే
      ఒకవేళ (transa == CblasConjNoTrans) సమూహం_transb = 2;
      ఒకవేళ (transa == CblasConjTrans)   సమూహం_transb = 3;
#అయితే_ముగింపు

#లేకపోతే
      సమూహం_m = m;
      సమూహం_n = n;
      సమూహం_k = k;

      సమూహం_lda = lda;
      సమూహం_ldb = ldb;
      సమూహం_ldc = ldc;
      
      ఒకవేళ (ta == 'N')   సమూహం_transa = 0;
      ఒకవేళ (ta == 'T')   సమూహం_transa = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (ta == 'C')   సమూహం_transa = 1;
#లేకపోతే
      ఒకవేళ (ta == 'C')   సమూహం_transa = 3;
#అయితే_ముగింపు
      ఒకవేళ (tb == 'N')   సమూహం_transb = 0;
      ఒకవేళ (tb == 'T')   సమూహం_transb = 1;
#నిర్వచితం_కాకపోతే COMPLEX
      ఒకవేళ (tb == 'C')   సమూహం_transb = 1;
#లేకపోతే
      ఒకవేళ (tb == 'C')   సమూహం_transb = 3;
#అయితే_ముగింపు
#అయితే_ముగింపు

      సమూహం_nrowa = సమూహం_m;
      ఒకవేళ (సమూహం_transa & 1) సమూహం_nrowa = సమూహం_k;
      సమూహం_nrowb = సమూహం_k;
      ఒకవేళ (సమూహం_transb & 1) సమూహం_nrowb = సమూహం_n;

      info=-1;
	
      ఒకవేళ (సమూహం_ldc < సమూహం_m) info = 13;
      ఒకవేళ (సమూహం_ldb < సమూహం_nrowb)  info = 10;
      ఒకవేళ (సమూహం_lda < సమూహం_nrowa)  info =  8;
      ఒకవేళ (సమూహం_k < 0)        info =  5;
      ఒకవేళ (సమూహం_n < 0)        info =  4;
      ఒకవేళ (సమూహం_m < 0)        info =  3;
      ఒకవేళ (సమూహం_transb < 0)        info =  2;
      ఒకవేళ (సమూహం_transa < 0)        info =  1;      
#ఒకవేళ defined(CBLAS)
    }
#అయితే_ముగింపు

    ఒకవేళ (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, పరిమాణం(ERROR_NAME));
      విడుదల(వాదనలు_array);
      తిరిగి;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    ఒకవేళ (సమూహం_m == 0 || సమూహం_n == 0) { విడుదల(వాదనలు_array); తిరిగి; }

    సమూహం_బహుళకం=బహుళకం;

#ఒకవేళ defined(SMP) || defined(SMALL_మాతృక_OPT)
    MNK = (ద్విగుణ) సమూహం_m * (ద్విగుణ) సమూహం_n * (ద్విగుణ) సమూహం_k;
#అయితే_ముగింపు

#నిర్వచితమైతే SMALL_మాతృక_OPT
#ఒకవేళ !defined(COMPLEX)
    ఒకవేళ (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_మాతృక_PERMIT(సమూహం_transa, సమూహం_transb, సమూహం_m, సమూహం_n, సమూహం_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#లేకపోతే
    ఒకవేళ (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_మాతృక_PERMIT(సమూహం_transa, సమూహం_transb, సమూహం_m, సమూహం_n, సమూహం_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#అయితే_ముగింపు
      సమూహం_routine=ఖాళీ;
#ఒకవేళ !defined(COMPLEX)
      ఒకవేళ(beta == 0.0){
	సమూహం_బహుళకం=బహుళకం | BLAS_SMALL_B0_OPT;
	సమూహం_small_మాతృక_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((సమూహం_transb<<2)|సమూహం_transa));
      }లేకపోతే{
	సమూహం_బహుళకం=బహుళకం | BLAS_SMALL_OPT;
	సమూహం_small_మాతృక_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((సమూహం_transb<<2)|సమూహం_transa));
      }
#లేకపోతే
      ఒకవేళ(beta[0] == 0.0 && beta[1] == 0.0){
	సమూహం_బహుళకం=బహుళకం | BLAS_SMALL_B0_OPT;
	సమూహం_small_మాతృక_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((సమూహం_transb<<2)|సమూహం_transa));
      }లేకపోతే{
	సమూహం_బహుళకం=బహుళకం | BLAS_SMALL_OPT;
	సమూహం_small_మాతృక_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((సమూహం_transb<<2)|సమూహం_transa));
      }

#అయితే_ముగింపు
      
    }లేకపోతే{
#అయితే_ముగింపు
      సమూహం_routine=(శూన్య*)(gemm[(సమూహం_transb<<2)|సమూహం_transa]);
#నిర్వచితమైతే SMALL_మాతృక_OPT
    }
#అయితే_ముగింపు    

    
      వాదనలు_array[i].m=సమూహం_m;
      వాదనలు_array[i].n=సమూహం_n;
      వాదనలు_array[i].k=సమూహం_k;
      వాదనలు_array[i].lda=సమూహం_lda;
      వాదనలు_array[i].ldb=సమూహం_ldb;
      వాదనలు_array[i].ldc=సమూహం_ldc;
#నిర్వచితమైతే COMPLEX
      వాదనలు_array[i].alpha=alpha;
      వాదనలు_array[i].beta=beta;
#లేకపోతే
      వాదనలు_array[i].alpha=&alpha;
      వాదనలు_array[i].beta=&beta;
#అయితే_ముగింపు

#ఒకవేళ defined(CBLAS)
      ఒకవేళ (క్రమం == Cblas_నిలువు_Major) {      
	వాదనలు_array[i].a=&(a[i*stridea]);
	వాదనలు_array[i].b=&(b[i*strideb]);
      }లేకపోతే ఒకవేళ(క్రమం == Cblas_అడ్డువరుస_Major){
	వాదనలు_array[i].a=&(b[i*strideb]);
	వాదనలు_array[i].b=&(a[i*stridea]);
      }
#లేకపోతే
	వాదనలు_array[i].a=&(a[i*stridea]);
	వాదనలు_array[i].b=&(b[i*strideb]);
#అయితే_ముగింపు
      
      వాదనలు_array[i].c= &c[i*stridec];
      
      వాదనలు_array[i].routine_బహుళకం=సమూహం_బహుళకం;
      వాదనలు_array[i].routine=సమూహం_routine;
#నిర్వచితమైతే SMALL_మాతృక_OPT
      ఒకవేళ (!సమూహం_routine)
        వాదనలు_array[i].routine=సమూహం_small_మాతృక_opt_routine;
#అయితే_ముగింపు      
  }

  ఒకవేళ(లెక్క>0) {
    GEMM_BATCH_THREAD(వాదనలు_array,లెక్క);
  }

  విడుదల(వాదనలు_array);
}
