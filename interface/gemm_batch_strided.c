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

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव "common.h"

शून्य openblas_warning(पूर्णाङ्क verbose, ध्रुव वर्ण * सन्देशः);

#यदि_अपरिभाषितम् COMPLEX
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "QGEMM_BATCH_STRIDED "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "DGEMM_BATCH_STRIDED "
#परिभाषय GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#परिभाषय ERROR_NAME "SGEMM_BATCH_STRIDED "
#परिभाषय GEMM_BATCH_THREAD sgemm_batch_thread
#यदि_अन्तः
#अन्यथा
#यदि_परिभाषितम् XDOUBLE
#परिभाषय ERROR_NAME "XGEMM_BATCH_STRIDED "
#नोचेद्यदि defined(DOUBLE)
#परिभाषय ERROR_NAME "ZGEMM_BATCH_STRIDED "
#परिभाषय GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#परिभाषय ERROR_NAME "CGEMM_BATCH_STRIDED "
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
शून्य NAME(वर्ण *transa, वर्ण *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           वर्ण ta = *transa;
   	   वर्ण tb = *transb;	   
	   blasint गणनम् = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#यदि !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#अन्यथा
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#यदि_अन्तः
#अन्यथा

शून्य CNAME(गणना CBLAS_क्रमः क्रमः, गणना CBLAS_TRANSPOSE transa, गणना CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#यदि_अपरिभाषितम् COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint गणनम्) {
#अन्यथा
	   शून्य * valpha,
	   शून्य * va, blasint lda, blasint stridea,
	   शून्य * vb, blasint ldb, blasint strideb,
	   शून्य * vbeta,
	   शून्य * vc, blasint ldc, blasint stridec, blasint गणनम्) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#यदि_अन्तः
#यदि_अन्तः
  BLASLONG गणः_m, गणः_n, गणः_k;
  BLASLONG गणः_lda, गणः_ldb, गणः_ldc;

  blas_arg_t * तर्काः_array=रिक्त;

  पूर्णाङ्क बहुलकम्=0, गणः_बहुलकम्=0;

  blasint i=0;

  पूर्णाङ्क गणः_transa, गणः_transb;
  BLASLONG गणः_nrowa, गणः_nrowb;
  blasint info;

  शून्य * गणः_routine=रिक्त;
#यदि_परिभाषितम् SMALL_MATRIX_OPT
  शून्य * गणः_small_matrix_opt_routine=रिक्त;
#यदि_अन्तः
  
#यदि defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुण MNK;
#यदि_अन्तः

  PRINT_DEBUG_CNAME;

  तर्काः_array=(blas_arg_t *)स्मृतिं_देहि(गणनम् * परिमाण(blas_arg_t));
  
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
  
  प्रति(i=0; i<गणनम्; i++) {

    गणः_m = गणः_n = गणः_k = 0;
    गणः_lda = गणः_ldb = गणः_ldc = 0;
    गणः_transa = -1;
    गणः_transb = -1;
    info   = 0;
   
#यदि defined(CBLAS)
    यदि (क्रमः == CblasColMajor) {
      गणः_m = m;
      गणः_n = n;
      गणः_k = k;

      गणः_lda = lda;
      गणः_ldb = ldb;
      गणः_ldc = ldc;

      यदि (transa == CblasNoTrans)     गणः_transa = 0;
      यदि (transa == CblasTrans)       गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa == CblasConjNoTrans) गणः_transa = 0;
      यदि (transa == CblasConjTrans)   गणः_transa = 1;
#अन्यथा
      यदि (transa == CblasConjNoTrans) गणः_transa = 2;
      यदि (transa == CblasConjTrans)   गणः_transa = 3;
#यदि_अन्तः
      यदि (transb == CblasNoTrans)     गणः_transb = 0;
      यदि (transb == CblasTrans)       गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transb == CblasConjNoTrans) गणः_transb = 0;
      यदि (transb == CblasConjTrans)   गणः_transb = 1;
#अन्यथा
      यदि (transb == CblasConjNoTrans) गणः_transb = 2;
      यदि (transb == CblasConjTrans)   गणः_transb = 3;
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
      
      गणः_m = n;
      गणः_n = m;
      गणः_k = k;

      गणः_lda = ldb;
      गणः_ldb = lda;
      गणः_ldc = ldc;

      यदि (transb == CblasNoTrans)     गणः_transa = 0;
      यदि (transb == CblasTrans)       गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transb == CblasConjNoTrans) गणः_transa = 0;
      यदि (transb == CblasConjTrans)   गणः_transa = 1;
#अन्यथा
      यदि (transb == CblasConjNoTrans) गणः_transa = 2;
      यदि (transb == CblasConjTrans)   गणः_transa = 3;
#यदि_अन्तः
      यदि (transa == CblasNoTrans)     गणः_transb = 0;
      यदि (transa == CblasTrans)       गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (transa == CblasConjNoTrans) गणः_transb = 0;
      यदि (transa == CblasConjTrans)   गणः_transb = 1;
#अन्यथा
      यदि (transa == CblasConjNoTrans) गणः_transb = 2;
      यदि (transa == CblasConjTrans)   गणः_transb = 3;
#यदि_अन्तः

#अन्यथा
      गणः_m = m;
      गणः_n = n;
      गणः_k = k;

      गणः_lda = lda;
      गणः_ldb = ldb;
      गणः_ldc = ldc;
      
      यदि (ta == 'N')   गणः_transa = 0;
      यदि (ta == 'T')   गणः_transa = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (ta == 'C')   गणः_transa = 1;
#अन्यथा
      यदि (ta == 'C')   गणः_transa = 3;
#यदि_अन्तः
      यदि (tb == 'N')   गणः_transb = 0;
      यदि (tb == 'T')   गणः_transb = 1;
#यदि_अपरिभाषितम् COMPLEX
      यदि (tb == 'C')   गणः_transb = 1;
#अन्यथा
      यदि (tb == 'C')   गणः_transb = 3;
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

    /* Every member shares M and N: an empty member means an empty batch. */
    यदि (गणः_m == 0 || गणः_n == 0) { मुक्त(तर्काः_array); प्रत्यागम; }

    गणः_बहुलकम्=बहुलकम्;

#यदि defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुण) गणः_m * (द्विगुण) गणः_n * (द्विगुण) गणः_k;
#यदि_अन्तः

#यदि_परिभाषितम् SMALL_MATRIX_OPT
#यदि !defined(COMPLEX)
    यदि (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(गणः_transa, गणः_transb, गणः_m, गणः_n, गणः_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#अन्यथा
    यदि (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(गणः_transa, गणः_transb, गणः_m, गणः_n, गणः_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#यदि_अन्तः
      गणः_routine=रिक्त;
#यदि !defined(COMPLEX)
      यदि(beta == 0.0){
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_B0_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((गणः_transb<<2)|गणः_transa));
      }अन्यथा{
	गणः_बहुलकम्=बहुलकम् | BLAS_SMALL_OPT;
	गणः_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((गणः_transb<<2)|गणः_transa));
      }
#अन्यथा
      यदि(beta[0] == 0.0 && beta[1] == 0.0){
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

    
      तर्काः_array[i].m=गणः_m;
      तर्काः_array[i].n=गणः_n;
      तर्काः_array[i].k=गणः_k;
      तर्काः_array[i].lda=गणः_lda;
      तर्काः_array[i].ldb=गणः_ldb;
      तर्काः_array[i].ldc=गणः_ldc;
#यदि_परिभाषितम् COMPLEX
      तर्काः_array[i].alpha=alpha;
      तर्काः_array[i].beta=beta;
#अन्यथा
      तर्काः_array[i].alpha=&alpha;
      तर्काः_array[i].beta=&beta;
#यदि_अन्तः

#यदि defined(CBLAS)
      यदि (क्रमः == CblasColMajor) {      
	तर्काः_array[i].a=&(a[i*stridea]);
	तर्काः_array[i].b=&(b[i*strideb]);
      }अन्यथा यदि(क्रमः == CblasRowMajor){
	तर्काः_array[i].a=&(b[i*strideb]);
	तर्काः_array[i].b=&(a[i*stridea]);
      }
#अन्यथा
	तर्काः_array[i].a=&(a[i*stridea]);
	तर्काः_array[i].b=&(b[i*strideb]);
#यदि_अन्तः
      
      तर्काः_array[i].c= &c[i*stridec];
      
      तर्काः_array[i].routine_बहुलकम्=गणः_बहुलकम्;
      तर्काः_array[i].routine=गणः_routine;
#यदि_परिभाषितम् SMALL_MATRIX_OPT
      यदि (!गणः_routine)
        तर्काः_array[i].routine=गणः_small_matrix_opt_routine;
#यदि_अन्तः      
  }

  यदि(गणनम्>0) {
    GEMM_BATCH_THREAD(तर्काः_array,गणनम्);
  }

  मुक्त(तर्काः_array);
}
