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

शून्य openblas_warning(पूर्णाङ्क verbose, ध्रुव वर्ण * संदेश);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#अन्यथा
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#नोचेद्यदि defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
स्थिर पूर्णाङ्क (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#यदि defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#अन्यथा
#define USE_SMALL_MATRIX_OPT 0
#endif

#यदि USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(table[idx]))
#अन्यथा
#define SMALL_KERNEL_ADDR(table, idx) ((शून्य *)(*(uintptr_t *)((वर्ण *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
स्थिर size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


स्थिर size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
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

#define ZGEMM_SMALL_KERNEL(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णाङ्क (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
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
	   blasint मोजणी = *matcount;
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
#endif
#अन्यथा

शून्य CNAME(गणना CBLAS_ORDER order, गणना CBLAS_TRANSPOSE transa, गणना CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint मोजणी) {
#अन्यथा
	   शून्य * valpha,
	   शून्य * va, blasint lda, blasint stridea,
	   शून्य * vb, blasint ldb, blasint strideb,
	   शून्य * vbeta,
	   शून्य * vc, blasint ldc, blasint stridec, blasint मोजणी) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=रिक्त;

  पूर्णाङ्क mode=0, group_mode=0;

  blasint i=0;

  पूर्णाङ्क group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  शून्य * group_routine=रिक्त;
#ifdef SMALL_MATRIX_OPT
  शून्य * group_small_matrix_opt_routine=रिक्त;
#endif
  
#यदि defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुण MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(मोजणी * परिमाण(blas_arg_t));
  
  यदि(args_array == रिक्त){
    openblas_warning(0, "memory alloc failed!\n");
    प्रत्यागम;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#नोचेद्यदि defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#अन्यथा
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#नोचेद्यदि defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  प्रति(i=0; i<मोजणी; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#यदि defined(CBLAS)
    यदि (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      यदि (transa == CblasNoTrans)     group_transa = 0;
      यदि (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      यदि (transa == CblasConjNoTrans) group_transa = 0;
      यदि (transa == CblasConjTrans)   group_transa = 1;
#अन्यथा
      यदि (transa == CblasConjNoTrans) group_transa = 2;
      यदि (transa == CblasConjTrans)   group_transa = 3;
#endif
      यदि (transb == CblasNoTrans)     group_transb = 0;
      यदि (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      यदि (transb == CblasConjNoTrans) group_transb = 0;
      यदि (transb == CblasConjTrans)   group_transb = 1;
#अन्यथा
      यदि (transb == CblasConjNoTrans) group_transb = 2;
      यदि (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      यदि (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      यदि (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      यदि (group_ldc < group_m) info = 13;
      यदि (group_ldb < group_nrowb)  info = 10;
      यदि (group_lda < group_nrowa)  info =  8;
      यदि (group_k < 0)        info =  5;
      यदि (group_n < 0)        info =  4;
      यदि (group_m < 0)        info =  3;
      यदि (group_transb < 0)        info =  2;
      यदि (group_transa < 0)        info =  1;

    }अन्यथा यदि (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      यदि (transb == CblasNoTrans)     group_transa = 0;
      यदि (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      यदि (transb == CblasConjNoTrans) group_transa = 0;
      यदि (transb == CblasConjTrans)   group_transa = 1;
#अन्यथा
      यदि (transb == CblasConjNoTrans) group_transa = 2;
      यदि (transb == CblasConjTrans)   group_transa = 3;
#endif
      यदि (transa == CblasNoTrans)     group_transb = 0;
      यदि (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      यदि (transa == CblasConjNoTrans) group_transb = 0;
      यदि (transa == CblasConjTrans)   group_transb = 1;
#अन्यथा
      यदि (transa == CblasConjNoTrans) group_transb = 2;
      यदि (transa == CblasConjTrans)   group_transb = 3;
#endif

#अन्यथा
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      यदि (ta == 'N')   group_transa = 0;
      यदि (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      यदि (ta == 'C')   group_transa = 1;
#अन्यथा
      यदि (ta == 'C')   group_transa = 3;
#endif
      यदि (tb == 'N')   group_transb = 0;
      यदि (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      यदि (tb == 'C')   group_transb = 1;
#अन्यथा
      यदि (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      यदि (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      यदि (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      यदि (group_ldc < group_m) info = 13;
      यदि (group_ldb < group_nrowb)  info = 10;
      यदि (group_lda < group_nrowa)  info =  8;
      यदि (group_k < 0)        info =  5;
      यदि (group_n < 0)        info =  4;
      यदि (group_m < 0)        info =  3;
      यदि (group_transb < 0)        info =  2;
      यदि (group_transa < 0)        info =  1;      
#यदि defined(CBLAS)
    }
#endif

    यदि (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, परिमाण(ERROR_NAME));
      मुक्त(args_array);
      प्रत्यागम;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    यदि (group_m == 0 || group_n == 0) { मुक्त(args_array); प्रत्यागम; }

    group_mode=mode;

#यदि defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुण) group_m * (द्विगुण) group_n * (द्विगुण) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#यदि !defined(COMPLEX)
    यदि (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#अन्यथा
    यदि (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=रिक्त;
#यदि !defined(COMPLEX)
      यदि(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#अन्यथा
      यदि(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }अन्यथा{
#endif
      group_routine=(शून्य*)(gemm[(group_transb<<2)|group_transa]);
#ifdef SMALL_MATRIX_OPT
    }
#endif    

    
      args_array[i].m=group_m;
      args_array[i].n=group_n;
      args_array[i].k=group_k;
      args_array[i].lda=group_lda;
      args_array[i].ldb=group_ldb;
      args_array[i].ldc=group_ldc;
#ifdef COMPLEX
      args_array[i].alpha=alpha;
      args_array[i].beta=beta;
#अन्यथा
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#यदि defined(CBLAS)
      यदि (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }अन्यथा यदि(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#अन्यथा
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      यदि (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  यदि(मोजणी>0) {
    GEMM_BATCH_THREAD(args_array,मोजणी);
  }

  मुक्त(args_array);
}
