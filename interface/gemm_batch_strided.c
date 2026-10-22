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

تهی openblas_warning(صحیح verbose, ثابت نویسه * پیام);

#ifndef COMPLEX
#ifdef XDOUBLE
#define ERROR_NAME "QGEMM_BATCH_STRIDED "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "DGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD dgemm_batch_thread
#وگرنه
#define ERROR_NAME "SGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD sgemm_batch_thread
#endif
#وگرنه
#ifdef XDOUBLE
#define ERROR_NAME "XGEMM_BATCH_STRIDED "
#وگرنهاگر defined(DOUBLE)
#define ERROR_NAME "ZGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD zgemm_batch_thread
#وگرنه
#define ERROR_NAME "CGEMM_BATCH_STRIDED "
#define GEMM_BATCH_THREAD cgemm_batch_thread
#endif
#endif
ایستا صحیح (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#اگر defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#define USE_SMALL_MATRIX_OPT 1
#وگرنه
#define USE_SMALL_MATRIX_OPT 0
#endif

#اگر USE_SMALL_MATRIX_OPT
#ifndef DYNAMIC_ARCH
#define SMALL_KERNEL_ADDR(table, idx) ((تهی *)(table[idx]))
#وگرنه
#define SMALL_KERNEL_ADDR(table, idx) ((تهی *)(*(uintptr_t *)((نویسه *)gotoblas + (size_t)(table[idx]))))
#endif


#ifndef COMPLEX
ایستا size_t gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


ایستا size_t gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#define GEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#define GEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#وگرنه

ایستا size_t zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

ایستا size_t zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#define ZGEMM_SMALL_KERNEL(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#define ZGEMM_SMALL_KERNEL_B0(idx) (صحیح (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
تهی NAME(نویسه *transa, نویسه *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           نویسه ta = *transa;
   	   نویسه tb = *transb;	   
	   blasint شمار = *matcount;
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
#وگرنه
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#endif
#وگرنه

تهی CNAME(شمارش CBLAS_ORDER order, شمارش CBLAS_TRANSPOSE transa, شمارش CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#ifndef COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint شمار) {
#وگرنه
	   تهی * valpha,
	   تهی * va, blasint lda, blasint stridea,
	   تهی * vb, blasint ldb, blasint strideb,
	   تهی * vbeta,
	   تهی * vc, blasint ldc, blasint stridec, blasint شمار) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#endif
#endif
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=پوچ;

  صحیح mode=0, group_mode=0;

  blasint i=0;

  صحیح group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  تهی * group_routine=پوچ;
#ifdef SMALL_MATRIX_OPT
  تهی * group_small_matrix_opt_routine=پوچ;
#endif
  
#اگر defined (SMP) || defined(SMALL_MATRIX_OPT)
  دوگانه MNK;
#endif

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)malloc(شمار * اندازه(blas_arg_t));
  
  اگر(args_array == پوچ){
    openblas_warning(0, "memory alloc failed!\n");
    بازگشت;
  }

#ifdef SMP
#ifndef COMPLEX
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#وگرنهاگر defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#وگرنه
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#endif
#وگرنه
#ifdef XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#وگرنهاگر defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#وگرنه
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#endif
#endif
#endif
  
  برای(i=0; i<شمار; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#اگر defined(CBLAS)
    اگر (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      اگر (transa == CblasNoTrans)     group_transa = 0;
      اگر (transa == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      اگر (transa == CblasConjNoTrans) group_transa = 0;
      اگر (transa == CblasConjTrans)   group_transa = 1;
#وگرنه
      اگر (transa == CblasConjNoTrans) group_transa = 2;
      اگر (transa == CblasConjTrans)   group_transa = 3;
#endif
      اگر (transb == CblasNoTrans)     group_transb = 0;
      اگر (transb == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      اگر (transb == CblasConjNoTrans) group_transb = 0;
      اگر (transb == CblasConjTrans)   group_transb = 1;
#وگرنه
      اگر (transb == CblasConjNoTrans) group_transb = 2;
      اگر (transb == CblasConjTrans)   group_transb = 3;
#endif
      group_nrowa = group_m;
      اگر (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      اگر (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      اگر (group_ldc < group_m) info = 13;
      اگر (group_ldb < group_nrowb)  info = 10;
      اگر (group_lda < group_nrowa)  info =  8;
      اگر (group_k < 0)        info =  5;
      اگر (group_n < 0)        info =  4;
      اگر (group_m < 0)        info =  3;
      اگر (group_transb < 0)        info =  2;
      اگر (group_transa < 0)        info =  1;

    }وگرنه اگر (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      اگر (transb == CblasNoTrans)     group_transa = 0;
      اگر (transb == CblasTrans)       group_transa = 1;
#ifndef COMPLEX
      اگر (transb == CblasConjNoTrans) group_transa = 0;
      اگر (transb == CblasConjTrans)   group_transa = 1;
#وگرنه
      اگر (transb == CblasConjNoTrans) group_transa = 2;
      اگر (transb == CblasConjTrans)   group_transa = 3;
#endif
      اگر (transa == CblasNoTrans)     group_transb = 0;
      اگر (transa == CblasTrans)       group_transb = 1;
#ifndef COMPLEX
      اگر (transa == CblasConjNoTrans) group_transb = 0;
      اگر (transa == CblasConjTrans)   group_transb = 1;
#وگرنه
      اگر (transa == CblasConjNoTrans) group_transb = 2;
      اگر (transa == CblasConjTrans)   group_transb = 3;
#endif

#وگرنه
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      اگر (ta == 'N')   group_transa = 0;
      اگر (ta == 'T')   group_transa = 1;
#ifndef COMPLEX
      اگر (ta == 'C')   group_transa = 1;
#وگرنه
      اگر (ta == 'C')   group_transa = 3;
#endif
      اگر (tb == 'N')   group_transb = 0;
      اگر (tb == 'T')   group_transb = 1;
#ifndef COMPLEX
      اگر (tb == 'C')   group_transb = 1;
#وگرنه
      اگر (tb == 'C')   group_transb = 3;
#endif
#endif

      group_nrowa = group_m;
      اگر (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      اگر (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      اگر (group_ldc < group_m) info = 13;
      اگر (group_ldb < group_nrowb)  info = 10;
      اگر (group_lda < group_nrowa)  info =  8;
      اگر (group_k < 0)        info =  5;
      اگر (group_n < 0)        info =  4;
      اگر (group_m < 0)        info =  3;
      اگر (group_transb < 0)        info =  2;
      اگر (group_transa < 0)        info =  1;      
#اگر defined(CBLAS)
    }
#endif

    اگر (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, اندازه(ERROR_NAME));
      آزاد(args_array);
      بازگشت;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    اگر (group_m == 0 || group_n == 0) { آزاد(args_array); بازگشت; }

    group_mode=mode;

#اگر defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (دوگانه) group_m * (دوگانه) group_n * (دوگانه) group_k;
#endif

#ifdef SMALL_MATRIX_OPT
#اگر !defined(COMPLEX)
    اگر (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#وگرنه
    اگر (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#endif
      group_routine=پوچ;
#اگر !defined(COMPLEX)
      اگر(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }وگرنه{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#وگرنه
      اگر(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }وگرنه{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#endif
      
    }وگرنه{
#endif
      group_routine=(تهی*)(gemm[(group_transb<<2)|group_transa]);
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
#وگرنه
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#endif

#اگر defined(CBLAS)
      اگر (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }وگرنه اگر(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#وگرنه
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#endif
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#ifdef SMALL_MATRIX_OPT
      اگر (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#endif      
  }

  اگر(شمار>0) {
    GEMM_BATCH_THREAD(args_array,شمار);
  }

  آزاد(args_array);
}
