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

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश "common.h"

व्योम openblas_warning(पूर्णांक verbose, स्थिर अक्षर * संदेश);

#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "QGEMM_BATCH_STRIDED "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "DGEMM_BATCH_STRIDED "
#घोषणा GEMM_BATCH_THREAD dgemm_batch_thread
#अन्यथा
#घोषणा ERROR_NAME "SGEMM_BATCH_STRIDED "
#घोषणा GEMM_BATCH_THREAD sgemm_batch_thread
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
#घोषणा ERROR_NAME "XGEMM_BATCH_STRIDED "
#या_अगर defined(DOUBLE)
#घोषणा ERROR_NAME "ZGEMM_BATCH_STRIDED "
#घोषणा GEMM_BATCH_THREAD zgemm_batch_thread
#अन्यथा
#घोषणा ERROR_NAME "CGEMM_BATCH_STRIDED "
#घोषणा GEMM_BATCH_THREAD cgemm_batch_thread
#पूर्ण_अगर
#पूर्ण_अगर
अटल पूर्णांक (*gemm[])(blas_arg_t *, BLASLONG *, BLASLONG *, IFLOAT *, IFLOAT *, BLASLONG) = {
  GEMM_NN, GEMM_TN, GEMM_RN, GEMM_CN,
  GEMM_NT, GEMM_TT, GEMM_RT, GEMM_CT,
  GEMM_NR, GEMM_TR, GEMM_RR, GEMM_CR,
  GEMM_NC, GEMM_TC, GEMM_RC, GEMM_CC,
};

#अगर defined(SMALL_MATRIX_OPT) && !defined(GEMM3M) && !defined(XDOUBLE)
#घोषणा USE_SMALL_MATRIX_OPT 1
#अन्यथा
#घोषणा USE_SMALL_MATRIX_OPT 0
#पूर्ण_अगर

#अगर USE_SMALL_MATRIX_OPT
#अगर_अघोषित DYNAMIC_ARCH
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(table[idx]))
#अन्यथा
#घोषणा SMALL_KERNEL_ADDR(table, idx) ((व्योम *)(*(uintptr_t *)((अक्षर *)gotoblas + (माप_प्रकार)(table[idx]))))
#पूर्ण_अगर


#अगर_अघोषित COMPLEX
अटल माप_प्रकार gemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, 0, 0,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, 0, 0,
};


अटल माप_प्रकार gemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, 0, 0,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, 0, 0,
};

#घोषणा GEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel_b0, (idx))
#घोषणा GEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, IFLOAT *, BLASLONG, FLOAT, IFLOAT *, BLASLONG, FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(gemm_small_kernel, (idx))
#अन्यथा

अटल माप_प्रकार zgemm_small_kernel[] = {
	GEMM_SMALL_KERNEL_NN, GEMM_SMALL_KERNEL_TN, GEMM_SMALL_KERNEL_RN, GEMM_SMALL_KERNEL_CN,
	GEMM_SMALL_KERNEL_NT, GEMM_SMALL_KERNEL_TT, GEMM_SMALL_KERNEL_RT, GEMM_SMALL_KERNEL_CT,
	GEMM_SMALL_KERNEL_NR, GEMM_SMALL_KERNEL_TR, GEMM_SMALL_KERNEL_RR, GEMM_SMALL_KERNEL_CR,
	GEMM_SMALL_KERNEL_NC, GEMM_SMALL_KERNEL_TC, GEMM_SMALL_KERNEL_RC, GEMM_SMALL_KERNEL_CC,
};

अटल माप_प्रकार zgemm_small_kernel_b0[] = {
	GEMM_SMALL_KERNEL_B0_NN, GEMM_SMALL_KERNEL_B0_TN, GEMM_SMALL_KERNEL_B0_RN, GEMM_SMALL_KERNEL_B0_CN,
	GEMM_SMALL_KERNEL_B0_NT, GEMM_SMALL_KERNEL_B0_TT, GEMM_SMALL_KERNEL_B0_RT, GEMM_SMALL_KERNEL_B0_CT,
	GEMM_SMALL_KERNEL_B0_NR, GEMM_SMALL_KERNEL_B0_TR, GEMM_SMALL_KERNEL_B0_RR, GEMM_SMALL_KERNEL_B0_CR,
	GEMM_SMALL_KERNEL_B0_NC, GEMM_SMALL_KERNEL_B0_TC, GEMM_SMALL_KERNEL_B0_RC, GEMM_SMALL_KERNEL_B0_CC,
};

#घोषणा ZGEMM_SMALL_KERNEL(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel, (idx))
#घोषणा ZGEMM_SMALL_KERNEL_B0(idx) (पूर्णांक (*)(BLASLONG, BLASLONG, BLASLONG, FLOAT *, BLASLONG, FLOAT , FLOAT, FLOAT *, BLASLONG, FLOAT *, BLASLONG)) SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, (idx))
#पूर्ण_अगर
#पूर्ण_अगर

OPENBLAS_EXPORT

#अगर_अघोषित CBLAS
व्योम NAME(अक्षर *transa, अक्षर *transb,
	   blasint * M, blasint * N, blasint * K,
	   FLOAT * Alpha,
	   IFLOAT * a, blasint * Lda,
	   blasint * stride_a,
	   IFLOAT *b, blasint * Ldb,
	   blasint * stride_b,
	   FLOAT * Beta,
	   FLOAT * c, blasint * Ldc, blasint * stride_c, blasint * matcount) {
 	  
           अक्षर ta = *transa;
   	   अक्षर tb = *transb;	   
	   blasint गिनती = *matcount;
	   blasint stridea= *stride_a;
	   blasint strideb= *stride_b;
	   blasint stridec= *stride_c;
	   blasint m=*M;
	   blasint n=*N;
	   blasint k=*K;
	   blasint lda=*Lda;
	   blasint ldb=*Ldb;
	   blasint ldc=*Ldc;
#अगर !defined(COMPLEX)	   
	   FLOAT alpha=*Alpha;
	   FLOAT beta=*Beta;
#अन्यथा
	   FLOAT *alpha=Alpha;
	   FLOAT *beta=Beta;
#पूर्ण_अगर
#अन्यथा

व्योम CNAME(क्रमागत CBLAS_ORDER order, क्रमागत CBLAS_TRANSPOSE transa, क्रमागत CBLAS_TRANSPOSE transb,
	   blasint m, blasint n, blasint k,
#अगर_अघोषित COMPLEX
	   FLOAT alpha,
	   IFLOAT * a, blasint lda, blasint stridea,
	   IFLOAT * b, blasint ldb, blasint strideb,
	   FLOAT beta,
	   FLOAT * c, blasint ldc, blasint stridec, blasint गिनती) {
#अन्यथा
	   व्योम * valpha,
	   व्योम * va, blasint lda, blasint stridea,
	   व्योम * vb, blasint ldb, blasint strideb,
	   व्योम * vbeta,
	   व्योम * vc, blasint ldc, blasint stridec, blasint गिनती) {

  FLOAT * alpha=(FLOAT *)valpha;
  FLOAT * beta=(FLOAT *)vbeta;
  FLOAT * a=(FLOAT*)va;
  FLOAT * b=(FLOAT*)vb;
  FLOAT * c=(FLOAT*)vc;
#पूर्ण_अगर
#पूर्ण_अगर
  BLASLONG group_m, group_n, group_k;
  BLASLONG group_lda, group_ldb, group_ldc;

  blas_arg_t * args_array=शून्य;

  पूर्णांक mode=0, group_mode=0;

  blasint i=0;

  पूर्णांक group_transa, group_transb;
  BLASLONG group_nrowa, group_nrowb;
  blasint info;

  व्योम * group_routine=शून्य;
#अगर_घोषित SMALL_MATRIX_OPT
  व्योम * group_small_matrix_opt_routine=शून्य;
#पूर्ण_अगर
  
#अगर defined (SMP) || defined(SMALL_MATRIX_OPT)
  द्विगुन MNK;
#पूर्ण_अगर

  PRINT_DEBUG_CNAME;

  args_array=(blas_arg_t *)दो_स्मृति(गिनती * माप(blas_arg_t));
  
  अगर(args_array == शून्य){
    openblas_warning(0, "memory alloc failed!\n");
    वापस;
  }

#अगर_घोषित SMP
#अगर_अघोषित COMPLEX
#अगर_घोषित XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_REAL;
#या_अगर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_REAL;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_REAL;
#पूर्ण_अगर
#अन्यथा
#अगर_घोषित XDOUBLE
  mode  =  BLAS_XDOUBLE | BLAS_COMPLEX;
#या_अगर defined(DOUBLE)
  mode  =  BLAS_DOUBLE  | BLAS_COMPLEX;
#अन्यथा
  mode  =  BLAS_SINGLE  | BLAS_COMPLEX;
#पूर्ण_अगर
#पूर्ण_अगर
#पूर्ण_अगर
  
  क्रम(i=0; i<गिनती; i++) {

    group_m = group_n = group_k = 0;
    group_lda = group_ldb = group_ldc = 0;
    group_transa = -1;
    group_transb = -1;
    info   = 0;
   
#अगर defined(CBLAS)
    अगर (order == CblasColMajor) {
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;

      अगर (transa == CblasNoTrans)     group_transa = 0;
      अगर (transa == CblasTrans)       group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (transa == CblasConjNoTrans) group_transa = 0;
      अगर (transa == CblasConjTrans)   group_transa = 1;
#अन्यथा
      अगर (transa == CblasConjNoTrans) group_transa = 2;
      अगर (transa == CblasConjTrans)   group_transa = 3;
#पूर्ण_अगर
      अगर (transb == CblasNoTrans)     group_transb = 0;
      अगर (transb == CblasTrans)       group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (transb == CblasConjNoTrans) group_transb = 0;
      अगर (transb == CblasConjTrans)   group_transb = 1;
#अन्यथा
      अगर (transb == CblasConjNoTrans) group_transb = 2;
      अगर (transb == CblasConjTrans)   group_transb = 3;
#पूर्ण_अगर
      group_nrowa = group_m;
      अगर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      अगर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      अगर (group_ldc < group_m) info = 13;
      अगर (group_ldb < group_nrowb)  info = 10;
      अगर (group_lda < group_nrowa)  info =  8;
      अगर (group_k < 0)        info =  5;
      अगर (group_n < 0)        info =  4;
      अगर (group_m < 0)        info =  3;
      अगर (group_transb < 0)        info =  2;
      अगर (group_transa < 0)        info =  1;

    }अन्यथा अगर (order == CblasRowMajor) {
      
      group_m = n;
      group_n = m;
      group_k = k;

      group_lda = ldb;
      group_ldb = lda;
      group_ldc = ldc;

      अगर (transb == CblasNoTrans)     group_transa = 0;
      अगर (transb == CblasTrans)       group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (transb == CblasConjNoTrans) group_transa = 0;
      अगर (transb == CblasConjTrans)   group_transa = 1;
#अन्यथा
      अगर (transb == CblasConjNoTrans) group_transa = 2;
      अगर (transb == CblasConjTrans)   group_transa = 3;
#पूर्ण_अगर
      अगर (transa == CblasNoTrans)     group_transb = 0;
      अगर (transa == CblasTrans)       group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (transa == CblasConjNoTrans) group_transb = 0;
      अगर (transa == CblasConjTrans)   group_transb = 1;
#अन्यथा
      अगर (transa == CblasConjNoTrans) group_transb = 2;
      अगर (transa == CblasConjTrans)   group_transb = 3;
#पूर्ण_अगर

#अन्यथा
      group_m = m;
      group_n = n;
      group_k = k;

      group_lda = lda;
      group_ldb = ldb;
      group_ldc = ldc;
      
      अगर (ta == 'N')   group_transa = 0;
      अगर (ta == 'T')   group_transa = 1;
#अगर_अघोषित COMPLEX
      अगर (ta == 'C')   group_transa = 1;
#अन्यथा
      अगर (ta == 'C')   group_transa = 3;
#पूर्ण_अगर
      अगर (tb == 'N')   group_transb = 0;
      अगर (tb == 'T')   group_transb = 1;
#अगर_अघोषित COMPLEX
      अगर (tb == 'C')   group_transb = 1;
#अन्यथा
      अगर (tb == 'C')   group_transb = 3;
#पूर्ण_अगर
#पूर्ण_अगर

      group_nrowa = group_m;
      अगर (group_transa & 1) group_nrowa = group_k;
      group_nrowb = group_k;
      अगर (group_transb & 1) group_nrowb = group_n;

      info=-1;
	
      अगर (group_ldc < group_m) info = 13;
      अगर (group_ldb < group_nrowb)  info = 10;
      अगर (group_lda < group_nrowa)  info =  8;
      अगर (group_k < 0)        info =  5;
      अगर (group_n < 0)        info =  4;
      अगर (group_m < 0)        info =  3;
      अगर (group_transb < 0)        info =  2;
      अगर (group_transa < 0)        info =  1;      
#अगर defined(CBLAS)
    }
#पूर्ण_अगर

    अगर (info >= 0) {
      BLASFUNC(xerbla)(ERROR_NAME, &info, माप(ERROR_NAME));
      मुक्त(args_array);
      वापस;
    }

    /* Every member shares M and N: an empty member means an empty batch. */
    अगर (group_m == 0 || group_n == 0) { मुक्त(args_array); वापस; }

    group_mode=mode;

#अगर defined(SMP) || defined(SMALL_MATRIX_OPT)
    MNK = (द्विगुन) group_m * (द्विगुन) group_n * (द्विगुन) group_k;
#पूर्ण_अगर

#अगर_घोषित SMALL_MATRIX_OPT
#अगर !defined(COMPLEX)
    अगर (*(FLOAT *)(&alpha) != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 *(FLOAT *)(&alpha), *(FLOAT *)(&beta))){
#अन्यथा
    अगर (((FLOAT *)(alpha))[0] != ZERO &&
        GEMM_SMALL_MATRIX_PERMIT(group_transa, group_transb, group_m, group_n, group_k,
                                 ((FLOAT *)(alpha))[0], ((FLOAT *)(alpha))[1],
                                 ((FLOAT *)(beta))[0], ((FLOAT *)(beta))[1])){
#पूर्ण_अगर
      group_routine=शून्य;
#अगर !defined(COMPLEX)
      अगर(beta == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(gemm_small_kernel, ((group_transb<<2)|group_transa));
      }
#अन्यथा
      अगर(beta[0] == 0.0 && beta[1] == 0.0){
	group_mode=mode | BLAS_SMALL_B0_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel_b0, ((group_transb<<2)|group_transa));
      }अन्यथा{
	group_mode=mode | BLAS_SMALL_OPT;
	group_small_matrix_opt_routine=SMALL_KERNEL_ADDR(zgemm_small_kernel, ((group_transb<<2)|group_transa));
      }

#पूर्ण_अगर
      
    }अन्यथा{
#पूर्ण_अगर
      group_routine=(व्योम*)(gemm[(group_transb<<2)|group_transa]);
#अगर_घोषित SMALL_MATRIX_OPT
    }
#पूर्ण_अगर    

    
      args_array[i].m=group_m;
      args_array[i].n=group_n;
      args_array[i].k=group_k;
      args_array[i].lda=group_lda;
      args_array[i].ldb=group_ldb;
      args_array[i].ldc=group_ldc;
#अगर_घोषित COMPLEX
      args_array[i].alpha=alpha;
      args_array[i].beta=beta;
#अन्यथा
      args_array[i].alpha=&alpha;
      args_array[i].beta=&beta;
#पूर्ण_अगर

#अगर defined(CBLAS)
      अगर (order == CblasColMajor) {      
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
      }अन्यथा अगर(order == CblasRowMajor){
	args_array[i].a=&(b[i*strideb]);
	args_array[i].b=&(a[i*stridea]);
      }
#अन्यथा
	args_array[i].a=&(a[i*stridea]);
	args_array[i].b=&(b[i*strideb]);
#पूर्ण_अगर
      
      args_array[i].c= &c[i*stridec];
      
      args_array[i].routine_mode=group_mode;
      args_array[i].routine=group_routine;
#अगर_घोषित SMALL_MATRIX_OPT
      अगर (!group_routine)
        args_array[i].routine=group_small_matrix_opt_routine;
#पूर्ण_अगर      
  }

  अगर(गिनती>0) {
    GEMM_BATCH_THREAD(args_array,गिनती);
  }

  मुक्त(args_array);
}
