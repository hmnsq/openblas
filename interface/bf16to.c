#అన్తర్భావ <stdio.h>
#అన్తర్భావ "common.h"
#నిర్వచితమైతే FUNCTION_PROFILE
#అన్తర్భావ "functable.h"
#అయితే_ముగింపు

#ఒకవేళ defined(DOUBLE_PREC)
#నిర్వచించు FLOAT_TYPE ద్విగుణ
#లేదావేళ defined(SINGLE_PREC)
#నిర్వచించు FLOAT_TYPE దశాంశం
#లేకపోతే
#అయితే_ముగింపు

OPENBLAS_EXPORT
#నిర్వచితం_కాకపోతే CBLAS
శూన్య NAME(blasint *N, bfloat16 *లోపల, blasint *INC_లోపల, FLOAT_TYPE *బహిర్గతం, blasint *INC_బహిర్గతం){
  BLASLONG n    = *N;
  BLASLONG inc_లోపల = *INC_లోపల;
  BLASLONG inc_బహిర్గతం = *INC_బహిర్గతం;

  PRINT_DEBUG_NAME;

  ఒకవేళ (n <= 0) తిరిగి;

  IDEBUG_ప్రారంభం;
  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (inc_లోపల < 0)   లోపల -= (n - 1) * inc_లోపల;
  ఒకవేళ (inc_బహిర్గతం < 0) బహిర్గతం -= (n - 1) * inc_బహిర్గతం;

#ఒకవేళ defined(DOUBLE_PREC)
  D_BF16_గమ్యం_K(n, లోపల, inc_లోపల, బహిర్గతం, inc_బహిర్గతం);
#లేదావేళ defined(SINGLE_PREC)
  S_BF16_గమ్యం_K(n, లోపల, inc_లోపల, బహిర్గతం, inc_బహిర్గతం);
#లేకపోతే
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(1, 2 * n, 2 * n);
  IDEBUG_చివర;
}
#లేకపోతే
శూన్య CNAME(blasint n, bfloat16 * లోపల, blasint inc_లోపల, FLOAT_TYPE * బహిర్గతం, blasint inc_బహిర్గతం){
  PRINT_DEBUG_CNAME;

  ఒకవేళ (n <= 0) తిరిగి;

  IDEBUG_ప్రారంభం;
  FUNCTION_PROFILE_ప్రారంభం();

  ఒకవేళ (inc_లోపల < 0)   లోపల -= (n - 1) * inc_లోపల;
  ఒకవేళ (inc_బహిర్గతం < 0) బహిర్గతం -= (n - 1) * inc_బహిర్గతం;

#ఒకవేళ defined(DOUBLE_PREC)
  D_BF16_గమ్యం_K(n, లోపల, inc_లోపల, బహిర్గతం, inc_బహిర్గతం);
#లేదావేళ defined(SINGLE_PREC)
  S_BF16_గమ్యం_K(n, లోపల, inc_లోపల, బహిర్గతం, inc_బహిర్గతం);
#లేకపోతే
#అయితే_ముగింపు

  FUNCTION_PROFILE_చివర(1, 2 * n, 2 * n);
  IDEBUG_చివర;
}
#అయితే_ముగింపు
