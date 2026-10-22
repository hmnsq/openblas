#اگر_تعریف_شده CBLAS

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <string.h>
#ضمن <stdarg.h>
#ضمن "common.h"

OPENBLAS_EXPORT

فراغ CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list آرگومانها;

   va_start(آرگومانها, form);

   إذا (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (دراز دراز)p, rout);
   vfprintf(stderr, form, آرگومانها);
   va_end(آرگومانها);
   خروج(-1);
}
#پایان_اگر

