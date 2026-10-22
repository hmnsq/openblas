#ifdef CBLAS

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <string.h>
#ضمن <stdarg.h>
#ضمن "common.h"

OPENBLAS_EXPORT

فراغ CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list حجوجي;

   va_start(حجوجي, form);

   إذا (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, حجوجي);
   va_end(حجوجي);
   exit(-1);
}
#endif

