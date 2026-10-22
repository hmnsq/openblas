#جے_تعریف_شدہ CBLAS

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <string.h>
#ضمن <stdarg.h>
#ضمن "common.h"

OPENBLAS_EXPORT

فراغ CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list دلیلاں;

   va_start(دلیلاں, form);

   إذا (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (لما لما)p, rout);
   vfprintf(stderr, form, دلیلاں);
   va_end(دلیلاں);
   باہر_نکلو(-1);
}
#جے_انت

