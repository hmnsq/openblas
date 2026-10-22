#که_تعریف_شوی CBLAS

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <string.h>
#ضمن <stdarg.h>
#ضمن "common.h"

OPENBLAS_EXPORT

فراغ CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list دلیلونه;

   va_start(دلیلونه, form);

   إذا (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (اوږد اوږد)p, rout);
   vfprintf(stderr, form, دلیلونه);
   va_end(دلیلونه);
   وتل(-1);
}
#د_که_پای

