#جيڪڏهن_وصف_ٿيل CBLAS

#ضمن <stdio.h>
#ضمن <stdlib.h>
#ضمن <string.h>
#ضمن <stdarg.h>
#ضمن "common.h"

OPENBLAS_EXPORT

فراغ CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list دليلون;

   va_start(دليلون, form);

   إذا (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (ڊگهو ڊگهو)p, rout);
   vfprintf(stderr, form, دليلون);
   va_end(دليلون);
   ٻاهر_نڪرو(-1);
}
#جيڪڏهن_پورو

