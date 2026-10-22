#اگر_تعریف_شدہ CBLAS

#شامل <stdio.h>
#شامل <stdlib.h>
#شامل <string.h>
#شامل <stdarg.h>
#شامل "common.h"

OPENBLAS_EXPORT

خالی CNAME(blasint p, حرف *rout, حرف *form, ...)
{
   va_list دلائل;

   va_start(دلائل, form);

   اگر (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (طویل طویل)p, rout);
   vfprintf(stderr, form, دلائل);
   va_end(دلائل);
   باہر_نکلو(-1);
}
#اگر_ختم

