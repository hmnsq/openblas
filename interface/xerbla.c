#ifdef CBLAS

#ಅನ್ತರ್ಭಾವ <stdio.h>
#ಅನ್ತರ್ಭಾವ <stdlib.h>
#ಅನ್ತರ್ಭಾವ <string.h>
#ಅನ್ತರ್ಭಾವ <stdarg.h>
#ಅನ್ತರ್ಭಾವ "common.h"

OPENBLAS_EXPORT

ಶೂನ್ಯ CNAME(blasint p, ವರ್ಣ *rout, ವರ್ಣ *form, ...)
{
   va_list args;

   va_start(args, form);

   ಯದಿ (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

