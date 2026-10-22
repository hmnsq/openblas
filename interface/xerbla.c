#ifdef CBLAS

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव <string.h>
#अन्तर्भाव <stdarg.h>
#अन्तर्भाव "common.h"

OPENBLAS_EXPORT

शून्य CNAME(blasint p, वर्ण *rout, वर्ण *form, ...)
{
   va_list args;

   va_start(args, form);

   जइ (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (long long)p, rout);
   vfprintf(stderr, form, args);
   va_end(args);
   exit(-1);
}
#endif

