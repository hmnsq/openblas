#यदि_परिभाषितम् CBLAS

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव <string.h>
#अन्तर्भाव <stdarg.h>
#अन्तर्भाव "common.h"

OPENBLAS_EXPORT

शून्य CNAME(blasint p, वर्ण *rout, वर्ण *form, ...)
{
   va_list तर्काः;

   va_start(तर्काः, form);

   यदि (p)
      fprintf(stderr, "Parameter %lld to routine %s was incorrect\n", (दीर्घ दीर्घ)p, rout);
   vfprintf(stderr, form, तर्काः);
   va_end(तर्काः);
   निर्गच्छ(-1);
}
#यदि_अन्तः

