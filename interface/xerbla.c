#ifdef CBLAS

#अन्तर्भाव <stdio.h>
#अन्तर्भाव <stdlib.h>
#अन्तर्भाव <string.h>
#अन्तर्भाव <stdarg.h>
#अन्तर्भाव "common.h"

OPENBLAS_EXPORT

शून्य CNAME(blasint p, अक्षर *rout, अक्षर *form, ...)
{
   बहु_सूची args;

   बहु_शुरू(args, form);

   जर (p)
      ख_लिखो(मानक_त्रुटि, "Parameter %lld to routine %s was incorrect\n", (दीर्घ दीर्घ)p, rout);
   भख_लिखो(मानक_त्रुटि, form, args);
   बहु_पूर्ण(args);
   निकास(-1);
}
#endif

