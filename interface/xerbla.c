#अगर_घोषित CBLAS

#समावेश <stdio.h>
#समावेश <stdlib.h>
#समावेश <string.h>
#समावेश <stdarg.h>
#समावेश "common.h"

OPENBLAS_EXPORT

व्योम CNAME(blasint p, अक्षर *rout, अक्षर *form, ...)
{
   बहु_सूची args;

   बहु_शुरू(args, form);

   अगर (p)
      ख_लिखो(मानक_त्रुटि, "Parameter %lld to routine %s was incorrect\n", (दीर्घ दीर्घ)p, rout);
   भख_लिखो(मानक_त्रुटि, form, args);
   बहु_पूर्ण(args);
   निकास(-1);
}
#पूर्ण_अगर

