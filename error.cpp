#include "error.h"


QString ErrorString;    // declare here, referenced by error.h

void Error(QString s)
{
    ErrorString = s;
    throw 1;
}
