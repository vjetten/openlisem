//---------------------------------------------------------------------------


#pragma hdrstop

#include "lismapminmax.h"
#include "csf.h"
#include "iface.h"
#include "cps.h"

 //---------------------------------------------------------------------------

#pragma package(smart_init)

//sets map series min and max when stopping run
void __fastcall TLisIFace::SetTimeseriesMinmax(char *filename)
{
    int i = 0;
    MAP *m;
    REAL8 minsv=+1e30, maxsv=-1e30;

    for (i = 1; i <= LastPCRTimestep; i ++)
    {
        REAL8 minv, maxv;
        char *fn = MakePCRFileName(filename,i);
        m = Mopen(fn, M_READ);
        if (m == NULL)
        {
           return;
        }

        RgetMaxVal(m, &maxv);
        RgetMinVal(m, &minv);
        if (maxsv < maxv) maxsv = maxv;
        if (minsv > minv) minsv = minv;
        Mclose(m);
    }

    for (i = 1; i <= LastPCRTimestep; i ++)
    {
        char *fn = MakePCRFileName(filename,i);
        m = Mopen(fn, M_WRITE);
        if (m == NULL)
        {

           return;
        }
        RputMaxVal(m, &maxsv);
        RputMinVal(m, &minsv);
        Mclose(m);
    }




}
