 //---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "lisstart.h"

USEFORM("iface.cpp", LisIFace);
USEFORM("lisstart.cpp", StartForm);
USEFORM("lishelp.cpp", HelpForm);
USEFORM("lisabout.cpp", AboutBox);
USEFORM("lisrunf.cpp", RunForm);
USEFORM("lisrainf.cpp", RainForm);
USEFORM("lisdirv.cpp", DirView);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
    try
    {
       Application->Initialize();

      Application->CreateForm(__classid(TLisIFace), &LisIFace);
       Application->CreateForm(__classid(TRunForm), &RunForm);
       Application->CreateForm(__classid(THelpForm), &HelpForm);
       Application->CreateForm(__classid(TAboutBox), &AboutBox);
       Application->CreateForm(__classid(TRainForm), &RainForm);
       Application->CreateForm(__classid(TDirView), &DirView);
       Application->Run();
        
    }
    catch (Exception &exception)
    {
        Application->ShowException(&exception);
    }
    return 0;
}
//---------------------------------------------------------------------------
