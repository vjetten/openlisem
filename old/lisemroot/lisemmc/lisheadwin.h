//------COMPILE AS win95-------------------------------------------------

//========comma delimited output===========
#define commaoutputstr "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\n"\
                       "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\"\n"
#define commaoutputstrMC "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\
                         ,\"mu0\",\"mu1\",\"mu2\",\"mu3\",\"mu4\",\"mu5\"\n"\
                         "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\",\
                          \"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\n"
#define commaoutputstrNUT "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\
                          ,\"mu0\",\"mu1\",\"mu2\",\"mu3\",\"mu4\",\"mu5\"\
                          ,\"Psol\",\"NH4sol\",\"NO3sol\",\"Psus\",\"NH4sus\",\"NO3sus\"\n"\
                          "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\"\
                          ,\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\
                          ,\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\n"

#define commaoutputstrBuffer "\"t\",\"P\",\"Water Vol\",\"Sed Vol\"\n" \
                             "\"min\", \"mm/h\", \"m3\", \"kg\"\n"

//VJ 050822 changed rainfall precision to 4 digits
#define commaformat4      "%.6g,%.4g,%.6g,%.6g\n"
#define commaformat       "%.6g,%.4g,%.6g,%.6g,%.6g\n"
#define commaformatMC     "%.6g,%.4g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n"
#define commaformatNUT    "%.6g,%.4g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n"

//========space delimited output===========
#define ncommaoutputstr "t   P    Q   Qsed Conc \n"     \
                        "min mm/h l/s kg/s g/l \n"
#define ncommaoutputstrMC "t   P     Q     Qsed   Conc  mu0   mu1   mu2   mu3   mu4   mu5\n"     \
                          "min mm/h  l/s   kg/s   g/l   g/l   g/l   g/l   g/l   g/l   g/l\n"
#define ncommaoutputstrNUT "t   P    Q   Qsed Conc  mu0   mu1   mu2   mu3   mu4   mu5  Psol  NH4sol  NO3sol Psus  NH4sus  NO3sus\n"\
                           "min mm/h l/s kg/s g/l   g/l   g/l   g/l   g/l   g/l   g/l  g/l   g/l     g/l    g/l   g/l     g/l\n"

#define ncommaoutputstrBuffer "t    P     Water Vol  Sed Vol\n" \
                              "min  mm/h  m3         kg\n"

//VJ 050822 changed rainfall precision to 4 digits
#define ncommaformat4      "%.6g %.4g %.6g %.6g\n"
#define ncommaformat       "%.6g %.4g %.6g %.6g %.6g\n"
#define ncommaformatMC     "%.6g %.4g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n"
#define ncommaformatNUT    "%.6g %.4g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n"
//VJ 030702 corrected number of times %.6g

//========timeplot headers
#define timeplotoutputstr "LISEM output\n 5 \n"\
                       "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"
#define timeplotoutputstrMC "LISEM output\n 11 \n"\
                       "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"\
                       "mu0 (g/l)\nmu1 (g/l)\nmu2 (g/l)\nmu3 (g/l)\nmu4 (g/l)\nmu5 (g/l)\n"
#define timeplotoutputstrNUT "LISEM output\n 11 \n"\
                          "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"\
                          "Psol (g/l)\nNH4sol (g/l)\nNO3sol (g/l)\nPsus (g/l)\nNH4sus (g/l)\nNO3sus (g/l)\n"
#define timeplotoutputstrBuffer "LISEM output\n 4 \n"\
                          "t (min)\nP (mm/h)\nWater Vol(m3)\nSed Vol (kg)\n"


//#define WRITE_END_RESULT       (writeEachTime || lastStep)
// obsolete

#define ReportV(S)    ReportText = S;

#define LisemWarning    LisemWarn

#define writeList(v, s) if (PERIODMAPS && PERIODCOUNT == PERIOD){\
                        strcpy(FileName, CatPath(s, RESPATH));\
                        write(v, MakeFileName(FileName,timestepindex));}

#define writepath(v, s) strcpy(FileName, CatPath(s, RESPATH));write(v, FileName)

#define writeTimeseries(v, s) if (SwitchWritePCRnames) write(v, MakePCRFileName(s,timeIndex+1));\
                              else write(v, MakeFileName(s,timestepindex))

//#define mapname(s) LisIFace->Lmapname(s)
#define mapname(s) Lmapname(s)

bool swatreBottomClosed = false;//(bool)SwitchImpermeable;

//---------------------------------------------------------------------------
extern "C" void myHandler(const char *msg)
{
// myHandler is linked to Errorhandler and MerrorHandler that is called from calc, cps, csf
// and replaces stderr (dos only) so when something goes wrong an exeption
// is thrown in the try...catch loop in lismain
    throw Exception(msg);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::DoHelp()
{
     HelpForm->ShowModal();
}
//---------------------------------------------------------------------------
void __fastcall LisThread::ShowHelp()
{
       Synchronize(DoHelp);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::LisemWarn(AnsiString s, AnsiString s1)
{
   WarningText = s+s1;
   Synchronize(LisemWarningV);
}
//---------------------------------------------------------------------------
static bool AnswerIsYes(char *answer)
        {return tolower(answer[0]) == 'y';}
//---------------------------------------------------------------------------

//------COMPILE AS win95-------------------------------------------------


