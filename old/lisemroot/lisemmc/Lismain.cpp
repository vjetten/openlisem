/******************************************************************************
Last changes: 10-12-2002, Victor

DOS COMPATIBILITY REMOVED

NOTE: no more static spatial vars, too much trouble cleaning them up when
      thread closes

NOTE: often error mV on LDD in kinematic wave caused by Qin that gets the
channels mask, so MVs where no channel !!!

Comments / changes are indicated as
//ADR
//VJ
//CW
//+++OBSOLETE+++

A NOTE ON FILE NAMES:

     something with sediment
Channelxxx      a channel variable
xxxDX           part cellsize of landuse in m (e.g. RoadWidthDX)
Totalxxx        is used for cumulative total for whole run, nonspatial
Sumxxx          is used for total of a map in a single timestep, nonspatial
SumTotxxx       is used for Total in run of spatial total of each timestep, spatial
Totxxx          is a temp total calculated on the way, spatial
xxxin           a variable before the kinematic wave
xxxout          a variable after the kinematic wave

Note:
1) variables that are declared with _spatial; _nonspatial; _spatialinpuy; _nonspatialinput
are put on a stack and can be declared anywhere in the code, and are visible anywhere in the code
when using "calc", because calc lookes on that stack for the variable. This circumvents the
global versu locAL declaration of variables. That is why single value parameters are also
declared with _nonspatial, instead of simply "float foo;" (non maps)

2) A large part of the code is written as calc statements. These are strnigs that are
interpreted (parsed) by the calc module. The strings are not checked at compile time
becausethe interpretation is done when the model runs. A writing error in calc(" something")
will only show up when you start to run the model.

3) basically the model is one big loop. To retain an understadable code parts of the code are
put in separate files, which are copied in the main text when the model compiles. These are
not functions, just simply pieces of code in seperate files. This makes the variable declarartion
simple, but it is not very flexible or elegant.
******************************************************************************/

#include <vcl.h>
#include <string.h>
#include <series.hpp>

#pragma hdrstop

#include "iface.h"

#pragma package(smart_init)

//------------------------------------------------------------------------------
extern "C" void myHandler(const char *msg)
{
// myHandler is linked to Errorhandler and MerrorHandler that is called from calc, cps, csf
// and replaces stderr (dos only) so when something goes wrong an exeption
// is thrown in the try...catch loop in lismain
    throw Exception(msg);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::LisemWarn(AnsiString s, AnsiString s1)
{
   WarningText = s+s1;
   Synchronize(LisemWarningV);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::SetTimeseriesMinmax(char *filename)
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
//---------------------------------------------------------------------------
// LET OP, ONLY SWITCHES, VALUES ARE READ BY GETFLOAT AND GETINT
void __fastcall LisThread::ParseInputData()
{
    char S[256], buf[256];
    int i, nrmaps = 0, j=0;
    bool go= false;


    SwitchCorrectMass = true;
    //VJ 080217:
    //NOTE: normal kinematic OF and kinematic channel have no mas balance correction
    //only used now in gully mass balance and wheeltrack water mass balance
    SwitchCorrectMassSED = true;
    // correct mass balanc of sediment after kin wave for all verions:
    // normal, gully, wheeltrack, nutrients, multiclass sediment etc.
    //applyDiagonalInKinematic = false;//()LisIFace->CheckDiagonalDX->Checked;
    ///OBSOLETE
      // never use diagonal in kin wave, causes huge mass balance errors
     SwatreInitialized = false;
     SwitchCrustPresent = false;
     SwitchGrassPresent = false;
     SwitchWheelPresent = false;
     SwitchCompactPresent = false;
     SwitchInfilGA2 = false;
    // initialize infil switches
     startbaseflowincrease = false;
     SwitchAllinChannel = true;
     SwitchGullyEqualWD = false;

//comments on added switches ported from the older code
//VJ 080214 include baseflow
//VJ 090225 throw all water and sed in channel in outletcell, in ftaction to channel
//VJ 040514 include buffers
//VJ 080423 Snowmelt
//VJ 040224 added if outlet no detachment
//VJ 040329 added equal erosion over width and depth in gully, set in paramgully
//VJ 060404 added gully infiltration
//VJ 040331 Included init gully dimensions
//VJ 050812 Included drainage from subsoil in G&A
//VJ 050822 three units to choose from, kg/cell, kg/m2, ton/ha
    buf[0] = '\0';
    j = 0;
    outflowFileName[0] = '\0';
   // outflowFileName2[0] = '\0';
   // outflowFileName3[0] = '\0';
   // outPointFileName[0] = '\0';

    for (j = 0; j < nrnamelist; j++)//do
    {
          char p1[128],p[128];
          int iii = namelist[j].iii;
          float vvv = namelist[j].vvv;
          strcpy(p1, namelist[j].name);
          strcpy(p, namelist[j].value);

          if (strcmp(p1, "LISEM Type")==0 && p)
          {
              SwitchWheelAsChannel = iii == LISEMWHEELTRACKS;
              SwitchMulticlass = iii == LISEMMULTICLASS;
              SwitchNutrients = iii == LISEMNUTRIENTS;
              SwitchGullies = iii == LISEMGULLIES;
          }
          if (strcmp(p1, "Map Directory")==0 && p) strcpy(PATH, p);
          if (strcmp(p1, "Result Directory")==0 && p) strcpy(RESPATH, p);
          if (strcmp(p1, "Table Directory")==0 && p) strcpy(tableDir, p);
          if (strcmp(p1, "Main results file")==0 && p) strcpy(resultFileName, CatPath(p, RESPATH));
//          if (strcmp(p1, "Total Runoff map")==0 && p) strcpy(totalRunoffFileName, CatPath(p, RESPATH));
//VJ 100115 total runoff map
          if (strcmp(p1, "Erosion map")==0 && p) strcpy(totalErosionFileName, CatPath(p, RESPATH));
          if (strcmp(p1, "Deposition map")==0 && p) strcpy(totalDepositionFileName, CatPath(p, RESPATH));
          if (strcmp(p1, "Soilloss map")==0 && p) strcpy(totalSoillossFileName,CatPath(p, RESPATH));
          if (strcmp(p1, "Filename point output")==0 && p) strcpy(outflowFileName, CatPath(p, RESPATH));
          /*
          if (strcmp(p1, "Outlet 1 file")==0 && p){
            strcpy(outflowFileName, CatPath(p, RESPATH));
            SwitchOutlet1 = p[0] != '\0';
          }
          if (strcmp(p1, "Outlet 2 file")==0 && p) {
            strcpy(outflowFileName2, CatPath(p, RESPATH));
            SwitchOutlet2 = p[0] != '\0';
          }
          if (strcmp(p1, "Outlet 3 file")==0 && p) {
            strcpy(outflowFileName3, CatPath(p, RESPATH));
            SwitchOutlet3 = p[0] != '\0';
          }
          */
          if (strcmp(p1, "Rainfall Directory")==0 && p) strcpy(buf, p);
          if (strcmp(p1, "Rainfall file")==0 && p) strcpy(rainFileName, CatPath(p, buf));
          if (strcmp(p1, "Snowmelt Directory")==0 && p) strcpy(buf, p);
          if (strcmp(p1, "Snowmelt file")==0 && p) strcpy(snowmeltFileName, CatPath(p, buf));

          memset(WRITETIME,0,100);
          if (strcmp(p1, "Output times")==0 && p)
          {
              int l = 0;
              char *q = strtok(p,","); if (q) { WRITETIME[l] = atof(q); l++; }
              while (q) {
                q = strtok(NULL,",");
                if (q) {
                   WRITETIME[l] = atof(q);
                   l++;
                 }
              }
          }

          if (strcmp(p1, "ClassMu")==0 && p)
          {
              char *q = strtok(p,","); if (q) { mu_cl[0] = atof(q); }
              q = strtok(NULL,","); if (q) { mu_cl[1] = atof(q); }
              q = strtok(NULL,","); if (q) { mu_cl[2] = atof(q); }
              q = strtok(NULL,","); if (q) { mu_cl[3] = atof(q); }
              q = strtok(NULL,","); if (q) { mu_cl[4] = atof(q); }
              q = strtok(NULL,","); if (q) { mu_cl[5] = atof(q); }
          }

          //options in the main code, order is not important
          if (strcmp(p1, "No Erosion simulation")==0)          SwitchNoErosion =        iii == 1;
          if (strcmp(p1, "Include main channels")==0)          SwitchIncludeChannel =   iii == 1;
          if (strcmp(p1, "Include channel infil")==0)          SwitchChannelInfil =     iii == 1;
          if (strcmp(p1, "Include channel baseflow")==0)       SwitchChannelBaseflow =  iii == 1;
          if (strcmp(p1, "Include snowmelt")==0)               SwitchSnowmelt =         iii == 1;
          if (strcmp(p1, "Alternative flow detachment")==0)    SwitchAltErosion =       iii == 1;
          if (strcmp(p1, "Simple depression storage")==0)      SwitchSimpleDepression = iii == 1;
          if (strcmp(p1, "Hard Surfaces")==0)                  SwitchHardsurface      = iii == 1;
          if (strcmp(p1, "Include buffers")==0)                SwitchBuffers =          iii == 1;
          if (strcmp(p1, "Include wheeltracks")==0)            SwitchInfilCompact =     iii == 1;
          if (strcmp(p1, "Include grass strips")==0)           SwitchInfilGrass =       iii == 1;
          if (strcmp(p1, "Include crusts")==0)                 SwitchInfilCrust =       iii == 1;
          if (strcmp(p1, "Impermeable sublayer")==0)           SwitchImpermeable =      iii == 1;
          if (strcmp(p1, "Matric head files")==0)              SwitchDumphead =         iii == 1;
          if (strcmp(p1, "Geometric mean Ksat")==0)            SwitchGeometricMean =    iii == 1;
          if (strcmp(p1, "Runoff maps in l/s/m")==0)           SwitchRunoffPerM =       iii == 1;
          if (strcmp(p1, "Timeseries as PCRaster")==0)         SwitchWritePCRnames =    iii == 1;
          if (strcmp(p1, "Timeplot as PCRaster")==0)           SwitchWritePCRtimeplot = iii == 1;
          if (strcmp(p1, "Regular runoff output")==0)          SwitchOutputTimeStep =   iii == 1;
          if (strcmp(p1, "User defined output")==0)            SwitchOutputTimeUser =   iii == 1;
          if (strcmp(p1, "No erosion at outlet")==0)           SwitchNoErosionOutlet =  iii == 1;
          if (strcmp(p1, "Subsoil drainage")==0)               SwitchDrainage =         iii == 1;
          if (strcmp(p1, "Gully infiltration")==0)             SwitchGullyInfil =       iii == 1;
          if (strcmp(p1, "Use initial gully dimensions")==0)   SwitchGullyInit =        iii == 1;
//VJ  091211
          if (strcmp(p1, "Report point output separate")==0)   SwitchSeparateOutput =   iii == 1;
          if (strcmp(p1, "Report point output for SOBEK")==0)   SwitchSOBEKOutput =   iii == 1;

          if (strcmp(p1, "SOBEK date string")==0)
          {
              char *q = strtok(p,",");
              strncpy(SOBEKdatestring, q, 10);

          }
//VJ 100116 Interception
          if (strcmp(p1, "Use canopy storage map")==0)   SwitchInterceptionLAI =        iii == 0;

          //outputmaps that are selected (1) or not (0)
          if (strcmp(p1, "CheckOutputMaps")==0)
          {
              char *q = strtok(p,",");SwitchMapoutRunoff= strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutConc  = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutWH    = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutWHC   = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutTC    = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutEros  = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutDepo  = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutV     = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutInf   = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutSs    = strcmp(q,"1") == 0;
              q = strtok(NULL,",");   SwitchMapoutChvol = strcmp(q,"1") == 0;
          }
          if (strcmp(S, "CheckOutputMapsMC")==0)
          {
                char *q = strtok(p,","); SwitchMapoutMC0 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC1 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC2 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC3 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC4 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC5 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutMC6 = strcmp(q,"1") == 0;
          }
          if (strcmp(S, "CheckOutputMapsNUT")==0)
          {
                char *q = strtok(p,","); SwitchMapoutPsol =   strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutPsus =   strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutPinf =   strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNH4sol = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNH4sus = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNH4inf = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNO3sol = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNO3sus = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNO3inf = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutPdep =   strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNH4dep = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNO3dep = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutPdet =   strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNH4det = strcmp(q,"1") == 0;
                q = strtok(NULL,",");   SwitchMapoutNO3det = strcmp(q,"1") == 0;
          }
          if (strcmp(S, "CheckOutputMapsGUL")==0)
          {
                char *q = strtok(p,",");SwitchMapoutGul0 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutGul1 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutGul2 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutGul3 = strcmp(q,"1") == 0;
                q = strtok(NULL,","); SwitchMapoutGul4 = strcmp(q,"1") == 0;
          }
     }
}

//------------------------------------------------------------------------------
void __fastcall LisThread::GetRunFile()
{
    char S[256], buf[256];
    FILE *fin = fopen(temprunname, "r");
    char IS = '=';
    int i, j, k;
    int lines = 0, vars = 0;
    char ch = '\0';

    do
    {
       ch = fgetc(fin);
       if (ch == '=')
          vars++;
       if (ch == '\n')
          lines++;
    }while (ch != EOF);
    lines++;
    vars++;
    fclose(fin);

    if (namelist)
    {
       free(namelist);
       namelist = NULL;
    }

    namelist = (nameList *)malloc(vars*sizeof(nameList));
    for (i = 0; i < vars; i++)
    {
        namelist[i].name[0] = '\0';
        namelist[i].value[0] = '\0';
        namelist[i].type = 0;
        namelist[i].vvv = 0;
        namelist[i].iii = 0;
        namelist[i].done = false;
    }
    nrnamelist = vars;
    
    fin = fopen(temprunname, "r");
    k = 0;
    for (i = 0; i < lines; i++)
    {
       memset(S,'\0',255);
       j = fscanf(fin,"%[^\n]\n", S);
       if (i == 0)
          runv3 = (strstr(S, "[LISEM for WINDOWS run file v3]") > 0);
       if (strchr(S, '=') > 0 && k < vars)
       {

          char *p = strtok(S,"=");
          if (p) strcpy(namelist[k].name, p);
          p = strtok(NULL,"=");
          if (p) strcpy(namelist[k].value, p);
          if (namelist[k].value[0] == '\0')
             namelist[k].type = 0; //empty
          else
          {
             int res = 0;
             namelist[k].type = 1; //string
             float _v;
             int _i;
             res = sscanf (namelist[k].value,"%f", &_v);   
             res = sscanf (namelist[k].value,"%d", &_i);
             namelist[k].vvv = _v;
             namelist[k].iii = _i;
             if (res > 0)
                namelist[k].type = 2; //number
          }
          k++;
       }
    }
    fclose(fin);
}
//---------------------------------------------------------------------------
char* __fastcall LisThread::Lmapname(char *vname, int nr)
{
   vname = strupr(vname);
   for (int i = 0; i < nrnamelist; i++)
     if (strcmp(vname, strupr(namelist[i].name)) == 0)
     {
        char *split = strrchr(namelist[i].value,DIR_PATH_DELIM_CHAR);
//        if (split == NULL)
        if (strlen(namelist[i].value) == 0)
        {
          char *p;
          p = (char *) malloc(128);
          sprintf(p,"map variable <%s> has no filename attached to it",vname);
          return(p);
         //return("NO MAPNAME DEFINED");
        }
        else
        {
          if (nr > 0 && !namelist[i].done && split == NULL)
          {
            strcpy(namelist[i].value, CatPath(namelist[i].value, PATH));
            namelist[i].done = true;
          }
          return(namelist[i].value);
/*         else
         if (nr == 1)
         {
            *split++;
            split[strlen(split)-4] = '\0';
            return(split);
         }
*/         
        }
      }
   return "wrong internal map ID";
//VJ 040220 added return error to check for coding error of filename
}
//---------------------------------------------------------------------------
float __fastcall LisThread::GetFloat(char *vname)
{
   vname = strupr(vname);
   for (int i = 0; i < nrnamelist; i++)
     if (strcmp(vname, strupr(namelist[i].name)) == 0)
        return (namelist[i].vvv);
}
//---------------------------------------------------------------------------
int __fastcall LisThread::GetInt(char *vname)
{
   vname = strupr(vname);
   for (int i = 0; i < nrnamelist; i++)
     if (strcmp(vname, strupr(namelist[i].name)) == 0)
        return (namelist[i].iii);
}
//---------------------------------------------------------------------------
// start of main program
void __fastcall LisThread::Execute()
{
     InitDone = false;
     // mark start of loop, is set to true after first timestep
     //used to initialize totals

     setmem(ErrorMessage,'\0',sizeof(ErrorMessage));
       // errormessage is used in the cps, calc and csf libraries
       // is used with the try...catch loops below
     exitOnError = 0;
     // avoid exiting with fatal error
     errorHandler = myHandler;
     // alternative for error.c in misc.lib
     MerrorHandler = myHandler;
     // alternative for mperror.c in csf.lib

     soilModelcrackV = NULL;
     soilModelwheelV = NULL;
     soilModelcrustV = NULL;
     soilModelcompactV = NULL;
     soilModelgrassV = NULL;
     // declaration pointers needed here for cleanup at end of this file, not very elegant

     PERIODMAPS = false;
     PERIODCOUNT = 1;

//VJ 050913 commandline param pestout
	  pestreport = 0;
     pestcounter = 0;
     SwitchPestout = LisIFace->CheckPestout;
     LastPCRTimestep = 0;
       // last timestep for setting min and max values of saved mapseries when run is stopped

     GetRunFile();
      // get variable list from temp runfile, ansi C

     ParseInputData();
      // parse variable list: get all switches (booleans)


//VJ 031218 added try..__finally loop to ensure cleanup of memory structures
// try loops are interrupted by errors, and error message is displayed
try  // first try for the __finally
{
  try  // second try loop to catch errors
  {
//------------------------------------------------------------------------------
//------ OLD MAIN CODE STARTS HERE----------------------------------------------
//------------------------------------------------------------------------------
begin_model{
//"begin_model{" and "}end_model" keep track of how deep you are for the map structures
//defined in csf.h

//******************************************************************************
//****** INPUT of runfile VARIABLES AND MAPS ***********************************
//******************************************************************************

//NOTE: some maps are still made even when that LISEM version is not selected,
//but in that case ONLY the pointer to the map is created with "spatial_input_if"
//to avoid errors. This only costs 4 bytes per map!

#include "lisparaminput.cpp"     //general vars needed in all lisem types
#include "lisparaminfil.cpp"     //infiltration variables
#include "lisparamchannel.cpp"   //channel vars
#include "lisparamwheel.cpp"     //wheeltrack related vars
#include "lisparammulti.cpp"     //multi grainsize version vars
#include "lisparamnutrients.cpp" //nutrient version vars
#include "lisparamgully.cpp"     //gully version vars


// -----------------------------------------------------------------------------
// ******** Initialisation global vars *****************************************
// -----------------------------------------------------------------------------

#include "lisparaminit.cpp"

// -----------------------------------------------------------------------------
// ******** open files for output **********************************************
// -----------------------------------------------------------------------------

#include "lisfileoutopen.cpp"

// -----------------------------------------------------------------------------
// ******** TOTALS FOR MASS BALANCE ETC ****************************************
// -----------------------------------------------------------------------------

#include "lisparamtotals.cpp"

// *****************************************************************************
// =============================================================================
// -------- START of timestep loop ---------------------------------------------
// =============================================================================
// *****************************************************************************

 //    LisIFace->Messages->Lines->Append("Starting Loop...");

     int stepNr = 0;    // to flag first step for initialization of some variables
     size_t timeIndex = 0; // for output of maps, used in writetimeseries

     pestreport = LisIFace->PestoutTimeinterval/DTMIN;

     for(timestepindex = STARTINTERVAL+DTMIN;        // start
                timestepindex <  ENDINTERVAL+DTMIN;                // length
                timestepindex += DTMIN)                      // next timestep
     begin{

     BOOL firsttime = (timestepindex <= STARTINTERVAL+1.5*DTMIN);
     //VJ 080628 timeseries input firsttime=true means open file
     //VJ 080804 changed to 1.5*dtmin to make sure it works

           //OBSOLETE
           //static bool lastStep = false;
           //PRECOND(! lastStep); // can never be true
           //lastStep = (index+DTMIN) >= ENDINTERVAL;

// *****************************************************************************
// ****** set up sub-gridcell land use *****************************************
// *****************************************************************************

       #include "lisPgridcell.cpp"

// *****************************************************************************
// ****** rainfall & interception **********************************************
// *****************************************************************************

       #include "lisPrainfall.cpp"
//dorainfall(index);
       // rainfall and interception

// *****************************************************************************
// ****** INFILTRATION *********************************************************
// *****************************************************************************

       #include "lisPinfiltration.cpp"

// *****************************************************************************
// ******* nutrient in solution balance ****************************************
// *****************************************************************************

       if (SwitchNutrients)
       {
         #include "lisPNutrients.cpp"
         // calculate nutrients in solution
       }

// *****************************************************************************
// ******* SURFACE storage in micro-depressions ********************************
// *****************************************************************************
//VJ 030702 wheeltracks before surface storage because result is a change of WH
        if (SwitchWheelPresent)
        {
            #include "lisPwheeltrack.cpp"
        }

//  #include "lisPSurfstor.cpp"  obsolete
//        #include "lisPSurfstor2.cpp"
        #include "lisPSurfstor3.cpp"

// *****************************************************************************
// ******** SPLASH detachment **************************************************
// *****************************************************************************

        #include "lisPsplash.cpp"

// *****************************************************************************
// ******** FLOW detachment **************************************************
// *****************************************************************************

        if(SwitchMulticlass)
        {
             #include "lisPMCNUTerosdepo.cpp"
             // erosion/deposition of multiclass sed. and suspended NUTs
        }
        else
        {
             #include "lisPerosdepo.cpp"
             // basic erosion deposition
             //#include "lisPerosdepobofu.cpp"
             // erosion deposition according to Rose Hairsine
        }

// *****************************************************************************
// ******** transport to the channel systems ***********************************
// *****************************************************************************

        #include "lisPfractionto.cpp"
        // fraction of all material that flows to channel, gullies, wheeltracks:
        // water, sediment, MC sediment, nutrients

// *****************************************************************************
// ******** channel flow calculations ******************************************
// *****************************************************************************

        if (SwitchIncludeChannel)
        {
            if (SwitchMulticlass)
            {
                #include "lisPMCchannel.cpp"
                if (SwitchNutrients)
                {
                    #include "lisPNUTchannel.cpp"
                    // order is important: NUTchannel uses vars from MCchannel
                }
            }
            else
            {
                #include "lisPchannel.cpp"
                // basic erosion/deposition channels
                //#include "lisPchannelbofu.cpp"
                // basic erosion deposition hairsine rose model
            }
        }
//VJ 030702 moved this to before surface storage        
/*
        if (SwitchWheelPresent)
        {
            #include "lisPwheeltrack.cpp"
        }
*/
// *****************************************************************************
// ******** Gully incision  ****************************************************
// *****************************************************************************

        if (SwitchGullies)
        {
            #include "lisPgully.cpp"
        }
        // NOTE: ORDER is important do this BEFORE overland flow:
        // gully outflow is combined with overlandflow as one kinematic wave


// *****************************************************************************
// ******** kinematic wave calculation for overland flow ***********************
// *****************************************************************************

        if (SwitchMulticlass)
        {
            #include "lisPMCkinematic.cpp"
            if (SwitchNutrients)
            {
                #include "lisPNUTkinematic.cpp"
            }
            // ORDER is important, vars are reused in NUT from MC
        }
        else
        {
            #include "lisPkinematic.cpp"
        }

// *****************************************************************************
// ******** water mass balance *************************************************
// *****************************************************************************

        #include "lisPmassbal.cpp"

// *****************************************************************************
// ******** sediment mass balance **********************************************
// *****************************************************************************

         if (!SwitchNoErosion)
         {
             if(!SwitchMulticlass)
             {
                #include "lisPmassbalsed.cpp"
             }
             else
             {
                #include "lisPMCmassbalsed.cpp"
                if (SwitchNutrients)
                {
                    #include "lisPNUTmassbalsed.cpp"
                }
             }
         }

// *****************************************************************************
// ******** OUTPUT to file and screen ******************************************
// *****************************************************************************

//VJ 031218 moved to before file and screen output to have the map series set to minmax
// even when interrupted                         
        if (RunDone || Terminated)
        {
           timestepindex = ENDINTERVAL+1;
 //          lastStep = true;
        }

        //********maps and outlet file update *************
        // all output to maps and tabels
        #include "lisfileoutput.cpp"

        //********SCREEN UPDATE *************
        // here model vars are linked to screen display vars
        #include "lisscreenwin.cpp"


        //-------------- INCREASE STEP -----------------------
        if (stepNr == 1)
           LisIFace->Messages->Lines->Clear();
        stepNr++;

     }end // end of the the timeloop

 }end_model
// *****************************************************************************
// =============================================================================
// -------- END of timestep loop -----------------------------------------------
// =============================================================================
// *****************************************************************************

// VJ 031218 ADDED try.._finally loop to ensure cleanup of memory structures

  } //try..catch when an error occurs the program jumps to here
  catch(Exception &E)
  {

  //VJ 030612 - included error message to text file
    FILE *ferror;
    AnsiString SDate = DateToStr(Date());
    AnsiString STime = TimeToStr(Time());
    AnsiString errorfname = "liserror.txt";
    Beep();
    E.Message = E.Message + "in function:" + AnsiString(ErrorMessage);

    ferror = fopen(errorfname.c_str(),"w");
    fprintf(ferror, "LISEM eror generated on %s - %s\n",SDate.c_str(), STime.c_str());
    fprintf(ferror, "%s\n",E.Message.c_str());
    fclose(ferror);

    // screen message
    Application->MessageBox(E.Message.c_str(), "LISEM ERROR",MB_OK+MB_ICONERROR);

    // free maps memory structure used so far before the exception occurred
    for (int i = 0; i < currDepth; i++)
       CpsLeaveBlock();

//continues below
  }
}  // second try..catch
__finally // do this even if an exception occurred:
{

    LisIFace->t_end = Time();
    runtimelength = (LisIFace->t_end-LisIFace->t_begin);
    LisIFace->runduration.sprintf("%8.3f",runtimelength*1440);
    LisIFace->Messages->Lines->Add("run time: "+(AnsiString)LisIFace->runduration + "  min");
    // calc runtime length

//VJ 031218 check dumpscreen to halt processes to finish dumpscreen
    // make a screendump
    LisIFace->CheckDumpScreen = true;
    for (;;)
    {
      Synchronize(DumpScreenSync);
      if (LisIFace->CheckDumpScreen)
         break;
    }

  // *****************************************************************************
  // ******** free mmory structures **********************************************
  // *****************************************************************************

    // clean up swatre memory structures
    if (SwatreInitialized)
    {
         FreeSwatreInfo();
         if (soilModelcrackV)
            CloseSwatre((SOIL_MODEL *)soilModelcrackV);
         if (soilModelwheelV)
            CloseSwatre((SOIL_MODEL *)soilModelwheelV);
         if (soilModelgrassV)
            CloseSwatre((SOIL_MODEL *)soilModelgrassV);
         if (soilModelcrustV)
            CloseSwatre((SOIL_MODEL *)soilModelcrustV);
         if (soilModelcompactV)
            CloseSwatre((SOIL_MODEL *)soilModelcompactV);
    }


    FreeTimeSerieInput((timeSerie[1]));

    FreeTimeSerieInput((timeSerie[2]));
    //VJ 080617 added = NULL to fix loop for missing map
    timeSerie[1] = NULL;
    timeSerie[2] = NULL;
     // clean up rainfall files

    RunDone = true; // tread is finished

  } // end loop try..__finally
}




