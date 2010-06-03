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
#include <series.hpp>

#pragma hdrstop

#include "iface.h"
#include "lisrunf.h"
#include "lismain.h"
#include "lishelp.h"
#include "mprolog.h"

#pragma package(smart_init)


//------------------------------------------------------------------------------

#include "lisheadwin.h"
// writeTimeseries etc definitions
// get map names

//------------------------------------------------------------------------------

void OldMain( void )
{
}

void __fastcall LisThread::MakeSwitchList()
{
    char S[256], buf[256];
    FILE *fin = fopen(temprunname, "r");
    char IS = '=';
    int i, nrmaps = 0, j=0;
    bool go= false;

    // MAIN CHOICES OF lisem TYPE
    SwitchWheelAsChannel = LisIFace->LisemType == LISEMWHEELTRACKS;
    SwitchMulticlass = LisIFace->LisemType == LISEMMULTICLASS;
    SwitchNutrients = LisIFace->LisemType == LISEMNUTRIENTS;
    SwitchGullies = LisIFace->LisemType == LISEMGULLIES;

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

    do
    {
      i = fscanf(fin,"%[^\n]\n", S);
      if (strchr(S,IS))
      {
          char *p1 = strtok(S,"=");
          char *p = strtok(NULL,"=");

          if (strcmp(p1, "Outlet 1 file")==0) SwitchOutlet1 = p[0] != '\0';
          if (strcmp(p1, "Outlet 2 file")==0) SwitchOutlet2 = p[0] != '\0';
          if (strcmp(p1, "Outlet 3 file")==0) SwitchOutlet3 = p[0] != '\0';

          if (strcmp(p1, "No Erosion simulation")==0)          SwitchNoErosion =        p == "1";
          if (strcmp(p1, "Include main channels")==0)          SwitchIncludeChannel =   p == "1";
          if (strcmp(p1, "Include channel infil")==0)          SwitchChannelInfil =     p == "1";
          if (strcmp(p1, "Include channel baseflow")==0)       SwitchChannelBaseflow =  p == "1";
          if (strcmp(p1, "Include snowmelt")==0)               SwitchSnowmelt =         p == "1";
          if (strcmp(p1, "Alternative flow detachment")==0)    SwitchAltErosion =       p == "1";
          if (strcmp(p1, "Alternative depression storage")==0) SwitchAltDepression =    p == "1";
          if (strcmp(p1, "Include buffers")==0)                SwitchBuffers =          p == "1";
          if (strcmp(p1, "Include wheeltracks")==0)            SwitchInfilCompact =     p == "1";
          if (strcmp(p1, "Include grass strips")==0)           SwitchInfilGrass =       p == "1";
          if (strcmp(p1, "Include crusts")==0)                 SwitchInfilCrust =       p == "1";
          if (strcmp(p1, "Impermeable sublayer")==0)           SwitchImpermeable =      p == "1";
          if (strcmp(p1, "Matric head files")==0)              SwitchDumphead =         p == "1";
          if (strcmp(p1, "Geometric mean Ksat")==0)            SwitchGeometricMean =    p == "1";
          if (strcmp(p1, "Runoff maps in l/s/m")==0)           SwitchRunoffPerM =       p == "1";
          if (strcmp(p1, "Timeseries as PCRaster")==0)         SwitchWritePCRnames =    p == "1";
          if (strcmp(p1, "Timeplot as PCRaster")==0)           SwitchWritePCRtimeplot = p == "1";
          if (strcmp(p1, "Regular runoff output")==0)          SwitchOutputTimeStep =   p == "1";
          if (strcmp(p1, "User defined output")==0)            SwitchOutputTimeUser =   p == "1";
          if (strcmp(p1, "No erosion at outlet")==0)           SwitchNoErosionOutlet =  p == "1";
          if (strcmp(p1, "Subsoil drainage")==0)               SwitchDrainage =         p == "1";
          if (strcmp(p1, "Gully infiltration")==0)             SwitchGullyInit =        p == "1";
          if (strcmp(p1, "Use initial gully dimensions")==0)   SwitchGullyInfil =       p == "1";
          if (strcmp(p1, "CheckOutputMaps")==0)
          {
              char *q = strtok(p,",");SwitchMapoutRunoff= q == "1";
              q = strtok(NULL,",");   SwitchMapoutConc  = q == "1";
              q = strtok(NULL,",");   SwitchMapoutWH    = q == "1";
              q = strtok(NULL,",");   SwitchMapoutRWH   = q == "1";
              q = strtok(NULL,",");   SwitchMapoutTC    = q == "1";
              q = strtok(NULL,",");   SwitchMapoutEros  = q == "1";
              q = strtok(NULL,",");   SwitchMapoutDepo  = q == "1";
              q = strtok(NULL,",");   SwitchMapoutV     = q == "1";
              q = strtok(NULL,",");   SwitchMapoutInf   = q == "1";
              q = strtok(NULL,",");   SwitchMapoutSs    = q == "1";
              q = strtok(NULL,",");   SwitchMapoutChvol = q == "1";
          }

          if (strcmp(S, "CheckOutputMapsMC")==0)
          {
                char *q = strtok(p,","); SwitchMapoutMC0 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC1 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC2 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC3 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC4 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC5 = q == "1";
                q = strtok(NULL,","); SwitchMapoutMC6 = q == "1";
          }
          if (strcmp(S, "CheckOutputMapsNUT")==0)
          {
                char *q = strtok(p,","); SwitchMapoutPsol =   q == "1";
                q = strtok(NULL,",");   SwitchMapoutPsus =   q == "1";
                q = strtok(NULL,",");   SwitchMapoutPinf =   q == "1";
                q = strtok(NULL,",");   SwitchMapoutNH4sol = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNH4sus = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNH4inf = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNO3sol = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNO3sus = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNO3inf = q == "1";
                q = strtok(NULL,",");   SwitchMapoutPdep =   q == "1";
                q = strtok(NULL,",");   SwitchMapoutNH4dep = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNO3dep = q == "1";
                q = strtok(NULL,",");   SwitchMapoutPdet =   q == "1";
                q = strtok(NULL,",");   SwitchMapoutNH4det = q == "1";
                q = strtok(NULL,",");   SwitchMapoutNO3det = q == "1";
          }
          if (strcmp(S, "CheckOutputMapsGUL")==0)
          {
                char *q = strtok(p,",");SwitchMapoutGul0 = q == "1";
                q = strtok(NULL,","); SwitchMapoutGul1 = q == "1";
                q = strtok(NULL,","); SwitchMapoutGul2 = q == "1";
                q = strtok(NULL,","); SwitchMapoutGul3 = q == "1";
                q = strtok(NULL,","); SwitchMapoutGul4 = q == "1";
          }
        }
    }while (i >= 0);
}

//------------------------------------------------------------------------------

void __fastcall LisThread::MakeMapnameList()
{
    char S[256], buf[256];
    FILE *fin = fopen(temprunname, "r");
    char IS = '=';
    int i, nrmaps = 0, j=0;
    bool go= false;

    do
    {
      i = fscanf(fin,"%[^\n]\n", S);
      if (strstr(S,"[map names]"))
         go = true;
      if (go && strchr(S,IS))
         nrmaps++;
    }while (i >= 0);
    // count nr of maps

    if (namelist)
    {
       free(namelist);
       namelist = NULL;
    }
    fclose(fin);

    nrnamelist = nrmaps;
    namelist = (nameList *)malloc(nrnamelist*sizeof(nameList));
    go = false;
    fin = fopen(temprunname, "r");
    do
    {
      i = fscanf(fin,"%[^\n]\n", S);
      if (strstr(S,"[map names]"))
         go = true;
      if (go && strchr(S,IS))
      {
         char *p;
         p = strtok(S,"=");
         strcpy(namelist[j].ID, strupr(p));
         p = strtok(NULL,"=");
         strcpy(namelist[j].fullname, strlwr(p));
         j++;
      }
    }while (i >= 0);

   fclose(fin);
}
//---------------------------------------------------------------------------
char* __fastcall LisThread::Lmapname(char *vname)
{
   vname = strupr(vname);
   for (int i = 0; i < nrnamelist; i++)
     if (strcmp(vname, namelist[i].ID) == 0)
     {
       AnsiString mapname = ExtractFileName(namelist[i].fullname);
       if (mapname.IsEmpty())
       {
         char *p;
         p = (char *) malloc(128);
         sprintf(p,"map variable <%s> has no filename attached to it",vname);
         return(p);
         //return("NO MAPNAME DEFINED");
       }
       else
         return(namelist[i].fullname);
      }
   return "wrong internal map ID";
//VJ 040220 added return error to check for coding error of filename
}
//---------------------------------------------------------------------------

// start of main program
void __fastcall LisThread::Execute()
//begin_model{
//"begin_model{" and "}end_model" keep track of how deep you are for the map structures
//defined in csf.h
{

     /************************
      all booleans here are options and selections from the interface, given
      to a boolean in the main loop. This is to separate the interface from
      the model as much as possible
      ***********************/
/*
      bool SwitchCorrectMass = true;
      //VJ 080217:
      //NOTE: normal kinematic OF and kinematic channel have no mas balance correction
      //only used now in gully mass balance and wheeltrack water mass balance

      bool SwitchCorrectMassSED = true;
      // correct mass balanc of sediment after kin wave for all verions:
      // normal, gully, wheeltrack, nutrients, multiclass sediment etc.

      //applyDiagonalInKinematic = false;//(BOOL)LisIFace->CheckDiagonalDX->Checked;
      ///OBSOLETE
        // never use diagonal in kin wave, causes huge mass balance errors

      bool SwatreInitialized = false;
      bool SwitchCrustPresent = false;
      bool SwitchGrassPresent = false;
      bool SwitchWheelPresent = false;
      bool SwitchCompactPresent = false;
      bool SwitchInfilGA2 = false;
        // infiltration booleans

// ****** Switches to try different options, link with interface
      bool SwitchIncludeChannel = LisIFace->CheckIncludeChannel->Checked;
   //VJ 080214 include baseflow
      bool SwitchChannelBaseflow = LisIFace->CheckChannelBaseflow->Checked;
      bool startbaseflowincrease = false; // needed to start increase of baseflo after OF start into channel
      bool SwitchChannelInfil = LisIFace->CheckChannelInfil->Checked;
//VJ 090225 throw all water and sed in channel in outletcell, in ftaction to channel
      bool SwitchAllinChannel = true; //LisIFace->CheckAllinChannel->Checked;

      bool SwitchKinwaveInfil = true;//LisIFace->CheckKinwaveInfil->Checked;
      bool SwitchNoErosion   =  LisIFace->CheckNoErosion->Checked;
      bool SwitchAltErosion   =  LisIFace->CheckAltErosion->Checked;
      bool SwitchAltDepression   =  LisIFace->CheckAltDepression->Checked;

//VJ 070909 add macropore network
// NOT IMPLEMENTED
//      bool SwitchMacroporeFlow   =  LisIFace->CheckMacroporeFlow->Checked;

// VJ 040514 include buffers
      bool SwitchBuffers   =  LisIFace->CheckBuffers->Checked;
// VJ 080423 Snowmelt
      bool SwitchSnowmelt   =  LisIFace->CheckSnowmelt->Checked;

      bool SwitchRunoffPerM = LisIFace->CheckRunoffPerM->Checked;
      bool SwitchInfilCompact = LisIFace->CheckInfilCompact->Checked;
      bool SwitchInfilCrust = LisIFace->CheckInfilCrust->Checked;
      bool SwitchInfilGrass = LisIFace->CheckInfilGrass->Checked;

      bool SwitchImpermeable =  LisIFace->CheckImpermeable->Checked;
      bool SwitchDumphead = LisIFace->CheckDumphead->Checked;
      bool SwitchGeometricMean = LisIFace->CheckGeometric->Checked;

      bool SwitchWheelAsChannel = LisIFace->CheckWheelAsChannel;
      bool SwitchMulticlass = LisIFace->CheckMulticlass;
      bool SwitchNutrients = LisIFace->CheckNutrients;
      bool SwitchGullies = LisIFace->CheckGullies;

      bool SwitchMapoutRunoff = LisIFace->CheckMapout0->Checked;
      bool SwitchMapoutConc = LisIFace->CheckMapout1->Checked;
      bool SwitchMapoutWH = LisIFace->CheckMapout2->Checked;
      bool SwitchMapoutRWH = LisIFace->CheckMapout3->Checked;
      bool SwitchMapoutTC = LisIFace->CheckMapout4->Checked;
      bool SwitchMapoutEros = LisIFace->CheckMapout5->Checked;
      bool SwitchMapoutDepo = LisIFace->CheckMapout6->Checked;
      bool SwitchMapoutV = LisIFace->CheckMapout7->Checked;
      bool SwitchMapoutInf = LisIFace->CheckMapout8->Checked;
      bool SwitchMapoutSs = LisIFace->CheckMapout9->Checked;
      bool SwitchMapoutChvol = LisIFace->CheckMapout10->Checked;

      bool SwitchMapoutMC0 = LisIFace->CheckMapoutMC0->Checked;
      bool SwitchMapoutMC1 = LisIFace->CheckMapoutMC1->Checked;
      bool SwitchMapoutMC2 = LisIFace->CheckMapoutMC2->Checked;
      bool SwitchMapoutMC3 = LisIFace->CheckMapoutMC3->Checked;
      bool SwitchMapoutMC4 = LisIFace->CheckMapoutMC4->Checked;
      bool SwitchMapoutMC5 = LisIFace->CheckMapoutMC5->Checked;
      bool SwitchMapoutMC6 = LisIFace->CheckMapoutMC6->Checked;

      bool SwitchMapoutPsol = LisIFace->CheckMapoutNut0->Checked;
      bool SwitchMapoutPsus = LisIFace->CheckMapoutNut1->Checked;
      bool SwitchMapoutPinf = LisIFace->CheckMapoutNut2->Checked;
      bool SwitchMapoutNH4sol = LisIFace->CheckMapoutNut3->Checked;
      bool SwitchMapoutNH4sus = LisIFace->CheckMapoutNut4->Checked;
      bool SwitchMapoutNH4inf = LisIFace->CheckMapoutNut5->Checked;
      bool SwitchMapoutNO3sol = LisIFace->CheckMapoutNut6->Checked;
      bool SwitchMapoutNO3sus = LisIFace->CheckMapoutNut7->Checked;
      bool SwitchMapoutNO3inf = LisIFace->CheckMapoutNut8->Checked;
      bool SwitchMapoutPdep = LisIFace->CheckMapoutNut9->Checked;
      bool SwitchMapoutNH4dep = LisIFace->CheckMapoutNut10->Checked;
      bool SwitchMapoutNO3dep = LisIFace->CheckMapoutNut11->Checked;
      bool SwitchMapoutPdet = LisIFace->CheckMapoutNut12->Checked;
      bool SwitchMapoutNH4det = LisIFace->CheckMapoutNut13->Checked;
      bool SwitchMapoutNO3det = LisIFace->CheckMapoutNut14->Checked;

      bool SwitchMapoutGul0 = LisIFace->CheckMapoutGul0->Checked;
      bool SwitchMapoutGul1 = LisIFace->CheckMapoutGul1->Checked;
      bool SwitchMapoutGul2 = LisIFace->CheckMapoutGul2->Checked;
      bool SwitchMapoutGul3 = LisIFace->CheckMapoutGul3->Checked;
      bool SwitchMapoutGul4 = LisIFace->CheckMapoutGul4->Checked;

      bool SwitchWritePCRnames = LisIFace->CheckWritePCRnames->Checked;
      bool SwitchWritePCRtimeplot = LisIFace->CheckWritePCRtimeplot->Checked;
//VJ 080621 better check for output
      bool SwitchOutlet1 = !LisIFace->E_OutletName->Text.IsEmpty();
      bool SwitchOutlet2 = !LisIFace->E_Outlet1Name->Text.IsEmpty();
      bool SwitchOutlet3 = !LisIFace->E_Outlet2Name->Text.IsEmpty();

      bool SwitchNoErosionOutlet = LisIFace->CheckNoErosionOutlet->Checked;
//VJ 040224 added if outlet no detachment
      bool SwitchGullyEqualWD = false;
//VJ 040329 added equal erosion over width and depth in gully, set in paramgully
      bool SwitchGullyInit = LisIFace->CheckGullyInit->Checked;
//VJ 060404 added gully infiltration
      bool SwitchGullyInfil = LisIFace->CheckGullyInfil->Checked;
//VJ 040331 Included init gully dimensions
      bool SwitchDrainage = LisIFace->CheckSubsoilDrainage->Checked;
//VJ 050812 Included drainage from subsoil in G&A
      char ErosionUnits = LisIFace->E_OutputUnits->ItemIndex;
//VJ 050822 three units to choose from, kg/cell, kg/m2, ton/ha
      bool SwitchPestout = LisIFace->CheckPestout;

     bool SwitchOutputTimeStep = LisIFace->E_OutputTimeStep->Checked;
     bool SwitchOutputTimeUser = LisIFace->E_OutputTimeUser->Checked;
*/
      /*********** end bool options ****************************/

     ksatCalibration = (double)LisIFace->CSpinEditKsat->Value/100.0;

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

     void *soilModelcrackV = NULL;
     void *soilModelwheelV = NULL;
     void *soilModelcrustV = NULL;
     void *soilModelcompactV = NULL;
     void *soilModelgrassV = NULL;
     // declaration pointers needed here for cleanup at end of this file, not very elegant

     double *mu_cl;
     mu_cl=(double *)malloc(sizeof(double)*6);
     // texture classes (size in mu), declaration needed here for cleanup

      //void *SortLDD = NULL;
      // currently not used!!

      FILE *fpoutflow1;
      FILE *fpoutflow2;
      FILE *fpoutflow3;
      FILE *BufferFout;
      // hydrograph output files declaration

      int INFIL_METHOD; //infiltration method number, see lisparaminfil.cpp

      bool PERIODMAPS = false;
      size_t PERIOD;
      size_t PERIODCOUNT = 1;

//VJ 050913 commandline param pestout
		float pestreport = 0;
      int pestcounter = 0;
      FILE *fpestout;

      float timestepindex = 0;

      MakeMapnameList();
      // make an input/output mapname list from the temp runfile written by the interface

      MakeSwitchList();
      // get all switches (booleans) from the temp runfile written by the interface



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
// ******** Infiltration initialisation ****************************************
// -----------------------------------------------------------------------------

//#include "lisparaminfil.cpp"

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
     size_t timeIndex = 0; // for output of maps 

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

       #include "lisPgridcell.inc"

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
//        if (LisIFace->done || Terminated)
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

    free(mu_cl);
      // free memory used for texture classes

    RunDone = true; // tread is finished

  } // end loop try..__finally
}
//}end_model
// "}end_model" is a define in csf.h



