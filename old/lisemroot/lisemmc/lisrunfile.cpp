//---------------------------------------------------------------------------
#include <vcl.h>
#include <dir.h>
#pragma hdrstop

#include "lisrunfile.h"
//#include "lisinit.h"
#include "iface.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
//---------------------------------------------------------------------------
void __fastcall TLisIFace::MakeNewRunfile(AnsiString name)
{
    TStringList *IniOp = new TStringList;

    strcpy(workdir, E_Workdir->Text.c_str());

    IniOp->Clear();
    IniOp->Add("[LISEM for WINDOWS run file v3]");
    IniOp->Add("");
    IniOp->Add("[LISEM main type]");
    IniOp->Add("LISEM Type="+(AnsiString)LisemType);
    IniOp->Add("");
    IniOp->Add("[Work Directory]");
    IniOp->Add("WorkDir="+(AnsiString)workdir);
    IniOp->Add("");
    IniOp->Add("[Input]");
    IniOp->Add("Map Directory="+E_MapDir->Text);
    IniOp->Add("Table Directory="+E_TableDir->Text);
    IniOp->Add("Rainfall Directory="+RainfallDir);
    IniOp->Add("Rainfall file="+E_RainfallName->Text);
    IniOp->Add("Incude Snowmelt="+(AnsiString)(short)CheckSnowmelt->Checked);
    IniOp->Add("Snowmelt Directory="+SnowmeltDir);
    IniOp->Add("Snowmelt file="+E_SnowmeltName->Text);
    IniOp->Add("");
    IniOp->Add("[Output]");
    IniOp->Add("Result Directory="+E_ResultDir->Text);
    IniOp->Add("Main results file="+E_TotalName->Text);
    IniOp->Add("Filename point output="+E_OutletName->Text);
    IniOp->Add("Report point output separate="+(AnsiString)(short)CheckSeparateOutput->Checked);
    IniOp->Add("Report point output for SOBEK="+(AnsiString)(short)CheckSOBEKOutput->Checked);
    IniOp->Add("SOBEK date string="+SOBEKDateString->Text);
//    IniOp->Add("Total Runoff map="+E_TotalRunoffName->Text);
    IniOp->Add("Erosion map="+E_ErosionName->Text);
    IniOp->Add("Deposition map="+E_DepositionName->Text);
    IniOp->Add("Soilloss map="+E_SoillossName->Text);
    IniOp->Add("");
    IniOp->Add("[Simulation times]");
    IniOp->Add("Begin time="+E_begintime->Text);
    IniOp->Add("End time="+E_Endtime->Text);
    IniOp->Add("Timestep="+E_Timestep->Text);
    IniOp->Add("");
    IniOp->Add("[General options]");
    IniOp->Add("No Erosion simulation="+(AnsiString)(short)CheckNoErosion->Checked);
    IniOp->Add("Include main channels="+(AnsiString)(short)CheckIncludeChannel->Checked);
//VJ 050704 include channel infiltration
    IniOp->Add("Include channel infil="+(AnsiString)(short)CheckChannelInfil->Checked);
//VJ 080217 include channel baseflow
    IniOp->Add("Include channel baseflow="+(AnsiString)(short)CheckChannelBaseflow->Checked);
//VJ 090225 include outlet cell
    IniOp->Add("All water and sediment to outlet="+(AnsiString)(short)CheckAllinChannel->Checked);
//VJ 080614 include snowmelt
    IniOp->Add("Include snowmelt="+(AnsiString)(short)CheckSnowmelt->Checked);
//VJ 040224 added if outlet no erosion or deposition
    IniOp->Add("No erosion at outlet="+(AnsiString)(short)CheckNoErosionOutlet->Checked);
//VJ 050301 alternative erosion
    IniOp->Add("Alternative flow detachment="+(AnsiString)(short)CheckAltErosion->Checked);
//VJ 090214 alternative depr storage
    IniOp->Add("Simple depression storage="+(AnsiString)(short)CheckSimpleDepression->Checked);
//VJ 100202 include hard surface switch
    IniOp->Add("Hard Surfaces="+(AnsiString)(short)CheckHardsurface->Checked);

    IniOp->Add("");
    IniOp->Add("[Interception]");
    //VJ 100116 Use canopy storage map
    IniOp->Add("Use canopy storage map="+(AnsiString)(short)!CheckInterceptionLAI->Checked);
   //VJ 100116 8 equations to choose
    IniOp->Add("Canopy storage equation="+(AnsiString)(short)E_InterceptionLAIType->ItemIndex);
    IniOp->Add("Stemflow fraction="+(AnsiString)(double)FractionStemflow->Value);

//VJ 040514 added if buffers yes or no
    IniOp->Add("");
    IniOp->Add("[Conservation]");
    IniOp->Add("Include grass strips="+(AnsiString)(short)CheckInfilGrass->Checked);
    IniOp->Add("Grassstrip Mannings n="+(AnsiString)(double)CalibrateGrassN->Value);
    IniOp->Add("Include buffers="+(AnsiString)(short)CheckBuffers->Checked);
    IniOp->Add("Sediment bulk density="+E_SedBulkDensity->Text);
    IniOp->Add("Include Sediment traps="+(AnsiString)(short)CheckSedtrap->Checked);

    IniOp->Add("");
    IniOp->Add("[Calibration]");
    IniOp->Add("Ksat calibration="+(AnsiString)(double)CalibrateKsat->Value);
    IniOp->Add("N calibration="+(AnsiString)(double)CalibrateN->Value);
    IniOp->Add("Channel Ksat calibration="+(AnsiString)(double)CalibrateChKsat->Value);
    IniOp->Add("Channel N calibration="+(AnsiString)(double)CalibrateChN->Value);
    IniOp->Add("Splash Delivery Ratio="+(AnsiString)(double)CalibrateSplashDelivery->Value);
    //+E_SplashDelivery->Text);

    IniOp->Add("");
    IniOp->Add("[Gully options]");
    IniOp->Add("Fcrit relation="+(AnsiString)(int)E_Fcritical->ItemIndex);
    IniOp->Add("Threshold gradient="+E_ThresholdGrad->Text);
    IniOp->Add("QW relation="+(AnsiString)(int)E_QWrelation->ItemIndex);
    IniOp->Add("QW param A="+E_QWparama->Text);
    IniOp->Add("QW param B="+E_QWparamb->Text);
    IniOp->Add("Gully infiltration="+(AnsiString)(short)CheckGullyInfil->Checked);
    IniOp->Add("Use initial gully dimensions="+(AnsiString)(short)CheckGullyInit->Checked);

    IniOp->Add("");
    IniOp->Add("[Infiltration]");
    IniOp->Add("Infil Method="+(AnsiString)E_InfilMethoda->ItemIndex);
    IniOp->Add("Include wheeltracks="+(AnsiString)(short)CheckInfilCompact->Checked);
    IniOp->Add("Include crusts="+(AnsiString)(short)CheckInfilCrust->Checked);
    IniOp->Add("Impermeable sublayer="+(AnsiString)(short)CheckImpermeable->Checked);
    IniOp->Add("Subsoil drainage="+(AnsiString)(short)CheckSubsoilDrainage->Checked);
    IniOp->Add("SWATRE internal minimum timestep="+E_SwatreDTSEC->Text);
    IniOp->Add("Matric head files="+(AnsiString)(short)CheckDumphead->Checked);
    IniOp->Add("Geometric mean Ksat="+(AnsiString)(short)CheckGeometric->Checked);

    IniOp->Add("");
    IniOp->Add("[Output maps]");
    IniOp->Add("Runoff maps in l/s/m="+(AnsiString)(short)CheckRunoffPerM->Checked);
    IniOp->Add("Timeseries as PCRaster="+(AnsiString)(short)CheckWritePCRnames->Checked);
    IniOp->Add("Timeplot as PCRaster="+(AnsiString)(short)CheckWritePCRtimeplot->Checked);
    IniOp->Add("Regular runoff output="+(AnsiString)(short)E_OutputTimeStep->Checked);
//VJ 050822 three units to choose from, kg/cell, kg/m2, ton/ha
    IniOp->Add("Erosion map units (0/1/2)="+(AnsiString)(short)E_OutputUnits->ItemIndex);
    if (E_OutputTimeStep->Checked)
       IniOp->Add("Output interval="+(AnsiString)E_OutputTimeSteps->Value);
    else
       IniOp->Add("Output interval=0");
    IniOp->Add("User defined output="+(AnsiString)(short)E_OutputTimeUser->Checked);
    if (!E_OutputTimeUser->Checked)
       IniOp->Add("Output times=");
    else
       IniOp->Add("Output times="+E_OutputTimes->Lines->CommaText);

    // output maps BASIC
    AnsiString OutCheck = "CheckOutputMaps=";
    int vi[20]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,11) == "CheckMapout" &&
           Components[i]->Name.Length() < 14)
    {
       TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
       int j = atoi(Components[i]->Name.c_str()+11);
       if (cb->Checked) vi[j]=1;
          else vi[j]=0;
    }
    //VJ 040218 changed 7 into 20, error
    for (int i = 0; i < 20; i++)
    if (vi[i]>-1)
      OutCheck = OutCheck+vi[i]+",";
    OutCheck.Delete(OutCheck.Length(),1);
    IniOp->Add(OutCheck);

        // output maps NUTRIENTS
    OutCheck = "CheckOutputMapsNUT=";
    for (int i = 0; i < 20; i++) vi[i]=-1;
//VJ corrected error here: vi[20]=-1;
    for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,14) == "CheckMapoutNut")
    {
       TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
       int j = atoi(Components[i]->Name.c_str()+14);
       if (cb->Checked)vi[j]=1;
          else vi[j]=0;
    }
    for (int i = 0; i < 20; i++)
    if (vi[i]>-1)
      OutCheck = OutCheck+vi[i]+",";
    OutCheck.Delete(OutCheck.Length(),1);
    IniOp->Add(OutCheck);

        // output maps MULTICLASS
    OutCheck = "CheckOutputMapsMC=";
    for (int i = 0; i < 20; i++) vi[i]=-1;
//VJ corrected error here: vi[20]=-1;
    for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,13) == "CheckMapoutMC")
    {
       TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
       int j = atoi(Components[i]->Name.c_str()+13);
       if (cb->Checked)vi[j]=1;
          else vi[j]=0;
    }
    for (int i = 0; i < 20; i++)
    if (vi[i]>-1)
      OutCheck = OutCheck+vi[i]+",";
    OutCheck.Delete(OutCheck.Length(),1);
    IniOp->Add(OutCheck);
    // output maps GULLIES
    OutCheck = "CheckOutputMapsGUL=";
    for (int i = 0; i < 20; i++) vi[i]=-1;
//VJ corrected error here: vi[20]=-1;
    for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,14) == "CheckMapoutGul")
    {
       TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
       int j = atoi(Components[i]->Name.c_str()+14);
       if (cb->Checked)vi[j]=1;
          else vi[j]=0;
    }
    for (int i = 0; i < 20; i++)
    if (vi[i]>-1)
      OutCheck = OutCheck+vi[i]+",";
    OutCheck.Delete(OutCheck.Length(),1);
    IniOp->Add(OutCheck);
    IniOp->Add("");
    IniOp->Add("[Texture classes]");
    AnsiString S = "ClassMu=";
    for (int i = 1; i<7; i++)
    {
       S = S +TextureClass->Cells[0][i];
       if ( i < 6) S = S +",";
    }
    IniOp->Add(S);
//    IniOp->Add("TCCal="+E_TCCalibration->Text);

    IniOp->Add("");
    IniOp->Add("[map names]");

    for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,4) == "Maps")
    {
         TStringGrid* S =(TStringGrid *)LisIFace->Components[i];
         IniOp->Add("");
         IniOp->Add("["+S->Name.SubString(5,128)+"]");
         for (int j = 1; j < S->RowCount; j++)
             if (!S->Cells[0][j].IsEmpty())
             IniOp->Add(S->Cells[5][j]+"="+S->Cells[1][j]);
    }
    /*
    S->Cells[0] = varname
    S->Cells[1] = filename
    S->Cells[2] = desc
    S->Cells[3] = directory
    S->Cells[4] = full path
    S->Cells[5] = ID up case
    */
    IniOp->SaveToFile(name);

    delete IniOp;
}
//---------------------------------------------------------------------------
char *ParseValue(TStringList *ts, char *text)
{
   int index = ts->IndexOfName(text);
   if (index < 0)
      return(NULL);
   char *p = strchr(ts->Strings[index].c_str(),'=');
   p++;
   return(p);
}
//---------------------------------------------------------------------------
bool __fastcall TLisIFace::ReadNewRunfile(AnsiString name)
{
    TStringList *IniOp = new TStringList;
    try
    {
        AnsiString _workdir, tmp;

        if (FileExists(name))
           IniOp->LoadFromFile(name);
        else
        {
           delete IniOp;
           return false;
        }

        bool runv3 = IniOp->Strings[0] == "[LISEM for WINDOWS run file v3]";

        _workdir = ExtractFilePath(name);

//        _workdir = IniOp->Values[(AnsiString)"WorkDir"];
        _workdir = CheckDir("Work Directory",_workdir);
        if (DirectoryExists(_workdir))
        {
           chdir(_workdir.c_str());
           strcpy(workdir, _workdir.c_str());
           E_Workdir->Text = _workdir;
        }

        CheckSeparateOutput->Checked = IniOp->Values[(AnsiString)"Report point output separate"] == "1";
        CheckSOBEKOutput->Checked = IniOp->Values[(AnsiString)"Report point output for SOBEK"] == "1";
        SOBEKDateString->Text = IniOp->Values["SOBEK date string"];

        CheckWritePCRnames->Checked = IniOp->Values[(AnsiString)"Timeseries as PCRaster"] == "1";
        CheckWritePCRtimeplot->Checked = IniOp->Values[(AnsiString)"Timeplot as PCRaster"] == "1";
        CheckRunoffPerM->Checked = IniOp->Values[(AnsiString)"Runoff maps in l/s/m"] == "1";

        CheckNoErosion->Checked = IniOp->Values[(AnsiString)"No Erosion simulation"] == "1";
        CheckIncludeChannel->Checked = IniOp->Values[(AnsiString)"Include main channels"] == "1";
//VJ 050704 include channel infiltration
        CheckChannelInfil->Checked = IniOp->Values[(AnsiString)"Include channel infil"] == "1";
//VJ 080217 include channel baseflow
        CheckChannelBaseflow->Checked = IniOp->Values[(AnsiString)"Include channel baseflow"] == "1";
//VJ 090225 include outlet cell
        CheckAllinChannel->Checked = IniOp->Values[(AnsiString)"All water and sediment to outlet"] == "1";
//VJ 080614 include snowmelt
        CheckSnowmelt->Checked = IniOp->Values[(AnsiString)"Include snowmelt"] == "1";

//VJ 030702 changed channel infil to kin wave infil
//        CheckKinwaveInfil->Checked = IniOp->Values[(AnsiString)"Kinematic wave infil"] == "1";
//VJ 040224 added if outlet no erosion or deposition
        CheckNoErosionOutlet->Checked = IniOp->Values[(AnsiString)"No erosion at outlet"] == "1";

        CheckInfilGrass->Checked = IniOp->Values[(AnsiString)"Include grass strips"] == "1";
        tmp = ParseValue(IniOp, "Grassstrip Mannings n");
        if (!tmp.IsEmpty()) CalibrateGrassN->Value = tmp.ToDouble();
//VJ 040214 added if buffers
        CheckBuffers->Checked = IniOp->Values[(AnsiString)"Include buffers"] == "1";
        CheckSedtrap->Checked = IniOp->Values[(AnsiString)"Include Sediment traps"] == "1";
        if (runv3)
         E_SedBulkDensity->Text = IniOp->Values[(AnsiString)"Sediment bulk density"];
        else
         E_SedBulkDensity->Text = "1800";

//VJ 050301 alternative erosion
	     CheckAltErosion->Checked = IniOp->Values[(AnsiString)"Alternative flow detachment"] == "1";
//VJ 090214 alternative depr storage
	     CheckSimpleDepression->Checked = IniOp->Values[(AnsiString)"Simple depression storage"] == "1";
        // to read old run files:
	     CheckSimpleDepression->Checked = IniOp->Values[(AnsiString)"Alternative depression storage"] == "1";
//VJ 100202 include hard surface switch
	     CheckHardsurface->Checked = IniOp->Values[(AnsiString)"Hard Surfaces"] == "1";

//VJ 100116 Use canopy storage map
        if (runv3)
           CheckInterceptionLAI->Checked = IniOp->Values[(AnsiString)"Use canopy storage map"] == "0";
        else
           CheckInterceptionLAI->Checked = true;
//VJ 100116 8 equations to choose
        if (runv3)
          E_InterceptionLAIType->ItemIndex = IniOp->Values[(AnsiString)"Canopy storage equation"].ToInt();
        else
          E_InterceptionLAIType->ItemIndex = 0;

        E_begintime->Text = IniOp->Values[(AnsiString)"Begin time"];
        E_Endtime->Text = IniOp->Values[(AnsiString)"End time"];
        E_Timestep->Text = IniOp->Values[(AnsiString)"Timestep"];

        if (IniOp->Values[(AnsiString)"Infil Method"].IsEmpty())
        {
           E_InfilMethoda->ItemIndex = IniOp->Values[(AnsiString)"Method"].ToInt();
        }
        else
           E_InfilMethoda->ItemIndex = IniOp->Values[(AnsiString)"Infil Method"].ToInt();

        E_SwatreDTSEC->Text = ParseValue(IniOp, "SWATRE internal minimum timestep");
        CheckDumphead->Checked = IniOp->Values[(AnsiString)"Matric head files"] == "1";
        CheckGeometric->Checked = IniOp->Values[(AnsiString)"Geometric mean Ksat"] == "1";

        ResetInfil();

        CheckInfilCompact->Checked = IniOp->Values[(AnsiString)"Include wheeltracks"] == "1";
        CheckInfilCrust->Checked = IniOp->Values[(AnsiString)"Include crusts"] == "1";
        CheckImpermeable->Checked = IniOp->Values[(AnsiString)"Impermeable sublayer"] == "1";
		  CheckSubsoilDrainage->Checked = IniOp->Values[(AnsiString)"Subsoil drainage="] == "1";

        CalibrateKsat->Value = IniOp->Values[(AnsiString)"Ksat calibration"].ToDouble();
    //    if (!runv3) CalibrateKsat->Value /= 100;
        tmp = IniOp->Values[(AnsiString)"N calibration"];
        if (!tmp.IsEmpty()) CalibrateN->Value = tmp.ToDouble();
        tmp = IniOp->Values[(AnsiString)"Channel Ksat calibration"];
        if (!tmp.IsEmpty()) CalibrateChKsat->Value = tmp.ToDouble();
        tmp = IniOp->Values[(AnsiString)"Channel N calibration"];
        if (!tmp.IsEmpty()) CalibrateChN->Value = tmp.ToDouble();
        tmp = IniOp->Values[(AnsiString)"Splash Delivery Ratio"];
        if (!tmp.IsEmpty()) CalibrateSplashDelivery->Value = tmp.ToDouble();

        if (CalibrateKsat->Value != 1.0 || CalibrateN->Value != 1.0 ||
           CalibrateChKsat->Value != 1.0 || CalibrateChN->Value != 1.0)
             JvLookOut1->ActivePage = LookOutPage4;

/******* DIR CHECKING **********/

//VJ 030625 changed
        OldMapDir = ParseValue(IniOp,"Map Directory");
        E_MapDir->Text = CheckDir("Maps Directory",OldMapDir);
        if (E_MapDir->Text == "DIR_NOT_EXIST")
           E_MapDir->Text = "";
        // do with oldmapdir to prevent endless checking of changes in E_MapDir
        E_TableDir->Text = ParseValue(IniOp,"Table Directory");
        E_TableDir->Text = CheckDir("Table Directory", E_TableDir->Text);
        if (E_TableDir->Text == "DIR_NOT_EXIST")
           E_TableDir->Text = "";

        RainfallDir = ParseValue(IniOp,"Rainfall Directory");
        RainfallDir = CheckDir("Rainfall Directory",RainfallDir);
        if (RainfallDir == "DIR_NOT_EXIST")
        {
           RainfallDir = "";
           E_RainfallName->Text = "";
        }
        else
        {
            E_RainfallName->Text = ParseValue(IniOp,"Rainfall file");
            E_RainfallName->Text = ExtractFileName(E_RainfallName->Text);
        }
        strcpy(RainfallNamePath, RainfallDir.c_str());
        strcat(RainfallNamePath,E_RainfallName->Text.c_str());
 //VJ 080614 include snowmelt
        CheckSnowmelt->Checked = IniOp->Values[(AnsiString)"Incude Snowmelt="] == "1";
        if (CheckSnowmelt->Checked)
        {
          SnowmeltDir = ParseValue(IniOp,"Snowmelt Directory");
          SnowmeltDir = CheckDir("Snowmelt Directory",SnowmeltDir);
          if (SnowmeltDir == "DIR_NOT_EXIST")
          {
             SnowmeltDir = "";
             E_SnowmeltName->Text = "";
             CheckSnowmelt->Checked =false;
          }
          else
          {
              E_SnowmeltName->Text = ParseValue(IniOp,"Snowmelt file");
              E_SnowmeltName->Text = ExtractFileName(E_SnowmeltName->Text);
          }
          strcpy(SnowmeltNamePath, SnowmeltDir.c_str());
          strcat(SnowmeltNamePath,E_SnowmeltName->Text.c_str());
        }


        E_ResultDir->Text = ExpandFileName(ParseValue(IniOp ,"Result Directory"));
        if (!DirectoryExists(E_ResultDir->Text))
        {
//VJ 031105 2.159 do not ask otherwise batch operations will hold on this question
// forcedirectory cannot handle emty strings, throws exception

          if(!E_ResultDir->Text.IsEmpty())
          {
            if (!batchrun)
            {
               AnsiString S = "Result Directory: \"" + E_ResultDir->Text + "\" does not exist, create it?";
               if (Application->MessageBox(S.c_str(), "LISEM Warning", MB_OKCANCEL+MB_ICONERROR) == IDOK)
                ForceDirectories(E_ResultDir->Text);
               else
                E_ResultDir->Text = "";
            }
            else
            if (!ForceDirectories(E_ResultDir->Text))
            {
        			CheckError(AnsiString("Result cannot be created: ")+E_ResultDir->Text);
            }
          }
        }
        // no further checking

        E_TotalName->Text = ParseValue(IniOp ,"Main results file");
        E_OutletName->Text = ParseValue(IniOp ,"Outlet 1 file");// leave old filename for old run files 
        //E_Outlet1Name->Text = ParseValue(IniOp ,"Outlet 2 file");
        //E_Outlet2Name->Text = ParseValue(IniOp ,"Outlet 3 file");

        E_OutletName->Text = ParseValue(IniOp ,"Filename point output"); // if newer fould do that
//        E_TotalRunoffName->Text = ParseValue(IniOp ,"Total Runoff map");
//VJ 100115 total runoff map
        E_ErosionName->Text = ParseValue(IniOp ,"Erosion map");
        E_DepositionName->Text = ParseValue(IniOp ,"Deposition map");
        E_SoillossName->Text = ParseValue(IniOp ,"Soilloss map");

//VJ 050822 three units to choose from, kg/cell, kg/m2, ton/ha
	     E_OutputUnits->ItemIndex = IniOp->Values[(AnsiString)"Erosion map units (0/1/2)"].ToIntDef(2);

        E_OutputTimeStep->Checked = IniOp->Values[(AnsiString)"Regular runoff output"] == "1";
        if (E_OutputTimeStep->Checked)
           E_OutputTimeSteps->Value = IniOp->Values[(AnsiString)"Output interval"].ToInt();

        char *p, *line;
        E_OutputTimeUser->Checked = IniOp->Values[(AnsiString)"User defined output"] == "1";
        if (E_OutputTimeUser->Checked)
        {
            line = strchr(IniOp->Strings[IniOp->IndexOfName("Output times")].c_str(),'=');
            line++;
            p = strtok(line, ",");
            E_OutputTimes->Clear();
            if (p)
               E_OutputTimes->Lines->Add(p);

            while (p)
            {
                p = strtok(NULL, ",");
                if (p)
                   E_OutputTimes->Lines->Add(p);
            }

       }

       //**** GULLY PARAMS ***********
       E_ThresholdGrad->Text = ParseValue(IniOp,"Threshold gradient");
       E_QWparama->Text = ParseValue(IniOp,"QW param A");
       E_QWparamb->Text = ParseValue(IniOp,"QW param B");
//VJ 050907 changed this to tointdef
       E_Fcritical->ItemIndex = IniOp->Values[(AnsiString)"Fcrit relation"].ToIntDef(0);
       E_QWrelation->ItemIndex = IniOp->Values[(AnsiString)"QW relation"].ToIntDef(0);
       CheckGullyInfil->Checked = IniOp->Values[(AnsiString)"Gully infiltration"] == "1";
       CheckGullyInit->Checked = IniOp->Values[(AnsiString)"Use initial gully dimensions"] == "1";

       //***** MULTICLASS SEDIMENT *****
       if (IniOp->IndexOfName("ClassMu") != -1)
       {
           line = strchr(IniOp->Strings[IniOp->IndexOfName("ClassMu")].c_str(),'=');
           int i = 1;
           line++;

           p = strtok(line, ",");
           if (p)
              TextureClass->Cells[0][1] = p;
           while (p)
           {
              p = strtok(NULL, ",");
              if (p)
                 TextureClass->Cells[0][++i] = p;
           }
       }


       // ********* MAP NAMES ***********
    /*
    S->Cells[0] = varname
    S->Cells[1] = filename
    S->Cells[2] = desc
    S->Cells[3] = directory
    S->Cells[4] = full path
    S->Cells[5] = ID up case
    */

       for (int i = 0; i < ComponentCount; i++)
       if (Components[i]->Name.SubString(0,4) == "Maps")
       {
         TStringGrid* S =(TStringGrid *)LisIFace->Components[i];
         for (int j = 1; j < S->RowCount; j++)
         {
           S->Cells[4][j] = IniOp->Values[S->Cells[5][j]];
           S->Cells[1][j] = ExtractFileName(S->Cells[4][j]);

//VJ 030606 added expandfilename
          S->Cells[3][j] = ExpandFileName(ExtractFileDir(S->Cells[4][j]));
           if (S->Cells[3][j].IsEmpty())
              S->Cells[3][j] = E_MapDir->Text;
         }
       }

       AnsiString CheckOut = ParseValue(IniOp ,"CheckOutputMaps");
       AnsiString CheckOutMC = ParseValue(IniOp ,"CheckOutputMapsMC");
       AnsiString CheckOutNUT = ParseValue(IniOp ,"CheckOutputMapsNUT");
       AnsiString CheckOutGUL = ParseValue(IniOp ,"CheckOutputMapsGUL");
       for (int i = 0; i < ComponentCount; i++)
       {
         if (Components[i]->Name.SubString(0,11) == "CheckMapout")
         {
            if(Components[i]->Name.Length() < 14)
            {
              TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
              int j = atoi(Components[i]->Name.c_str()+11);
              cb->Checked = atoi(CheckOut.c_str()+2*j) == 1;
            }
            if (Components[i]->Name.SubString(0,14) == "CheckMapoutGul")
            {
              TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
              int j = atoi(Components[i]->Name.c_str()+14);
              cb->Checked = atoi(CheckOutGUL.c_str()+2*j) == 1;
            }
            if (Components[i]->Name.SubString(0,13) == "CheckMapoutMC")
            {
              TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
              int j = atoi(Components[i]->Name.c_str()+13);
              cb->Checked = atoi(CheckOutMC.c_str()+2*j) == 1;
            }
            if (Components[i]->Name.SubString(0,14) == "CheckMapoutNut")
            {
              TCheckBox *cb = (TCheckBox *)LisIFace->Components[i];
              int j = atoi(Components[i]->Name.c_str()+14);
              cb->Checked = atoi(CheckOutNUT.c_str()+2*j) == 1;
            }

         }
       }


       if (CheckMulticlass)
            TextureClassCheck();

       InitOutMapNames();

    }
    catch(...)
    {
    	delete IniOp;
        return false;
    }

    delete IniOp;
    return true;
}



