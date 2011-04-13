/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
 \file LisUIrunfile.cpp
 \brief functions to read and parse the runfile
 *
 * GetRunfile: read runfile and put all variables into the 'namelist' structure
 * ParseInputData: namelist structure is parsed to fill the interface
 *
 * updateModelData: update variables with changes in interface, called by save file
 *
 */

#include "lisemqt.h"
#include "model.h"
#include "global.h"

//---------------------------------------------------------------------------
//VJ 110107 corrections
// fill namelist with the actual runfile data but correct for old runfiles
// so that faulty data or obsolete vars are ignored
void lisemqt::GetRunfile()
{
   QFile fin(op.runfilename);
   if (!fin.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, "openLISEM",
                           QString("Cannot read file %1:\n%2.")
                           .arg(op.runfilename)
                           .arg(fin.errorString()));
      return;
   }

   oldRunfile = false;
   int i = 0;
   while (!fin.atEnd())
   {
      QString S = fin.readLine().trimmed();

      if (i == 0 && !S.contains("openLISEM"))
         oldRunfile = true;
      i++;

      if (S.contains("="))
      {

         QStringList SL = S.split(QRegExp("="));

         for (int j = 0; j < nrnamelist; j++)
         {
            if (namelist[j].name == SL[0].trimmed())
            {
               namelist[j].value = SL[1].trimmed();
               break;
            }
         }
         // VJ 110107 only read those variables that are in the coded definition to avoid
         //old junk that is read
      }
      // namelist now contains the actual runfile data
   }
}
//---------------------------------------------------------------------------
void lisemqt::ParseInputData()
{
   int j=0;
   bool dummyrain, dummysnow;

   // get all the options/checks
   for (j = 0; j < nrnamelist; j++)  //VJ 110107 changed to nrnamelist
   {
      int iii = namelist[j].value.toInt();
      QString p1 = namelist[j].name;
      QString p = namelist[j].value;
      bool check = iii == 1;
      if (p1.contains("["))
         continue;

      // main lisem types
      /*
          if (p1.compare("LISEM Type")==0)
          {
              SwitchWheelAsChannel = iii == LISEMWHEELTRACKS;
              SwitchMulticlass = iii == LISEMMULTICLASS;
              SwitchNutrients = iii == LISEMNUTRIENTS;
              SwitchGullies = iii == LISEMGULLIES;
          }
   */

      //options in the main code, order is not important
      if (p1.compare("No Erosion simulation")==0)          checkNoErosion->setChecked(check);
      if (p1.compare("Include main channels")==0)          checkIncludeChannel->setChecked(check);
      if (p1.compare("Include channel infil")==0)          checkChannelInfil->setChecked(check);
      if (p1.compare("Include channel baseflow")==0)       checkChannelBaseflow->setChecked(check);
      if (p1.compare("Include tile drains")==0)            checkIncludeTiledrains->setChecked(check);
      //if (p1.compare("All water and sediment to outlet")==0) checkAllinChannel->setChecked(check);

      if (p1.compare("Include Rainfall")==0)               dummyrain = check;//checkRainfall->setChecked(check);
      if (p1.compare("Include Snowmelt")==0)               dummysnow = check;//checkSnowmelt->setChecked(check);
      if (p1.compare("Alternative flow detachment")==0)    checkAltErosion->setChecked(check);
      if (p1.compare("Simple depression storage")==0)      checkSimpleDepression->setChecked(check);
      if (p1.compare("Hard Surfaces")==0)                  checkHardsurface->setChecked(check);
      if (p1.compare("Include buffers")==0)                checkBuffers->setChecked(check);
      if (p1.compare("Include Sediment traps")==0)         checkSedtrap->setChecked(check);
      if (p1.compare("Include wheeltracks")==0)            checkInfilCompact->setChecked(check);
      if (p1.compare("Include grass strips")==0)           checkInfilGrass->setChecked(check);
      if (p1.compare("Include crusts")==0)                 checkInfilCrust->setChecked(check);
      if (p1.compare("Impermeable sublayer")==0)           checkImpermeable->setChecked(check);
      //		if (p1.compare("Matric head files")==0)              checkDumphead->setChecked(check);
      if (p1.compare("Geometric mean Ksat")==0)            checkGeometric->setChecked(check);
      if (p1.compare("2nd back substitution")==0)          checkBacksubstitution->setChecked(check);
      //   if (p1.compare("Runoff maps in l/s/m")==0)           checkRunoffPerM->setChecked(check);
      if (p1.compare("Timeseries as PCRaster")==0)         checkWritePCRnames->setChecked(check);
      if (p1.compare("Timeplot as PCRaster")==0)           checkWritePCRtimeplot->setChecked(check);
      if (p1.compare("Regular runoff output")==0)          checkOutputTimeStep->setChecked(check);
      if (p1.compare("User defined output")==0)            checkOutputTimeUser->setChecked(check);
      if (p1.compare("No erosion at outlet")==0)           checkNoErosionOutlet->setChecked(check);
      //    if (p1.compare("Subsoil drainage")==0)               checkDrainage->setChecked(check);
      //    if (p1.compare("Gully infiltration")==0)             checkGullyInfil->setChecked(check);
      //    if (p1.compare("Use initial gully dimensions")==0)   checkGullyInit->setChecked(check);
      if (p1.compare("Report point output separate")==0)   checkSeparateOutput->setChecked(check);
      if (p1.compare("Report point output for SOBEK")==0)  checkSOBEKOutput->setChecked(check);
      if (p1.compare("SOBEK date string")==0)              SOBEKdatestring->setText(p);
      if (p1.compare("Sediment bulk density")==0)          E_BulkDens->setText(p);
      if (p1.compare("Use canopy storage map")==0)			  radioButton_9->setChecked(check);
      //checkInterceptionLAI->setChecked(!check);
      if (p1.compare("Canopy storage equation")==0)
      {
         //InterceptionEqNr = iii;
         switch (iii) {
         case 0 : radioButton_1->setChecked(true); break;
         case 1 : radioButton_2->setChecked(true); break;
         case 2 : radioButton_3->setChecked(true); break;
         case 3 : radioButton_4->setChecked(true); break;
         case 4 : radioButton_5->setChecked(true); break;
         case 5 : radioButton_6->setChecked(true); break;
         case 6 : radioButton_7->setChecked(true); break;
         case 7 : radioButton_8->setChecked(true); break;
         case 8 : radioButton_9->setChecked(true); break;
         }
      }
      if (p1.compare("Infil Method")==0 || p1.compare("Method")==0) //<= very old runfile
      {
         namelist[j].name = QString("Infil Method");
         uiInfilMethod = 0;
         switch(iii){
         case INFIL_SWATRE : uiInfilMethod = 1; break;
         case INFIL_GREENAMPT : uiInfilMethod = 2; break;
         case INFIL_GREENAMPT2 : uiInfilMethod = 2; checkInfil2layer->setChecked(true); break;
         case INFIL_SMITH : uiInfilMethod = 3; break;
         case INFIL_SMITH2 : uiInfilMethod = 3; checkInfil2layer->setChecked(true); break;
         case INFIL_KSAT : uiInfilMethod = 4; break;
         }
         E_InfiltrationMethod->setCurrentIndex(uiInfilMethod);
      }

      if (p1.compare("Ksat calibration")==0)
      {
         double val = namelist[j].value.toDouble();
         if (oldRunfile)
         {
            val/=100;
            QMessageBox::warning(this,"openLISEM",QString("Old runfile detected: calibration value Ksat changed from % to fraction:\n"
                                                          "Ksat calibration divided by 100, check 'Calibration' options in main menu."));
         }
         E_CalibrateKsat->setValue(val);
      }
      if (p1.compare("N calibration")==0)            E_CalibrateN->setValue(namelist[j].value.toDouble());
      if (p1.compare("Channel Ksat calibration")==0) E_CalibrateChKsat->setValue(namelist[j].value.toDouble());
      if (p1.compare("Channel N calibration")==0)    E_CalibrateChN->setValue(namelist[j].value.toDouble());
      if (p1.compare("Splash Delivery Ratio")==0)    E_SplashDelibery->setValue(namelist[j].value.toDouble());
      if (p1.compare("Stemflow fraction")==0)        E_StemflowFraction->setValue(namelist[j].value.toDouble());
      if (p1.compare("Canopy Openess")==0)        E_CanopyOpeness->setValue(namelist[j].value.toDouble());
      // VJ 110209 canopy openess, factor Aston as user input

      if (p1.compare("Output interval")==0)   printinterval->setValue(max(1,iii));

      if (p1.compare("Erosion map units (0/1/2)")==0)
      {
         int units = p.toInt();
         if (units == 0)
            checkUnits_tonha->setChecked(true);
         if (units == 1)
            checkUnits_kgcell->setChecked(true);
         if (units == 2)
            checkUnits_kgm2->setChecked(true);
      }


      if (p1.compare("CheckOutputMaps")==0)
      {
         outputcheck = p.split(",");
         checkBox_OutRunoff->setChecked(bool(outputcheck.at(0).toInt() == 1));
         checkBox_OutWH->setChecked(bool(outputcheck.at(2).toInt() == 1));
         checkBox_OutWHC->setChecked(bool(outputcheck.at(3).toInt() == 1));
         checkBox_OutInf->setChecked(bool(outputcheck.at(8).toInt() == 1));
         checkBox_OutV->setChecked(bool(outputcheck.at(7).toInt() == 1));
         checkBox_OutDet->setChecked(bool(outputcheck.at(5).toInt() == 1));
         checkBox_OutDep->setChecked(bool(outputcheck.at(6).toInt() == 1));
         checkBox_OutConc->setChecked(bool(outputcheck.at(1).toInt() == 1));
         checkBox_OutTC->setChecked(bool(outputcheck.at(4).toInt() == 1));
         checkBox_OutSurfStor->setChecked(bool(outputcheck.at(9).toInt() == 1));
         checkBox_OutChanVol->setChecked(bool(outputcheck.at(10).toInt() == 1));
         // checkboxes normal output map series, numbering according to original LISEM
      }
   }

   if (!dummyrain && ! dummysnow)
      QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");

   checkRainfall->setChecked(dummyrain);
   checkSnowmelt->setChecked(dummysnow);

   // get directory and file names
   for (j = 0; j < nrnamelist; j++)//VJ 110107 changed to nrnamelist
   {
      QString p1 = namelist[j].name;
      QString p = namelist[j].value;

      if (p1.compare("Begin time")==0) E_BeginTime->setText(p);
      if (p1.compare("End time")==0) E_EndTime->setText(p);
      if (p1.compare("Timestep")==0) E_Timestep->setText(p);
      if (p1.compare("SWATRE internal minimum timestep")==0)
      {
         swatreDT = p.toDouble();
         E_SWATREDtsecFraction->setValue(swatreDT/E_Timestep->text().toDouble());
      }

      // input ourput dirs and file names
      if (p1.compare("Map Directory")==0) E_MapDir->setText(CheckDir(p));
      if (p1.compare("Result Directory")==0) E_ResultDir->setText(CheckDir(p));

      if (p1.compare("Main results file")==0) E_MainTotals->setText(p);
      if (p1.compare("Filename point output")==0) E_PointResults->setText(p);
      if (p1.compare("Filename landunit output")==0) E_LandunitResults->setText(p);
      // resultDir is added in report operation

      if (p1.compare("Rainfall Directory")==0) RainFileDir = CheckDir(p);
      if (p1.compare("Rainfall file")==0)
      {
         E_RainfallName->setText(p);
         RainFileName = /*rainFileDir + */E_RainfallName->text();
      }

      if (p1.compare("Erosion map")==0) E_DetachmentMap->setText(p);
      if (p1.compare("Deposition map")==0) E_DepositionMap->setText(p);
      if (p1.compare("Soilloss map")==0) E_SoillossMap->setText(p);
      // resultDir is added in report operation
      //NO checking
      //		if (checkSnowmelt->isChecked())
      //	{
      if (p1.compare("Snowmelt Directory")==0) SnowmeltFileDir = CheckDir(p);
      if (p1.compare("Snowmelt file")==0)
      {
         E_SnowmeltName->setText(p);
         SnowmeltFileName = /*SnowmeltFileDir + */E_SnowmeltName->text();
      }
      //}

      if (p1.compare("Table Directory")==0)
      {
         SwatreTableDir = CheckDir(p);
         if (SwatreTableDir.isEmpty())
            SwatreTableDir = E_MapDir->text();
         E_SwatreTableDir->setText(SwatreTableDir);
      }

      if (p1.compare("Table File")==0)
      {
         SwatreTableName = p;
         E_SwatreTableName->setText(SwatreTableName);
      }
   }

   if (SwatreTableName.isEmpty())
   {
      SwatreTableName = E_MapDir->text() + QString("profile.inp");
      E_SwatreTableName->setText(SwatreTableName);
   }

   //****====------====****//
   // get all map names, DEFmaps contains default map names and descriptions
   // adapt the DEFmaps list with names from the run file
   // this is to display the correct names in the interface
   for (j = mapstartnr; j < nrnamelist; j++)  //VJ 110107 changed to nrnamelist
   {
      for (int i = 0; i < DEFmaps.size(); i++)
      {
         QStringList S = DEFmaps.at(i).split(";",QString::SkipEmptyParts);
         if (S.contains(namelist[j].name))
         {
            QFileInfo fil(namelist[j].value);
            S[2] = fil.fileName();  //VJ bug fix from 4 to 2
            namelist[j].value = fil.fileName();
            // replace namelist string with filename only
            // some runfiles have the complete pathname
            DEFmaps.replace(i, S.join(";") );
         }
      }
   }

   // strip pathname from output filename
   for (int j = 0; j < nrnamelist; j++)
      if (namelist[j].name.startsWith("OUT"))
      {
         QFileInfo fil(namelist[j].value);
         namelist[j].value = fil.fileName();
      }


   // fill the mapList structure with all map names fom the runfile
   // if there are new variables that are not in the run file
   // the maplist contains he default names already
   // this is to get the correct names for the model run
   fillNamelistMapnames(false);
   for (int k = 0; k < nrmaplist; k++)
       mapList[k].dir = E_MapDir->text();
   // dir necessary?

   //RunAllChecks();
   //obsolete: is done elsewhere
}
//---------------------------------------------------------------------------
QString lisemqt::CheckDir(QString p)
{
   /* TODO mulitplatform: fromNativeSeparators etc*/
   p.replace("/","\\");
   if (!p.endsWith("\\"))
      p = p + "\\";
   if (!QDir(p).exists())
      p.clear();

   return p;
}
//---------------------------------------------------------------------------
// change runfile strings with current interface options, called by savefile
// savefile is called just before the model is run with tmp runfile
void lisemqt::updateModelData()
{
   if(!checkRainfall->isChecked() && !checkSnowmelt->isChecked())
      QMessageBox::warning(this,"openLISEM","No rainfall or snowmelt, running on empty!");


   for (int j = 0; j < nrnamelist; j++)
   {
      QString p1 = namelist[j].name;
      QString p;
      if (p1.compare("No Erosion simulation")==0) 			  namelist[j].value.setNum((int)checkNoErosion->isChecked());
      if (p1.compare("Include main channels")==0) 			  namelist[j].value.setNum((int)checkIncludeChannel->isChecked());
      if (p1.compare("Include channel infil")==0)          namelist[j].value.setNum((int)checkChannelInfil->isChecked());
      if (p1.compare("Include channel baseflow")==0)       namelist[j].value.setNum((int)checkChannelBaseflow->isChecked());
      if (p1.compare("Include tile drains")==0)            namelist[j].value.setNum((int)checkIncludeTiledrains->isChecked());
      if (p1.compare("Include Rainfall")==0)               namelist[j].value.setNum((int)checkRainfall->isChecked());
      if (p1.compare("Include Snowmelt")==0)               namelist[j].value.setNum((int)checkSnowmelt->isChecked());
      if (p1.compare("Alternative flow detachment")==0)    namelist[j].value.setNum((int)checkAltErosion->isChecked());
      if (p1.compare("Simple depression storage")==0)      namelist[j].value.setNum((int)checkSimpleDepression->isChecked());
      if (p1.compare("Hard Surfaces")==0)                  namelist[j].value.setNum((int)checkHardsurface->isChecked());
      if (p1.compare("Include buffers")==0)                namelist[j].value.setNum((int)checkBuffers->isChecked());
      if (p1.compare("Include Sediment traps")==0)         namelist[j].value.setNum((int)checkSedtrap->isChecked());
      if (p1.compare("Include wheeltracks")==0)            namelist[j].value.setNum((int)checkInfilCompact->isChecked());
      if (p1.compare("Include grass strips")==0)           namelist[j].value.setNum((int)checkInfilGrass->isChecked());
      if (p1.compare("Include crusts")==0)                 namelist[j].value.setNum((int)checkInfilCrust->isChecked());
      if (p1.compare("Impermeable sublayer")==0)           namelist[j].value.setNum((int)checkImpermeable->isChecked());
      //if (p1.compare("Matric head files")==0)              namelist[j].value.setNum((int)checkDumphead->isChecked());
      if (p1.compare("Geometric mean Ksat")==0)            namelist[j].value.setNum((int)checkGeometric->isChecked());
      if (p1.compare("2nd back substitution")==0)          namelist[j].value.setNum((int)checkBacksubstitution->isChecked());
      if (p1.compare("Timeseries as PCRaster")==0)         namelist[j].value.setNum((int)checkWritePCRnames->isChecked());
      if (p1.compare("Timeplot as PCRaster")==0)           namelist[j].value.setNum((int)checkWritePCRtimeplot->isChecked());
      if (p1.compare("Regular runoff output")==0)          namelist[j].value.setNum((int)checkOutputTimeStep->isChecked());
      if (p1.compare("User defined output")==0)            namelist[j].value.setNum((int)checkOutputTimeUser->isChecked());
      if (p1.compare("No erosion at outlet")==0)           namelist[j].value.setNum((int)checkNoErosionOutlet->isChecked());
      if (p1.compare("Report point output separate")==0)   namelist[j].value.setNum((int)checkSeparateOutput->isChecked());
      if (p1.compare("Report point output for SOBEK")==0)  namelist[j].value.setNum((int)checkSOBEKOutput->isChecked());
      if (p1.compare("SOBEK date string")==0)              namelist[j].value = SOBEKdatestring->text();
      if (p1.compare("Sediment bulk density")==0)          namelist[j].value = E_BulkDens->text();
      //if (p1.compare("Use canopy storage map")==0)   	     namelist[j].value.setNum((int)!checkInterceptionLAI->isChecked());
      if (p1.compare("Canopy storage equation")==0)
      {
         int i = 0;
         if(radioButton_1->isChecked()) i = 0;
         if(radioButton_2->isChecked()) i = 1;
         if(radioButton_3->isChecked()) i = 2;
         if(radioButton_4->isChecked()) i = 3;
         if(radioButton_5->isChecked()) i = 4;
         if(radioButton_6->isChecked()) i = 5;
         if(radioButton_7->isChecked()) i = 6;
         if(radioButton_8->isChecked()) i = 7;
         if(radioButton_9->isChecked()) i = 8;
         namelist[j].value.setNum(i);
      }

      if (p1.compare("Begin time")==0) namelist[j].value = E_BeginTime->text();
      if (p1.compare("End time")==0)   namelist[j].value = E_EndTime->text();
      if (p1.compare("Timestep")==0)   namelist[j].value = E_Timestep->text();
      if (p1.compare("Map Directory")==0)    namelist[j].value = E_MapDir->text();
      if (p1.compare("Result Directory")==0) namelist[j].value = E_ResultDir->text();
      if (p1.compare("Main results file")==0) namelist[j].value = E_MainTotals->text();
      if (p1.compare("Filename point output")==0) namelist[j].value = E_PointResults->text();
      if (p1.compare("Filename landunit output")==0) namelist[j].value = E_LandunitResults->text();
      if (p1.compare("Rainfall Directory")==0) namelist[j].value = RainFileDir;
      if (p1.compare("Rainfall file")==0) namelist[j].value = E_RainfallName->text();
      if (p1.compare("Erosion map")==0) namelist[j].value = E_DetachmentMap->text();
      if (p1.compare("Deposition map")==0) namelist[j].value = E_DepositionMap->text();
      if (p1.compare("Soilloss map")==0) namelist[j].value = E_SoillossMap->text();
      if (p1.compare("Snowmelt Directory")==0) namelist[j].value = SnowmeltFileDir;
      if (p1.compare("Snowmelt file")==0) namelist[j].value = E_SnowmeltName->text();
      if (p1.compare("Ksat calibration")==0) namelist[j].value = E_CalibrateKsat->text();
      if (p1.compare("N calibration")==0) namelist[j].value = E_CalibrateN->text();
      if (p1.compare("Channel Ksat calibration")==0) namelist[j].value = E_CalibrateChKsat->text();
      if (p1.compare("Channel N calibration")==0) namelist[j].value = E_CalibrateChN->text();
      if (p1.compare("Splash Delivery Ratio")==0) namelist[j].value = E_SplashDelibery->text();      
      if (p1.compare("Stemflow fraction")==0) namelist[j].value = E_StemflowFraction->text();
      if (p1.compare("Canopy Openess")==0) namelist[j].value = E_CanopyOpeness->text();
      // VJ 110209 canopy openess, factor Aston as user input

      if (p1.compare("Output interval")==0) namelist[j].value = printinterval->cleanText();
      if (p1.compare("Regular runoff output")==0) namelist[j].value.setNum(1);
      if (p1.compare("User defined output")==0) namelist[j].value.setNum(0);
      if (p1.compare("Output times")==0) namelist[j].value.setNum(0);
      //TODO fix output stuff

      if (p1.compare("Table Directory")==0) namelist[j].value = E_SwatreTableDir->text();//setTextSwatreTableDir;
      if (p1.compare("Table File")==0) namelist[j].value = E_SwatreTableName->text();//SwatreTableName;
      if (p1.compare("SWATRE internal minimum timestep")==0)
      {
         double fraction = E_SWATREDtsecFraction->value();
         swatreDT = E_Timestep->text().toDouble()*fraction;
         namelist[j].value.setNum(swatreDT,'g',6);
      }

      if (p1.compare("Infil Method")==0)
      {
         switch(uiInfilMethod)
         {
         case 0 : namelist[j].value.setNum(INFIL_NONE); break;
         case 1 : namelist[j].value.setNum(INFIL_SWATRE);break;
         case 2 : if(checkInfil2layer->isChecked()) namelist[j].value.setNum(INFIL_GREENAMPT2);
            else namelist[j].value.setNum(INFIL_GREENAMPT);break;
         case 3 : if(checkInfil2layer->isChecked()) namelist[j].value.setNum(INFIL_SMITH2);
            else namelist[j].value.setNum(INFIL_SMITH); break;
         case 4: namelist[j].value.setNum(INFIL_KSAT); break;
         }
      }

      if (p1.compare("Erosion map units (0/1/2)")==0)
      {
         if (checkUnits_tonha->isChecked())
            namelist[j].value.setNum(0);
         if (checkUnits_kgcell->isChecked())
            namelist[j].value.setNum(1);
         if (checkUnits_kgm2->isChecked())
            namelist[j].value.setNum(2);
      }
      //VJ 110110 added

      if (p1.compare("CheckOutputMaps")==0)
      {
         outputcheck.clear();
         if (			checkBox_OutRunoff->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			  checkBox_OutConc->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			    checkBox_OutWH->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (	    	   checkBox_OutWHC->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (	          checkBox_OutTC->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			   checkBox_OutDet->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			   checkBox_OutDep->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			     checkBox_OutV->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (			   checkBox_OutInf->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (	    checkBox_OutSurfStor->isChecked()) outputcheck << "1"; else outputcheck << "0";
         if (       checkBox_OutChanVol->isChecked()) outputcheck << "1"; else outputcheck << "0";
         namelist[j].value = outputcheck.join(",");
      }
      //namelist[j].value = p;
   }

   //get all actual mapnames from the mapList structure
   fillNamelistMapnames(true);

}
//---------------------------------------------------------------------------
