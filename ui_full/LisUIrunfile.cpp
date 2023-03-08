/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
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

#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "global.h"
#include "option.h"

//---------------------------------------------------------------------------
//! fill namelist with the actual runfile data but correct for old runfiles
//! so that faulty data or obsolete vars are ignored
//! This function reads the runfile and checks against default names and descriptions
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

    //   currentDir = QFileInfo(op.runfilename).path();//absoluteFilePath();
    //   QDir::setCurrent(currentDir);

    // read all lines in the runfile BUT
    // each line is compared to the hardcoded name in lisemqt::defaultRunFile()
    // so old junk or misspelled stuff is simply IGNORED
    // namelist now contains the actual runfile data
    oldRunfile = false;
    saveRunFileOnce=false;
    int i = 0;
    //int found = 0;

    while (!fin.atEnd())
    {
        QString S = fin.readLine().trimmed();

        if (i == 0 && !S.contains("openLISEM"))
            oldRunfile = true;
        if (i == 0 && !S.contains("[openLISEM runfile version 6.0]"))
            saveRunFileOnce = true;

        i++;
        if (S.contains("="))
        {
            QStringList SL = S.split(QRegExp("="));

            for (int j = 0; j < nrnamelist; j++) {
                if (namelist[j].name == SL[0].trimmed()) {
                    namelist[j].value = SL[1].trimmed();
                    namelist[j].gotit = true;
                    break;
                }
            }
        }
    }

    for (int i = 0; i < nrnamelist; i++) {
 //      if (!namelist[i].value.isEmpty() && namelist[i].gotit)
   //       qDebug() << namelist[i].name << namelist[i].value << namelist[i].gotit;
        if (!namelist[i].value.isEmpty() && !namelist[i].gotit)
            saveRunFileOnce = true;
    }
}
//---------------------------------------------------------------------------
//! ParseInputData : interpret runfile text and fill interface variables
void lisemqt::ParseInputData()
{
    int j=0;
    int dummykinwave = 1;
    // get all the options/checks

    resetAll();
    bool ETmaps = false;
    bool Rainmaps = false;

    QLocale loc = QLocale::system(); // current locale
    QChar pnt = loc.decimalPoint();

    for (j = 0; j < nrnamelist; j++)
    {

        QString p = namelist[j].value;
        QString p1 = namelist[j].name;
        bool ok;

        int iii = loc.toInt(p);
        double valc = loc.toDouble(p,&ok); //namelist[j].value.toDouble(&ok);
        if (!ok && pnt == ',') {
            QString p2 = p.replace(",",".");
            valc = p2.toDouble();
        }
        if (!ok && pnt == '.') {
            QString p2 = p.replace(",",".");
            valc = p2.toDouble();
        }

    //  qDebug() << p1 <<  p  << iii <<  valc;

        bool check = iii == 1;
        if (p1.contains("["))
            continue;
        if (p1.compare("Event based")==0)                   checkEventBased->setChecked(check);
        if (p1.compare("Use Rainfall maps")==0)             Rainmaps = check;
        if (p1.compare("Rainfall ID interpolation")==0)     checkIDinterpolation->setChecked(check);
        if (p1.compare("IDI factor")==0)                    E_IDIfactor->setValue(valc);
        if (p1.compare("Rainfall Bias Correction")==0)      E_biasCorrectionP->setValue(valc);

        if (p1.compare("Include ET")==0)                    checkIncludeET->setChecked(check);
        if (p1.compare("Use ET maps")==0)                   ETmaps = check;
        if (p1.compare("Daily ET")==0)                      checkDailyET->setChecked(check);
        if (p1.compare("ET Bias Correction")==0)            E_biasCorrectionET->setValue(valc);
        if (p1.compare("Rainfall ET threshold")==0)         E_rainfallETA_threshold->setValue(valc);

        //if (p1.compare("Include Snowmelt")==0)            checkSnowmelt->setChecked(check);

        if (p1.compare("Nr user Cores")==0) nrUserCores->setValue(iii);
        if (p1.compare("Result datetime")==0) checkAddDatetime->setChecked(check);
        if (p1.compare("Timeplot as PCRaster")==0)           checkWritePCRaster->setChecked(!check);
        if (p1.compare("Report point output separate")==0)   checkSeparateOutput->setChecked(check);
        //       if (p1.compare("Report point output for SOBEK")==0)  checkWriteSOBEK->setChecked(check);
        if (p1.compare("Report discharge units")==0)
        {
            int units = p.toInt();
            if (units == 0)
                checkUnits_ls->setChecked(true);
            if (units == 1)
                checkUnits_m3s->setChecked(true);
        }

        if (p1.compare("Report digits out")==0)              E_DigitsOut->setValue(iii);
        if (p1.compare("Report format GTiff")==0)            checkFormatGtiff->setChecked(check);
        if (p1.compare("End run report")==0)                 checkEndRunReport->setChecked(check);
      //  if (p1.compare("Add timestamp")==0)                  checkOutputTimestamp->setChecked(check);
      //  if (p1.compare("Variable Timestep")==0)              checkVariableTimestep->setChecked(check);
        if (p1.compare("Include Satellite Image")==0)        checksatImage->setChecked(check);
        if (p1.compare("Output interval")==0)                printinterval->setValue(std::max(1,iii));
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

        // GENERAL
        if (p1.compare("Include main channels")==0)          checkIncludeChannel->setChecked(check);
        if (p1.compare("Include channel infil")==0)          checkChannelInfil->setChecked(check);
        if (p1.compare("Include channel baseflow")==0)       checkChannelBaseflow->setChecked(check);
        if (p1.compare("Include stationary baseflow")==0)       checkStationaryBaseflow->setChecked(check);
      //  if (p1.compare("Adjust channel crosssection")==0)     checkChannelAdjustCHW->setChecked(check);
        if (p1.compare("GW flow explicit")==0)              checkGWflowexplicit->setChecked(check);
        if (p1.compare("GW flow SWAT")==0)                  checkGWflowSWAT->setChecked(check);
        if (p1.compare("GW recharge factor")==0)            GW_recharge->setValue(valc);
        if (p1.compare("GW flow factor")==0)                GW_flow->setValue(valc);
        if (p1.compare("GW river inflow factor")==0)           GW_inflow->setValue(valc);
        if (p1.compare("GW slope factor")==0)               GW_slope->setValue(valc);
      //  if (p1.compare("GW lag factor")==0)               GW_lag->setValue(valc);
        if (p1.compare("GW deep percolation")==0)           GW_deep->setValue(valc);
        if (p1.compare("GW threshold factor")==0)           GW_threshold->setValue(valc);
        if (p1.compare("GW initial level")==0)              GW_initlevel->setValue(valc);

        if (p1.compare("Include channel culverts")==0)      checkChannelCulverts->setChecked(check);
        if (p1.compare("Include channel inflow")==0)        checkChannelInflow->setChecked(check);
        if (p1.compare("Include Erosion simulation")==0)    checkDoErosion->setChecked(check);
        if (p1.compare("Include road system")==0)           checkRoadsystem->setChecked(check);
        if (p1.compare("Include storm drains")==0)          checkStormDrains->setChecked(check);
        if (p1.compare("Hard Surfaces")==0)                 checkHardsurface->setChecked(check);
        if (p1.compare("Include house storage")==0)         checkHouses->setChecked(check);
        if (p1.compare("Include buildings")==0)             checkHouses->setChecked(check);

        if (p1.compare("Add buildings to DEM")==0)          checkAddBuildingDEM->setChecked(check);
        if (p1.compare("Add building fraction")==0)         E_AddBuildingFraction->setValue(valc);
        if (p1.compare("Include raindrum storage")==0)      checkRaindrum->setChecked(check);
        if (p1.compare("Include tile drains")==0)           checkIncludeTiledrains->setChecked(check);

        // INTERCEPTION
        if (p1.compare("Canopy storage equation")==0)
        {
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
        if (p1.compare("Include litter interception")==0)    checkIncludeLitter->setChecked(check);
        if (p1.compare("Litter interception storage")==0)    E_LitterSmax->setValue(valc);
        //if (p1.compare("Stemflow fraction")==0)        E_StemflowFraction->setValue(valc);
        if (p1.compare("Canopy Openess")==0)           E_CanopyOpeness->setValue(valc);
        // VJ 110209 canopy openess, factor Aston as user input

        // INFILTRATION
        if (p1.compare("Infil Method")==0) {
            switch(iii)
            {
            case INFIL_NONE : E_InfiltrationMethod->setCurrentIndex(0);break;
            case INFIL_SWATRE : E_InfiltrationMethod->setCurrentIndex(1);break;
            case INFIL_GREENAMPT : E_InfiltrationMethod->setCurrentIndex(2);break;
            case INFIL_SMITH : E_InfiltrationMethod->setCurrentIndex(3); break;
            }
        }


        if (p1.compare("Include compacted")==0)              checkInfilCompact->setChecked(check);
        if (p1.compare("Include grass strips")==0)           checkInfilGrass->setChecked(check);
        if (p1.compare("Grassstrip Mannings n")==0)           E_GrassStripN->setValue(valc);
        if (p1.compare("Include crusts")==0)                 checkInfilCrust->setChecked(check);
        if (p1.compare("Impermeable sublayer")==0)           checkInfilImpermeable->setChecked(check);
        if (p1.compare("Two layer")==0)                      checkInfil2layer->setChecked(check);
        if (p1.compare("Geometric mean Ksat")==0)            checkGeometric->setChecked(check);
        //	  if (p1.compare("Matric head files")==0)              checkDumphead->setChecked(check);
        if (p1.compare("Sediment trap Mannings n")==0)           E_SedTrapN->setValue(valc);

        //VJ 111120 water repellency
//        if (p1.compare("Use Water Repellency")==0)      checkWaterRepellency->setChecked(check);
//        if (p1.compare("Water Repellency A")==0)        E_waterRep_a->setValue(valc);
//        if (p1.compare("Water Repellency B")==0)        E_waterRep_b->setValue(valc);
//        if (p1.compare("Water Repellency C")==0)        E_waterRep_c->setValue(valc);
//        if (p1.compare("Water Repellency D")==0)        E_waterRep_d->setValue(valc);


        // FLOW
        if (p1.compare("Include flow barriers")==0)          checkFlowBarriers->setChecked(check);
        if (p1.compare("Flow barrier table filename")==0)    line_FlowBarriers->setText(p);
        if (p1.compare("Include buffers")==0)          checkBuffers->setChecked(check);
        if (p1.compare("Minimum reported flood height")==0)  E_floodMinHeight->setValue(valc);
        if (p1.compare("Flooding courant factor")==0)        E_courantFactor->setValue(valc);
        if (p1.compare("Flood solution")==0)
        {
            switch (iii) {
            case 0 : rb_swof2->setChecked(true); break;
            case 1 : rb_swofnomuscl->setChecked(true); break;
            case 2 : rb_swofmuscl->setChecked(true); break;
            }
        }
        if (p1.compare("Routing Kin Wave 2D")==0)            dummykinwave = iii;
        if (p1.compare("Flow Boundary 2D")==0)               E_FlowBoundary->setValue(iii);
        if (p1.compare("Correct DEM")==0)                    checkCorrectDem->setChecked(check);
        if (p1.compare("Use 2D Diagonal flow")==0)           check2DDiagonalFlow->setChecked(check);
        if (p1.compare("Use 2D Diagonal flow new")==0)       check2DDiagonalFlowNew->setChecked(check);
        if (p1.compare("Pit Value")==0)                      E_pitValue->setValue(valc);
        if (p1.compare("Flood initial level map")==0)        checkFloodInitial->setChecked(check);
        if (p1.compare("Timestep flood")==0)                 E_TimestepMinFlood->setValue(valc);
    //    if (p1.compare("Use SWOF watersheds")==0)            checkSWOFWatersheds->setChecked(check);


        // EROSION
        if (p1.compare("Splash equation")==0)                E_splashEquation->setValue(iii);
        if (p1.compare("Detachment efficiency")==0)          E_EfficiencyDET->setCurrentIndex(iii-1);
        if (p1.compare("Detachment efficiency channel")==0)  E_EfficiencyDETCH->setCurrentIndex(iii-1);
        if (p1.compare("Direct efficiency channel")==0)      E_EfficiencyDirect->setValue(valc);
        if (p1.compare("Settling Velocity")==0)              E_settlingVelocity->setCurrentIndex(iii-1);
        if (p1.compare("Include Sediment traps")==0)         checkSedtrap->setChecked(check);
        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
//            if (p.contains(",")){
//                p.replace(",",";");
//            }
            param = p.split(";",Qt::SkipEmptyParts);
            //qDebug() << p << param;
            if (param.count() > 4) {
                bool ok;
                int i = param[0].toInt(&ok, 10);
                param.clear();
                if (i == 1) param << "1"; else param << "0";
                param <<  "28.3" << "0.52" << "0.042";
            }
            radioButtonKE1->setChecked(param[0].toInt());
            spinKEparameterA1->setValue(param[1].toDouble());
            spinKEparameterB1->setValue(param[2].toDouble());
            spinKEparameterC1->setValue(param[3].toDouble());
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
//            if (p.contains(",")){
//                p.replace(",",";");
//            }
            param = p.split(";",Qt::SkipEmptyParts);
            if (param.count() > 3) {
                bool ok;
                int i = param[0].toInt(&ok, 10);
                param.clear();
                if (i == 1) param << "1"; else param << "0";
                param <<  "8.95" << "8.44";
            }
            radioButtonKE2->setChecked(param[0].toInt());
            spinKEparameterA2->setValue(param[1].toDouble());
            spinKEparameterB2->setValue(param[2].toDouble());
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
//            if (p.contains(",")){
//                p.replace(",",";");
//            }
            param = p.split(";",Qt::SkipEmptyParts);
            if (param.count() > 3) {
                bool ok;
                int i = param[0].toInt(&ok, 10);
                param.clear();
                if (i == 1) param << "1"; else param << "0";
                param <<  "7.6" << "0.22";
            }
            radioButtonKE3->setChecked(param[0].toInt());
            spinKEparameterA3->setValue(param[1].toDouble());
            spinKEparameterB3->setValue(param[2].toDouble());
        }
        if (p1.compare("KE time based")==0)                   checkKETimebased->setChecked(check);

        if (p1.compare("Use material depth")==0)              checkMaterialDepth->setChecked(check);
        if (p1.compare("No detachment boundary")==0)          checkNoSedBoundary->setChecked(check);
        if (p1.compare("Use 2 phase flow")==0)                checkSed2Phase->setChecked(check);
        if (p1.compare("River BL method")==0)                 E_RBLMethod->setCurrentIndex(iii-1);
        if (p1.compare("River SS method")==0)                 E_RSSMethod->setCurrentIndex(iii-1);
        if (p1.compare("Include River diffusion")==0)         checkDiffusionCH->setChecked(check);
        if (p1.compare("Flooding BL method")==0)              E_BLMethod->setCurrentIndex(iii-1);
        if (p1.compare("Flooding SS method")==0)              E_SSMethod->setCurrentIndex(iii-1);
        if (p1.compare("Include diffusion")==0)               checkDiffusion->setChecked(check);
        if (p1.compare("Sigma diffusion")==0)                 E_SigmaDiffusion->setValue(valc);

//        if (p1.compare("Use grain size distribution")==0)     checkSedMultiGrain->setChecked(check);
//        if (p1.compare("Estimate grain size distribution")==0)checkEstimateGrainSizeDistribution->setChecked(check);
//        if (p1.compare("Read grain distribution maps")==0)    checkReadGrainSizeDistribution->setChecked(check);
//        if (p1.compare("Number of grain size classes (simulated)")==0)
//                                                              E_NumberClasses->setValue(iii);
//        if (p1.compare("Grain size class maps")==0)   {
//            if (p.contains(",")) p.replace(",",";");
//            E_GrainSizes->setText(p);
//        }
        //ADVANCED
        if (p1.compare("Advanced Options")==0)                 checkAdvancedOptions->setChecked(check);
      //  if (p1.compare("Calculate erosion inside 2D loop")==0)  checkErosionInsideLoop->setChecked(check);
        if (p1.compare("Use linked list")==0)        checkLinkedList->setChecked(check);
        if (p1.compare("Flooding SWOF flux limiter")==0)     E_FloodFluxLimiter->setValue(iii);
        if (p1.compare("Flooding SWOF Reconstruction")==0)   E_FloodReconstruction->setValue(iii);
        if (p1.compare("Use time avg V")==0)                 checkTimeavgV->setChecked(check);
        if (p1.compare("Flooding mixing coefficient")==0)    E_mixingFactor->setValue(valc);
        if (p1.compare("Flooding runoff partitioning")==0)   E_runoffPartitioning->setValue(valc);
        if (p1.compare("Flood max iterations")==0)           E_FloodMaxIter->setValue(iii);
        if (p1.compare("Min WH flow")==0)                    E_minWHflow->setText(p);
        if (p1.compare("Use gravity flow")==0)               checkGravityToChannel->setChecked(check);
     //   if (p1.compare("Angle flow to channel")==0)          E_angleToChannel->setValue(valc);
       // if (p1.compare("Use fixed angle")==0)                checkFixedAngle->setChecked(check);
        if (p1.compare("Use Channel Kinwave dt")==0)         checkKinWaveChannel->setChecked(check);
        if (p1.compare("Channel KinWave dt")==0)             E_ChannelKinWaveDt->setValue(valc);
        if (p1.compare("Use Channel Max V")==0)         checkChanMaxVelocity->setChecked(check);
        if (p1.compare("Channel Max V")==0)             E_chanMaxVelocity->setValue(valc);

        //CALIBRATION
        if (p1.compare("Smax calibration")==0)         E_CalibrateSmax->setValue(valc);
        if (p1.compare("Ksat calibration")==0)         E_CalibrateKsat->setValue(valc);
        if (p1.compare("Ksat2 calibration")==0)         E_CalibrateKsat2->setValue(valc);
        if (p1.compare("Grain Size calibration D50")==0)   E_CalibrateD50->setValue(valc);
        if (p1.compare("Grain Size calibration D90")==0)   E_CalibrateD90->setValue(valc);
        if (p1.compare("N calibration")==0)            E_CalibrateN->setValue(valc);
        if (p1.compare("Theta calibration")==0)        E_CalibrateTheta->setValue(valc);
        if (p1.compare("Psi calibration")==0)          E_CalibratePsi->setValue(valc);
        if (p1.compare("Channel Ksat calibration")==0) E_CalibrateChKsat->setValue(valc);
        if (p1.compare("Channel N calibration")==0)    E_CalibrateChN->setValue(valc);
        if (p1.compare("Channel tortuosity")==0)    E_CalibrateChTor->setValue(valc);
        if (p1.compare("Cohesion calibration")==0)     E_CalibrateCOH->setValue(valc);
        if (p1.compare("Cohesion Channel calibration")==0)    E_CalibrateCHCOH->setValue(valc);
        if (p1.compare("Ucr Channel calibration")==0)    E_CalibrateCHUcr->setValue(valc);
        if (p1.compare("SV calibration")==0)    E_CalibrateCHSV->setValue(valc);
        if (p1.compare("Aggregate stability calibration")==0)    E_CalibrateAS->setValue(valc);
        if (p1.compare("Splash Delivery Ratio")==0)    E_SplashDelibery->setValue(valc);
        if (p1.compare("Particle Cohesion of Deposited Layer")==0) E_DepositedCohesion->setValue(valc);
        if (p1.compare("Sediment bulk density")==0)          E_BulkDens->setValue(valc);


        // STANDARD OUTPUT FILES
        if (p1.compare("OutRunoff")==0)         checkBox_OutRunoff->setChecked(check);
        if (p1.compare("OutWH")==0)             checkBox_OutWH->setChecked(check);
        if (p1.compare("OutV")==0)              checkBox_OutV->setChecked(check);
        if (p1.compare("OutInterception")==0)  checkBox_OutInterception->setChecked(check);
        if (p1.compare("OutSurfStor")==0)       checkBox_OutSurfStor->setChecked(check);
        if (p1.compare("OutInf")==0)            checkBox_OutInf->setChecked(check);
        if (p1.compare("OutTileDrain")==0)      checkBox_OutTiledrain->setChecked(check);
        if (p1.compare("OutTileVolume")==0)         checkBox_OutTileVol->setChecked(check);
        if (p1.compare("OutTheta")==0)         checkBox_OutTheta->setChecked(check);
        if (p1.compare("OutDet")==0)     checkBox_OutDet->setChecked(check);
        if (p1.compare("OutDep")==0)     checkBox_OutDep->setChecked(check);
        if (p1.compare("OutTC")==0)      checkBox_OutTC->setChecked(check);
        if (p1.compare("OutConc")==0)    checkBox_OutConc->setChecked(check);
        if (p1.compare("OutSed")==0)     checkBox_OutSed->setChecked(check);
        if (p1.compare("OutSL")==0)      checkBox_OutSL->setChecked(check);
        if (p1.compare("OutSedSS")==0)     checkBox_OutSedSS->setChecked(check);
        if (p1.compare("OutSedBL")==0)     checkBox_OutSedBL->setChecked(check);

   }

    // ###################################

    radioETfile->setChecked(!ETmaps);
    radioETsatfile->setChecked(ETmaps);
    radioRainFile->setChecked(!Rainmaps);
    radioRainSatFile->setChecked(Rainmaps);

    groupAdvanced->setVisible(checkAdvancedOptions->isChecked());

    on_checkIncludeET_toggled(checkIncludeET->isChecked());

    if (checkSedtrap->isChecked())
        on_checkSedtrap_clicked();
    if (checkInfilGrass->isChecked())
        on_checkInfilGrass_clicked();

    E_SigmaDiffusion->setEnabled(checkDiffusion->isChecked());

    checkOverlandFlow1D->setChecked(dummykinwave == 1);
    checkOverlandFlow2Ddyn->setChecked(dummykinwave == 3);
    checkOverlandFlow2Dkindyn->setChecked(dummykinwave == 4);
    setFloodTab(true);//dummykinwave > 1);
    setErosionTab();

    // first guess
    E_WorkDir = QFileInfo(E_runFileList->currentText()).dir().absolutePath();
    QDir dir(E_WorkDir);
    if (dir.cdUp())
        E_WorkDir = dir.absolutePath()+"/";
    // workdir is now parent of runfile directory
   // qDebug() << E_WorkDir;

    QString daystart, minstart, dayend, minend;
    for (j = 0; j < nrnamelist; j++)
    {
        QString p1 = namelist[j].name;
        QString p = namelist[j].value;
//qDebug() << j << namelist[j].name << namelist[j].value;
        if (p1.compare("Begin time day")==0) daystart = p;//E_BeginTimeDay->setText(p);
        if (p1.compare("Begin time")==0) minstart = p;//E_BeginTimeMin->setText(p);
        if (p1.compare("End time day")==0)   dayend = p;//E_EndTimeDay->setText(p);
        if (p1.compare("End time")==0)   minend = p;//E_EndTimeMin->setText(p);
        if (p1.compare("Timestep")==0) E_Timestep->setText(p);

        // input output dirs and file names
        if (p1.compare("Map Directory")==0)
        {
            E_MapDir->setText(CheckDir(p));

            if (QFileInfo(E_MapDir->text()).exists())
            {
                E_WorkDir = E_MapDir->text();
                QDir dir(E_WorkDir);
                if (dir.cdUp())
                    E_WorkDir = dir.absolutePath()+"/";
                // workdir is now parent of maps directory
            }

            if (E_MapDir->text().isEmpty() && QFileInfo(E_WorkDir).exists())
            {
                E_MapDir->setText(E_WorkDir);
                if (!QFileInfo(E_MapDir->text()).exists())
                    E_MapDir->setText(E_WorkDir);
            }
        }
        if (p1.compare("Result Directory")==0)
        {
            if (doBatchmode)
                E_ResultDir->setText(CheckDir(p, true));
            else
                E_ResultDir->setText(CheckDir(p, false));
            if (!QFileInfo(E_ResultDir->text()).exists() && QFileInfo(E_WorkDir).exists())
                E_ResultDir->setText(E_WorkDir + "res/");
        }

        if (p1.compare("Main results file")==0) E_MainTotals->setText(p);
        if (p1.compare("Total Series file")==0) E_SeriesTotals->setText(p);
        if (p1.compare("Filename point output")==0) E_PointResults->setText(p);
   //     if (p1.compare("Filename landunit output")==0) E_LandunitResults->setText(p);
        // resultDir is added in report operation

        if (p1.compare("Rainfall Directory")==0) RainFileDir = CheckDir(p);
        if (p1.compare("Rainfall file")==0)  RainFileName = p;

        if (p1.compare("Rainfall Map Directory")==0) RainSatFileDir = CheckDir(p);
        if (p1.compare("Rainfall maplist name")==0)  RainSatFileName = p;

        if (p1.compare("ET Directory")==0) ETFileDir = CheckDir(p);
        if (p1.compare("ET file")==0) ETFileName = p;

        if (p1.compare("ET Map Directory")==0) ETSatFileDir = CheckDir(p);
        if (p1.compare("ET maplist name")==0) ETSatFileName = p;

        if (p1.compare("Snowmelt Directory")==0) SnowmeltFileDir = CheckDir(p);
        if (p1.compare("Snowmelt file")==0) SnowmeltFileName = p;

//        if (p1.compare("Snowmelt Map Directory")==0) SnowmeltSatFileDir = CheckDir(p);
//        if (p1.compare("Snowmelt maplist name")==0) {
//            SnowmeltSatFileName = p;
//            E_SnowmeltSatName->setText(SnowmeltSatFileDir + SnowmeltSatFileName);
//            if (!QFileInfo(E_SnowmeltSatName->text()).exists())
//            {
//                SnowmeltSatFileDir = QString(E_WorkDir + "rain/");
//                E_SnowmeltSatName->setText(SnowmeltSatFileDir + p);
//            }
//        }

//        if (p1.compare("Discharge inflow directory")==0) DischargeinDir = CheckDir(p);
//        if (p1.compare("Discharge inflow file")==0)
//        {
//            E_DischargeInName->setText(DischargeinDir + p);
//            DischargeinFileName = p;
//            if (!QFileInfo(E_DischargeInName->text()).exists())
//            {
//                DischargeinDir = QString(E_WorkDir + "rain/");
//                E_DischargeInName->setText(DischargeinDir + p);
//            }
//        }

        if (p1.compare("satImage Directory")==0) satImageFileDir = CheckDir(p);
        if (p1.compare("satImage file")==0) satImageFileName = p;

        if (p1.compare("Rainfall map")==0) E_RainfallMap->setText(p);
        if (p1.compare("Interception map")==0) E_InterceptionMap->setText(p);
        if (p1.compare("Infiltration map")==0) E_InfiltrationMap->setText(p);
        if (p1.compare("Runoff map")==0) E_RunoffMap->setText(p);
        if (p1.compare("Channel discharge map")==0) E_ChannelQtotm3Map->setText(p);

        if (p1.compare("Erosion map")==0) E_DetachmentMap->setText(p);
        if (p1.compare("Deposition map")==0) E_DepositionMap->setText(p);
        if (p1.compare("Soilloss map")==0) E_SoillossMap->setText(p);
        if (p1.compare("Channel detachment map")==0) E_ChanDetachmentMap->setText(p);
        if (p1.compare("Channel deposition map")==0) E_ChanDepositionMap->setText(p);

        if (p1.compare("Flood time map")==0) E_FloodTimeMap->setText(p);
        if (p1.compare("Flood start time")==0) E_FloodFEW->setText(p);
        if (p1.compare("Channel Max Q")==0) E_ChannelMaxQ->setText(p);
        if (p1.compare("Channel Max WH")==0) E_ChannelMaxWH->setText(p);
        if (p1.compare("Max Velocity")==0) E_FloodmaxVMap->setText(p);
        if (p1.compare("Max Momentum")==0) E_FloodmaxVHMap->setText(p);
        if (p1.compare("Flood stats")==0) E_FloodStats->setText(p);

        if (p1.compare("Storm Drain map")==0) E_stormDrainMap->setText(p);
        if (p1.compare("Storm Drain Vol map")==0) E_stormDrainVolMap->setText(p);

        if (uiInfilMethod == 1 && p1.compare("Table Directory")==0)
        {
            SwatreTableDir = CheckDir(p);
            if (SwatreTableDir.isEmpty())
                SwatreTableDir = E_MapDir->text();
            E_SwatreTableDir->setText(SwatreTableDir);
        }

        if (uiInfilMethod == 1 && p1.compare("Table File")==0)
        {
            SwatreTableName = p;
            if (SwatreTableName.isEmpty())
                SwatreTableName = E_MapDir->text() + QString("profile.inp");
            E_SwatreTableName->setText(SwatreTableName);
        }
    }

    E_rainsatName->setText(RainSatFileDir + RainSatFileName);
    //qDebug() << RainSatFileName << RainSatFileDir;
    if (!QFileInfo(E_rainsatName->text()).exists())
    {
        RainSatFileDir = QString(E_WorkDir + "rain/");
        E_rainsatName->setText(RainSatFileDir + RainSatFileName);
    }
    E_RainfallName->setText(RainFileDir + RainFileName);
    if (!QFileInfo(E_RainfallName->text()).exists())
    {
        RainFileDir = QString(E_WorkDir + "rain/");
        E_RainfallName->setText(RainFileDir + RainFileName);
    }


    if (checkIncludeET->isChecked()) {
        E_ETName->setText(ETFileDir + ETFileName);
        if (!QFileInfo(E_ETName->text()).exists() && !E_ETName->text().isEmpty())
        {
            ETFileDir = QString(E_WorkDir + "rain/");
            E_ETName->setText(ETFileDir + ETFileName);
        }

        E_ETsatName->setText(ETSatFileDir + ETSatFileName);
        if (!QFileInfo(E_ETsatName->text()).exists() && !E_ETsatName->text().isEmpty())
        {
            ETSatFileDir = QString(E_WorkDir + "rain/");
            E_ETsatName->setText(ETSatFileDir + ETSatFileName);
        }
    }

//    E_SnowmeltName->setText(SnowmeltFileDir + SnowmeltFileName);
//    if (!QFileInfo(E_SnowmeltName->text()).exists())
//    {
//        SnowmeltFileDir = QString(E_WorkDir + "rain/");
//        E_SnowmeltName->setText(ETFileDir + SnowmeltFileName);
//    }

    E_satImageName->setText(satImageFileDir +satImageFileName);
    if (!QFileInfo(E_satImageName->text()).exists())
    {
        satImageFileDir = "";//QString(E_WorkDir + "maps/");
        satImageFileName = "";
        E_satImageName->setText("");
        //E_satImageName->setText(satImageFileDir + satImageFileName);
    }

    int days = daystart.toInt();
    int mins = minstart.toInt();
    int daye = dayend.toInt();
    int mine = minend.toInt();

    days = std::max(1,std::min(days, 366));
    daye = std::max(1,std::min(daye, 366));
//qDebug() << days << mins << daye << mine;
    if (!checkEventBased->isChecked()) {
        if (mins > 1440) {
           days = mins/1440 + 1;
           mins = mins % 1440;
        }
        if (mine > 1440) {
            daye = mine/1440 + 1;
            mine = mine % 1440;
        }
    }

    E_BeginTimeDay->setText(QString("%1:%2").arg(days,3, 10, QLatin1Char('0')).arg(mins,4, 10, QLatin1Char('0')));
    E_EndTimeDay->setText(QString("%1:%2").arg(daye,3, 10, QLatin1Char('0')).arg(mine,4, 10, QLatin1Char('0')));

    on_checkIncludeChannel_clicked();
    on_checkMaterialDepth_clicked();

    //****====------====****//

    // get all map names, DEFmaps contains default map names and descriptions
    // adapt the DEFmaps list with names from the run file
    // this is to display the correct names in the interface
    for (j = mapstartnr; j < nrnamelist; j++)
    {
        for (int i = 0; i < DEFmaps.size(); i++)
        {
            QStringList S = DEFmaps.at(i).split(";",Qt::SkipEmptyParts);
            if (S.contains(namelist[j].name))
            {
                if(namelist[j].value == "chansedmixdeth.map") namelist[j].value = "chansedmixdepth.map";
                if(namelist[j].value == "sedmixdeth.map") namelist[j].value = "sedmixdepth.map";
                //qDebug() << j << namelist[j].value;

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

//for (int j = 0; j < nrnamelist; j++)
//    qDebug() << namelist[j].name << namelist[j].value;
    // fill the mapList structure with all map names fom the runfile
    // if there are new variables that are not in the run file
    // the maplist contains the default names already
    // this is to get the correct names for the model run
    fillNamelistMapnames(false);
    for (int k = 0; k < nrmaplist; k++)
        mapList[k].dir = E_MapDir->text();

    //RunAllChecks();
    //obsolete: is done in void lisemqt::on_E_runFileList_currentIndexChanged(int)

}
//---------------------------------------------------------------------------
QString lisemqt::CheckDir(QString p, bool makeit)
{
    /* TODO mulitplatform: fromNativeSeparators etc*/
    QString path;

    if (p.isEmpty() || p == "/")
        return(p);

    path = QDir(p).fromNativeSeparators(p);
    path = QDir(path).absoluteFilePath(path);
    if (!path.endsWith("/"))
        path = path + '/';

    if (!QDir(path).exists())
    {
        if (makeit)
        {
            QDir(path).mkpath(path);
            //qDebug() << "NOTE: Result dir created !";
        }
        else
        {
            QMessageBox::warning(this,"openLISEM",QString("The following directory does not exist:\n%1\nUsing the work directory, check your pathnames").arg(path));
            path.clear();
        }
    }

    return path;
}
//---------------------------------------------------------------------------
//! This function changes the runfile with the currentinterface options
//! It is is called by save file, befor ethe actual saving
//! It is ALSO called just before the model is run to create a tmp runfile for the model to read
void lisemqt::updateModelData()
{
//    if(!checkRainfall->isChecked() && !checkSnowmelt->isChecked())
//        QMessageBox::warning(this,"openLISEM","No rainfall or snowmelt, running on empty!");

    int days = E_BeginTimeDay->text().split(":")[0].toInt();
    int mins = E_BeginTimeDay->text().split(":")[1].toInt();
    int daye = E_EndTimeDay->text().split(":")[0].toInt();
    int mine = E_EndTimeDay->text().split(":")[1].toInt();
   // qDebug() << days << mins << daye << mine;

    if (!checkEventBased->isChecked()) {
    if (mins > 1440) {
       days = mins/1440 + 1;
       mins = mins % 1440;
    }
    if (mine > 1440) {
        daye = mine/1440 + 1;
        mine = mine % 1440;
    }
    }
 //   qDebug() << days << mins << daye << mine;

    for (int j = 0; j < nrnamelist; j++)
    {
        QString p1 = namelist[j].name;

        if (p1.compare("Event based")==0)              namelist[j].value.setNum((int)checkEventBased->isChecked());
        if (p1.compare("Use Rainfall maps")==0)        namelist[j].value.setNum((int)radioRainSatFile->isChecked());
        if (p1.compare("Rainfall ID interpolation")==0)namelist[j].value.setNum((int)checkIDinterpolation->isChecked());
        if (p1.compare("Include ET")==0)               namelist[j].value.setNum((int)checkIncludeET->isChecked());
        if (p1.compare("Use ET maps")==0)              namelist[j].value.setNum((int)radioETsatfile->isChecked());
        if (p1.compare("Daily ET")==0)                 namelist[j].value.setNum((int)checkDailyET->isChecked());

        //if (p1.compare("Include Snowmelt")==0)               namelist[j].value.setNum((int)checkSnowmelt->isChecked());

        if (p1.compare("Nr user Cores")==0) namelist[j].value.setNum(nrUserCores->value());
        // erosion
        if (p1.compare("Include Erosion simulation")==0)      namelist[j].value.setNum((int)checkDoErosion->isChecked());

        //channels
        if (p1.compare("Include main channels")==0)          namelist[j].value.setNum((int)checkIncludeChannel->isChecked());
        if (p1.compare("Include channel infil")==0)          namelist[j].value.setNum((int)checkChannelInfil->isChecked());
        if (p1.compare("Include channel baseflow")==0)       namelist[j].value.setNum((int)checkChannelBaseflow->isChecked());
        if (p1.compare("Include stationary baseflow")==0)    namelist[j].value.setNum((int)checkStationaryBaseflow->isChecked());
      //  if (p1.compare("Adjust channel crosssection")==0)    namelist[j].value.setNum((int)checkChannelAdjustCHW->isChecked());
        if (p1.compare("Include channel culverts")==0)       namelist[j].value.setNum((int)checkChannelCulverts->isChecked());
        if (p1.compare("Include channel inflow")==0)         namelist[j].value.setNum((int)checkChannelInflow->isChecked());
        // groundwater
        if (p1.compare("GW flow explicit")==0)               namelist[j].value.setNum((int)checkGWflowexplicit->isChecked());
        if (p1.compare("GW flow SWAT")==0)               namelist[j].value.setNum((int)checkGWflowSWAT->isChecked());
        if (p1.compare("GW recharge factor")==0)             namelist[j].value = GW_recharge->text();
        if (p1.compare("GW flow factor")==0)                 namelist[j].value = GW_flow->text();
        if (p1.compare("GW river inflow factor")==0)                 namelist[j].value = GW_inflow->text();
        if (p1.compare("GW slope factor")==0)                namelist[j].value = GW_slope->text();
   //     if (p1.compare("GW lag factor")==0)                  namelist[j].value = GW_lag->text();
        if (p1.compare("GW deep percolation")==0)            namelist[j].value = GW_deep->text();
        if (p1.compare("GW threshold factor")==0)            namelist[j].value = GW_threshold->text();
   //     if (p1.compare("GW initial level")==0)               namelist[j].value = GW_initlevel->text();

        if (p1.compare("Include flow barriers")==0)          namelist[j].value.setNum((int)checkFlowBarriers->isChecked());
        if (p1.compare("Include buffers")==0)                namelist[j].value.setNum((int) checkBuffers->isChecked());
        if (p1.compare("Flow barrier table filename")==0)    namelist[j].value = line_FlowBarriers->text();

        if (p1.compare("Include litter interception")==0)    namelist[j].value.setNum((int)checkIncludeLitter->isChecked());
        if (p1.compare("Litter interception storage")==0)    namelist[j].value = E_LitterSmax->text();

        if (p1.compare("Routing Kin Wave 2D")==0)
        {
            if (checkOverlandFlow1D->isChecked())  namelist[j].value = "1";
            if (checkOverlandFlow2Ddyn->isChecked())  namelist[j].value = "3";
            if (checkOverlandFlow2Dkindyn->isChecked())  namelist[j].value = "4";
        }
        if (p1.compare("Flow Boundary 2D")==0)               namelist[j].value = E_FlowBoundary->text();
        if (p1.compare("Flooding courant factor")==0)        namelist[j].value = E_courantFactor->text();
        if (p1.compare("Use gravity flow")==0)               namelist[j].value.setNum((int)checkGravityToChannel->isChecked());
      //  if (p1.compare("Angle flow to channel")==0)          namelist[j].value = E_angleToChannel->text();
        if (p1.compare("Include diffusion")==0)              namelist[j].value.setNum((int)checkDiffusion->isChecked());
        if (p1.compare("Sigma diffusion")==0)                namelist[j].value = E_SigmaDiffusion->text();
        if (p1.compare("Include River diffusion")==0)              namelist[j].value.setNum((int)checkDiffusion->isChecked());
        if (p1.compare("Flooding SWOF flux limiter")==0)     namelist[j].value = E_FloodFluxLimiter->text();
        if (p1.compare("Flooding SWOF Reconstruction")==0)   namelist[j].value = E_FloodReconstruction->text();
        if (p1.compare("Minimum reported flood height")==0)  namelist[j].value = E_floodMinHeight->text();
        if (p1.compare("Flooding mixing coefficient")==0)    namelist[j].value = E_mixingFactor->text();
        if (p1.compare("Flooding runoff partitioning")==0)   namelist[j].value = E_runoffPartitioning->text();
        if (p1.compare("Flood initial level map")==0)        namelist[j].value.setNum((int)checkFloodInitial->isChecked());
        if (p1.compare("Pit Value")==0)                      namelist[j].value = E_pitValue->text();
   //     if (p1.compare("Calculate erosion inside 2D loop")==0) namelist[j].value.setNum((int)checkErosionInsideLoop->isChecked());
        if (p1.compare("Use linked list")==0)                namelist[j].value.setNum((int)checkLinkedList->isChecked());
        if (p1.compare("Use Channel Kinwave dt")==0)         namelist[j].value.setNum((int)checkKinWaveChannel->isChecked());
        if (p1.compare("Channel KinWave dt")==0)             namelist[j].value = E_ChannelKinWaveDt->text();
        if (p1.compare("Use Channel Max V")==0)              namelist[j].value.setNum((int)checkChanMaxVelocity->isChecked());
        if (p1.compare("Channel Max V")==0)                  namelist[j].value = E_chanMaxVelocity->text();

        if (p1.compare("Flood max iterations")==0)           namelist[j].value = E_FloodMaxIter->text();
        if (p1.compare("Min WH flow")==0)                    namelist[j].value = E_minWHflow->text();

        if (p1.compare("Timestep flood")==0)                 namelist[j].value = E_TimestepMinFlood->text();
        if (p1.compare("Use time avg V")==0)                 namelist[j].value.setNum((int) checkTimeavgV->isChecked());
        if (p1.compare("Correct DEM")==0)                    namelist[j].value.setNum((int) checkCorrectDem->isChecked());
        if (p1.compare("Use 2D Diagonal flow")==0)           namelist[j].value.setNum((int) check2DDiagonalFlow->isChecked());
        if (p1.compare("Use 2D Diagonal flow new")==0)       namelist[j].value.setNum((int) check2DDiagonalFlowNew->isChecked());

    //    if (p1.compare("Use SWOF watersheds")==0)            namelist[j].value.setNum((int) checkSWOFWatersheds->isChecked());

    //    if (p1.compare("Use fixed angle")==0)                namelist[j].value.setNum((int) checkFixedAngle->isChecked());
     //   if (p1.compare("Variable Timestep")==0)              namelist[j].value.setNum((int) checkVariableTimestep->isChecked());

        if (p1.compare("Flood solution")==0)
        {
            int i = 0;
            if(rb_swof2->isChecked()) i = 0;
            if(rb_swofnomuscl->isChecked()) i = 1;
            if(rb_swofmuscl->isChecked()) i = 2;
            namelist[j].value.setNum(i);
        }

        if (p1.compare("Splash equation")==0)                namelist[j].value = E_splashEquation->text();
        if (p1.compare("Detachment efficiency")==0)          namelist[j].value = QString::number(E_EfficiencyDET->currentIndex()+1);
        if (p1.compare("Detachment efficiency channel")==0)  namelist[j].value = QString::number(E_EfficiencyDETCH->currentIndex()+1);
        if (p1.compare("Direct efficiency channel")==0)      namelist[j].value = E_EfficiencyDirect->text();
        if (p1.compare("Settling Velocity")==0)              namelist[j].value = QString::number(E_settlingVelocity->currentIndex()+1);

        if (p1.compare("Use 2 phase flow")==0)              namelist[j].value.setNum((int) checkSed2Phase->isChecked());

        if (p1.compare("Use material depth")==0)             namelist[j].value.setNum((int)checkMaterialDepth->isChecked());
        if (p1.compare("No detachment boundary")==0)         namelist[j].value.setNum((int)checkNoSedBoundary->isChecked());

        if (p1.compare("Flooding BL method")==0)             namelist[j].value = QString::number(E_BLMethod->currentIndex()+1);
        if (p1.compare("Flooding SS method")==0)             namelist[j].value = QString::number(E_SSMethod->currentIndex()+1);
        if (p1.compare("River BL method")==0)                namelist[j].value = QString::number(E_RBLMethod->currentIndex()+1);
        if (p1.compare("River SS method")==0)                namelist[j].value = QString::number(E_RSSMethod->currentIndex()+1);

//        if (p1.compare("Use grain size distribution")==0)    namelist[j].value.setNum((int)checkSedMultiGrain->isChecked());
//        if (p1.compare("Estimate grain size distribution")==0)namelist[j].value.setNum((int)checkEstimateGrainSizeDistribution->isChecked());
//        if (p1.compare("Read grain distribution maps")==0)    namelist[j].value.setNum((int)checkReadGrainSizeDistribution->isChecked());
//        if (p1.compare("Number of grain size classes (simulated)")==0)  namelist[j].value = E_NumberClasses->text();
//        if (p1.compare("Grain size class maps")==0)   namelist[j].value = E_GrainSizes->text();

        //tile drains
        if (p1.compare("Include tile drains")==0)            namelist[j].value.setNum((int)checkIncludeTiledrains->isChecked());

        //houses
        if (p1.compare("Include house storage")==0)          namelist[j].value.setNum((int)checkHouses->isChecked());
        if (p1.compare("Include buildings")==0)          namelist[j].value.setNum((int)checkHouses->isChecked());
        if (p1.compare("Add buildings to DEM")==0)           namelist[j].value.setNum((int)checkAddBuildingDEM->isChecked());
        if (p1.compare("Add building fraction")==0)           namelist[j].value = E_AddBuildingFraction->text();
        if (p1.compare("Include raindrum storage")==0)       namelist[j].value.setNum((int)checkRaindrum->isChecked());
        if (p1.compare("Include road system")==0)            namelist[j].value.setNum((int)checkRoadsystem->isChecked());
        if (p1.compare("Hard Surfaces")==0)                  namelist[j].value.setNum((int)checkHardsurface->isChecked());
        if (p1.compare("Include storm drains")==0)           namelist[j].value.setNum((int)checkStormDrains->isChecked());

        if (p1.compare("Include Satellite Image")==0)        namelist[j].value.setNum((int)checksatImage->isChecked());

        if (p1.compare("Include Sediment traps")==0)         namelist[j].value.setNum((int)checkSedtrap->isChecked());
        if (p1.compare("Include compacted")==0)            namelist[j].value.setNum((int)checkInfilCompact->isChecked());
        if (p1.compare("Include grass strips")==0)           namelist[j].value.setNum((int)checkInfilGrass->isChecked());
        if (p1.compare("Grassstrip Mannings n")==0)          namelist[j].value = E_GrassStripN->text();
        if (p1.compare("Sediment Trap Mannings n")==0)          namelist[j].value = E_SedTrapN->text();

        if (p1.compare("Include crusts")==0)                 namelist[j].value.setNum((int)checkInfilCrust->isChecked());
        if (p1.compare("Impermeable sublayer")==0)          namelist[j].value.setNum((int)checkInfilImpermeable->isChecked());
        if (p1.compare("Two layer")==0)                     namelist[j].value.setNum((int)checkInfil2layer->isChecked());
        //if (p1.compare("Matric head files")==0)              namelist[j].value.setNum((int)checkDumphead->isChecked());
        if (p1.compare("Geometric mean Ksat")==0)            namelist[j].value.setNum((int)checkGeometric->isChecked());
        if (p1.compare("Timeplot as PCRaster")==0)           namelist[j].value.setNum(checkWritePCRaster->isChecked() ? 0 : 1);
        if (p1.compare("Report point output separate")==0)   namelist[j].value.setNum((int)checkSeparateOutput->isChecked());
        if (p1.compare("Report digits out")==0)             namelist[j].value = E_DigitsOut->text();

        if (p1.compare("Report format GTiff")==0)             namelist[j].value.setNum((int)checkFormatGtiff->isChecked());
        if (p1.compare("End run report")==0)                namelist[j].value.setNum((int)checkEndRunReport->isChecked());

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
            if(radioButton_9->isChecked()) i = 8;  // 8 = storage map, not LAI eq
            namelist[j].value.setNum(i);
        }

        //VJ 110705 KE equations
        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
            param << (radioButtonKE1->isChecked()?"1":"0") << spinKEparameterA1->text() << spinKEparameterB1->text() << spinKEparameterC1->text();
            namelist[j].value = param.join(";");
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
            param << (radioButtonKE2->isChecked()?"1":"0") << spinKEparameterA2->text() << spinKEparameterB2->text();
            namelist[j].value = param.join(";");
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
            param << (radioButtonKE3->isChecked()?"1":"0") << spinKEparameterA3->text() << spinKEparameterB3->text();
            namelist[j].value = param.join(";");
        }
        if (p1.compare("KE time based")==0)      namelist[j].value.setNum((int)checkKETimebased->isChecked());

        if (p1.compare("Begin time day")==0) namelist[j].value = QString("%1").arg(days,3,10, QLatin1Char('0'));//E_BeginTimeDay->text();
        if (p1.compare("Begin time")==0) namelist[j].value = QString("%1").arg(mins,4,10, QLatin1Char('0'));//E_BeginTimeMin->text();
        if (p1.compare("End time day")==0)   namelist[j].value = QString("%1").arg(daye,3,10, QLatin1Char('0'));//E_EndTimeDay->text();
        if (p1.compare("End time")==0)   namelist[j].value = QString("%1").arg(mine,4,10, QLatin1Char('0'));//E_EndTimeMin->text();
        if (p1.compare("Timestep")==0)   namelist[j].value = E_Timestep->text();
     //   if (p1.compare("Work Directory")==0)    namelist[j].value = E_WorkDir;//->text();
        if (p1.compare("Map Directory")==0)    namelist[j].value = E_MapDir->text();
        if (p1.compare("Result Directory")==0) namelist[j].value = E_ResultDir->text();
        if (p1.compare("Main results file")==0) namelist[j].value = E_MainTotals->text();
        if (p1.compare("Total series file")==0) namelist[j].value = E_SeriesTotals->text();
        if (p1.compare("Filename point output")==0) namelist[j].value = E_PointResults->text();
     //   if (p1.compare("Filename landunit output")==0) namelist[j].value = E_LandunitResults->text();

        if (p1.compare("Rainfall file")==0) namelist[j].value = RainFileName;
        if (p1.compare("Rainfall Directory")==0) namelist[j].value = RainFileDir;
        if (p1.compare("Rainfall maplist name")==0) namelist[j].value = RainSatFileName;
        if (p1.compare("Rainfall Map Directory")==0) namelist[j].value = RainSatFileDir;
        if (p1.compare("Rainfall Bias Correction")==0) namelist[j].value = E_biasCorrectionP->text();
        if (p1.compare("IDI factor")==0)            namelist[j].value = E_IDIfactor->text();

        if (p1.compare("ET file") ==0)        namelist[j].value = ETFileName;
        if (p1.compare("ET Directory") ==0)   namelist[j].value = ETFileDir;
        if (p1.compare("ET maplist name") ==0) namelist[j].value = ETSatFileName;
        if (p1.compare("ET Map Directory") ==0) namelist[j].value = ETSatFileDir;
        if (p1.compare("ET Bias Correction")==0) namelist[j].value = E_biasCorrectionET->text();
        if (p1.compare("Rainfall ET threshold")==0) namelist[j].value = E_rainfallETA_threshold->text();

      //  if (p1.compare("Snowmelt Directory")==0) namelist[j].value = SnowmeltFileDir;
      //  if (p1.compare("Snowmelt file")==0) namelist[j].value = SnowmeltFileName;

        if (p1.compare("satImage Directory")==0) namelist[j].value = satImageFileDir;
        if (p1.compare("satImage file")==0) namelist[j].value = satImageFileName;
        if (p1.compare("Advanced Options")==0) namelist[j].value.setNum((int)checkAdvancedOptions->isChecked());

        if (p1.compare("Rainfall map")==0) namelist[j].value = E_RainfallMap->text();
        if (p1.compare("Interception map")==0) namelist[j].value = E_InterceptionMap->text();
        if (p1.compare("Infiltration map")==0) namelist[j].value = E_InfiltrationMap->text();
        if (p1.compare("Runoff map")==0) namelist[j].value = E_RunoffMap->text();
        // if (p1.compare("Runoff fraction map")==0) namelist[j].value = E_RunoffFractionMap->text();
        if (p1.compare("Channel discharge map")==0) namelist[j].value = E_ChannelQtotm3Map->text();
        if (p1.compare("WH max level map")==0) namelist[j].value = E_WHmaxMap->text();

        //if (p1.compare("Flood level map")==0) namelist[j].value = E_FloodlevelMap->text();
        if (p1.compare("Flood time map")==0) namelist[j].value = E_FloodTimeMap->text();
        if (p1.compare("Flood start time")==0) namelist[j].value = E_FloodFEW->text();
        if (p1.compare("Max Velocity")==0) namelist[j].value = E_FloodmaxVMap->text();
        if (p1.compare("Max Momentum")==0) namelist[j].value = E_FloodmaxVHMap->text();
        if (p1.compare("Channel Max Q")==0) namelist[j].value = E_ChannelMaxQ->text();
        if (p1.compare("Channel Max WH")==0) namelist[j].value = E_ChannelMaxWH->text();
        if (p1.compare("Flood stats")==0) namelist[j].value = E_FloodStats->text();
        if (p1.compare("Storm Drain map")==0) namelist[j].value = E_stormDrainMap->text();
        if (p1.compare("Storm Drain Vol map")==0) namelist[j].value = E_stormDrainVolMap->text();

        if (p1.compare("Erosion map")==0) namelist[j].value = E_DetachmentMap->text();
        if (p1.compare("Deposition map")==0) namelist[j].value = E_DepositionMap->text();
        if (p1.compare("Soilloss map")==0) namelist[j].value = E_SoillossMap->text();
        if (p1.compare("Channel detachment map")==0) namelist[j].value = E_ChanDetachmentMap->text();
        if (p1.compare("Channel deposition map")==0) namelist[j].value = E_ChanDepositionMap->text();

        if (p1.compare("Grain Size calibration D50")==0)   namelist[j].value = E_CalibrateD50->text();
        if (p1.compare("Grain Size calibration D90")==0)   namelist[j].value = E_CalibrateD90->text();
        if (p1.compare("Smax calibration")==0) namelist[j].value = E_CalibrateSmax->text();
        if (p1.compare("Ksat calibration")==0) namelist[j].value = E_CalibrateKsat->text();
        if (p1.compare("Ksat2 calibration")==0) namelist[j].value = E_CalibrateKsat2->text();
        if (p1.compare("N calibration")==0) namelist[j].value = E_CalibrateN->text();
        if (p1.compare("Theta calibration")==0) namelist[j].value = E_CalibrateTheta->text();
        if (p1.compare("Psi calibration")==0) namelist[j].value = E_CalibratePsi->text();
        if (p1.compare("Channel Ksat calibration")==0) namelist[j].value = E_CalibrateChKsat->text();
        if (p1.compare("Channel N calibration")==0) namelist[j].value = E_CalibrateChN->text();
        if (p1.compare("Channel tortuosity")==0) namelist[j].value = E_CalibrateChTor->text();
        if (p1.compare("Cohesion calibration")==0) namelist[j].value = E_CalibrateCOH->text();
        if (p1.compare("Cohesion Channel calibration")==0) namelist[j].value = E_CalibrateCHCOH->text();
        if (p1.compare("Ucr Channel calibration")==0) namelist[j].value = E_CalibrateCHUcr->text();
        if (p1.compare("SV calibration")==0) namelist[j].value = E_CalibrateCHSV->text();
        if (p1.compare("Aggregate stability calibration")==0) namelist[j].value = E_CalibrateAS->text();
        if (p1.compare("Splash Delivery Ratio")==0) namelist[j].value = E_SplashDelibery->text();
        if (p1.compare("Particle Cohesion of Deposited Layer")==0) namelist[j].value = E_DepositedCohesion->text();
     //   if (p1.compare("Stemflow fraction")==0) namelist[j].value = E_StemflowFraction->text();
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
            double fraction = 0.2;//E_SWATREDtsecFraction->value();
            swatreDT = E_Timestep->text().toDouble()*fraction;
            namelist[j].value.setNum(swatreDT,'g',6);
        }

        if (p1.compare("Infil Method")==0)
        {
            switch(E_InfiltrationMethod->currentIndex())
            {
            case 0 : namelist[j].value.setNum(INFIL_NONE);break;
            case 1 : namelist[j].value.setNum(INFIL_SWATRE);break;
            case 2 : namelist[j].value.setNum(INFIL_GREENAMPT);break;
            case 3 : namelist[j].value.setNum(INFIL_SMITH); break;
            }
        }

        //VJ 111120 water repellency
//        if (p1.compare("Use Water Repellency")==0)      namelist[j].value.setNum((int)checkWaterRepellency->isChecked());
//        if (p1.compare("Water Repellency A")==0)        namelist[j].value.setNum(E_waterRep_a->value(),'g',6);
//        if (p1.compare("Water Repellency B")==0)        namelist[j].value.setNum(E_waterRep_b->value(),'g',6);
//        if (p1.compare("Water Repellency C")==0)        namelist[j].value.setNum(E_waterRep_c->value(),'g',6);
//        if (p1.compare("Water Repellency D")==0)        namelist[j].value.setNum(E_waterRep_d->value(),'g',6);

        if (p1.compare("Report discharge units")==0)
        {
            if (checkUnits_ls->isChecked())
                namelist[j].value.setNum(0);
            if (checkUnits_m3s->isChecked())
                namelist[j].value.setNum(1);
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

        // make a string for all output maps


        if (p1.compare("OutRunoff")==0)         namelist[j].value.setNum((int)checkBox_OutRunoff->isChecked());
        if (p1.compare("OutWH")==0)             namelist[j].value.setNum((int)checkBox_OutWH->isChecked());
        if (p1.compare("OutV")==0)              namelist[j].value.setNum((int)checkBox_OutV->isChecked());
        if (p1.compare("OutInterception")==0)  namelist[j].value.setNum((int)checkBox_OutInterception->isChecked());
        if (p1.compare("OutSurfStor")==0)       namelist[j].value.setNum((int)checkBox_OutSurfStor->isChecked());
        if (p1.compare("OutInf")==0)            namelist[j].value.setNum((int)checkBox_OutInf->isChecked());
        if (p1.compare("OutTileDrain")==0)      namelist[j].value.setNum((int)checkBox_OutTiledrain->isChecked());
        if (p1.compare("OutTileVolume")==0)     namelist[j].value.setNum((int)checkBox_OutTileVol->isChecked());
        if (p1.compare("OutTheta")==0)         namelist[j].value.setNum((int)checkBox_OutTheta->isChecked());

        if (p1.compare("OutDet")==0)     namelist[j].value.setNum((int)checkBox_OutDet->isChecked());
        if (p1.compare("OutDep")==0)     namelist[j].value.setNum((int)checkBox_OutDep->isChecked());
        if (p1.compare("OutTC")==0)      namelist[j].value.setNum((int)checkBox_OutTC->isChecked());
        if (p1.compare("OutConc")==0)    namelist[j].value.setNum((int)checkBox_OutConc->isChecked());
        if (p1.compare("OutSed")==0)     namelist[j].value.setNum((int)checkBox_OutSed->isChecked());
        if (p1.compare("OutSL")==0)      namelist[j].value.setNum((int)checkBox_OutSL->isChecked());
        if (p1.compare("OutSedSS")==0)     namelist[j].value.setNum((int)checkBox_OutSedSS->isChecked());
        if (p1.compare("OutSedBL")==0)     namelist[j].value.setNum((int)checkBox_OutSedBL->isChecked());

        if (p1.compare("Result datetime")==0) namelist[j].value.setNum((int)checkAddDatetime->isChecked());
     //   if (p1.compare("Add timestamp")==0)   namelist[j].value.setNum((int)checkOutputTimestamp->isChecked());

    }
    //get all actual mapnames from the mapList structure
    fillNamelistMapnames(true);

    currentDir = E_WorkDir;
    QDir::setCurrent(currentDir);

    if (saveRunFileOnce) {
        savefile(op.runfilename);
        saveRunFileOnce = false;
        QMessageBox::warning(this,"openLISEM",QString("The run file has changed: ") +
            QString("obsolete options are removed and missing options use default values. ") +
            QString("The new run files has your choices where applicable."));
    }

}
