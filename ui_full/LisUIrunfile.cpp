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

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
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
        }

    }

    if (optionList.count() > 0)
    {
        for (i=0; i< optionList.count(); i++)
        {


            if (optionList[i].contains("="))
            {
                QString S = optionList[i];
                S.remove('[');
                S.remove(']');
                QStringList SL = S.split(QRegExp("="));

                for (int j = 0; j < nrnamelist; j++)
                {
                    if (namelist[j].name == SL[0].trimmed())
                    {
                        namelist[j].value = SL[1].trimmed();
                        break;
                    }
                }
            }
        }

    }
}
//---------------------------------------------------------------------------
//! ParseInputData : interpret runfile text and fill interface variables
void lisemqt::ParseInputData()
{
    int j=0;
    bool dummyrain = false;
    bool dummysnow = false;
    bool dummy2layerinfil = false;
    bool dummyErosion = false;

    bool seterosionold = false;
    // get all the options/checks
    for (j = 0; j < nrnamelist; j++)  //VJ 110107 changed to nrnamelist
    {
        int iii = namelist[j].value.toInt();
        double val = namelist[j].value.toDouble();
        QString p1 = namelist[j].name;
        QString p = namelist[j].value;
        bool check = iii == 1;
        if (p1.contains("["))
            continue;

        //options in the main code, order is not important
        if (p1.compare("No Erosion simulation")==0){

            checkDoErosion->setChecked(!check);
            dummyErosion = !check;
            if(!check)
            {
                seterosionold = true;
            }

        }
        if (p1.compare("Include Erosion simulation")==0){

            if(!seterosionold)
            {
                dummyErosion = check || seterosionold;
                checkDoErosion->setChecked(check || seterosionold);
            }
        }
        if (p1.compare("Include main channels")==0)          checkIncludeChannel->setChecked(check);
        if (p1.compare("Include channel infil")==0)          checkChannelInfil->setChecked(check);
        if (p1.compare("Include channel baseflow")==0)       checkChannelBaseflow->setChecked(check);

    //    if (p1.compare("Include channel flooding")==0)       checkChannelFlood->setChecked(check);
        if (p1.compare("Include road system")==0)            checkRoadsystem->setChecked(check);

        //houses
        if (p1.compare("Include house storage")==0)          checkHouses->setChecked(check);
        if (p1.compare("Include raindrum storage")==0)       checkRaindrum->setChecked(check);
        // flooding
       // if (p1.compare("Flood method explicit")==0)        dummyFloodExplicit = check;

        if (p1.compare("Include litter interception")==0)    checkIncludeLitter->setChecked(check);
        if (p1.compare("Litter interception storage")==0)    E_LitterSmax->setValue(val);

        //Surface Flow
        if (p1.compare("Enable Solid Phase")==0)                                    checkSolidPhase->setChecked(check);
        if (p1.compare("Enable Entrainment")==0)                                    checkEntrainment->setChecked(check);
        if (p1.compare("Enable Deposition")==0)                                     checkBox_Deposition->setChecked(check);
        if (p1.compare("Enable Compaction")==0)                                     checkBox_Compaction->setChecked(check);
        if (p1.compare("Enable Seismic") == 0)                                      checkSeismic->setChecked(check);
        if (p1.compare("Enable Slope Stability")==0)                                checkSlopeStability->setChecked(check);
        if (p1.compare("Enable Upslope Forcing") == 0)                              checkUpslopeForcing->setChecked(check);
        if (p1.compare("Enable Downslope Forcing") == 0)                              checkDownslopeForcing->setChecked(check);
        if (p1.compare("Include Bedrock Layer") == 0)                              checkBedRockLayer->setChecked(check);
        if (p1.compare("Enable Slope Failure")==0)                                  checkSlopeFailure->setChecked(check);
        if (p1.compare("Minimum Safety Factor")==0)                                 E_PRESF->setValue(val);
        if (p1.compare("Maximum Safety Factor")==0)                                 E_POSTSF->setValue(val);
        if (p1.compare("Maximum safety factor for display")==0)                     E_SFMaxDisplay->setValue(val);
        if (p1.compare("Entrainment Coefficient")==0)                               E_EntrainmentCoefficient->setValue(val);
        if (p1.compare("Minimum Entrainment Height")==0)                            E_EntrainmentMinimumWaterHeight->setValue(val);
        if (p1.compare("Minimum Failure Height")==0)                                E_MinimumFailureHeight->setValue(val);
        //if (p1.compare("Enable Levees")==0)                                         checkLevees->setChecked(check);
        if (p1.compare("Enable Barriers")==0)                                       checkBox_Barriers->setChecked(check);
        if (p1.compare("Enable Flow Barriers")==0)                                  checkBox_FlowBarriers->setChecked(check);
        if (p1.compare("Flow barrier table filename")==0)                           lineEdit_FlowBarriers->setText(p);
        if (p1.compare("Enable Inflow")==0)                                           checkBox_Inflow->setChecked(check);
        if (p1.compare("Inflow table filename")==0)                           lineEdit_Inflow->setText(p);
        if (p1.compare("Include Initial FluidSolid Mixture")==0)                    checkBox_UFInitial->setChecked(check);
        if (p1.compare("Include Forced FluidSolid Mixture")==0)                     checkBox_UFForced->setChecked(check);
        if (p1.compare("Incldue Maximum ChannelVolume")==0)                         checkBox_ChannelMaxVolume->setChecked(check);
        if (p1.compare("Incldue Maximum Volume")==0)                                checkBox_MaxVolume->setChecked(check);


        if (p1.compare("Suspended Viscosity")==0)                                    checkBox_SuspendedViscosity->setChecked(check);
        if (p1.compare("Lax Multiplier")==0)                                        E_LaxMultiplier->setValue(val);
        if (p1.compare("Friction force correction")==0)                             E_FrictionCorrection->setValue(val);
        if (p1.compare("Erosion Cohesion Calibration")==0)                          E_CalibrateESC->setValue(val);
        if (p1.compare("Erosion Grain Size Calibration")==0)                        E_CalibrateEGS->setValue(val);

        if (p1.compare("Limit Cores")==0)                                           checkBox_CoreLimit->setChecked(check);
        if (p1.compare("Core Limit")==0)                                            E_CoreLimit->setValue(val);

        if (p1.compare("Solid-Fluid Drag Coefficient")==0)                          E_UFDrag->setValue(val);
        if (p1.compare("Use HLL2")==0)                                              checkBox_HLL2->setChecked(check);
        if (p1.compare("Kinematic Timestep Power")==0)                              E_KinematicTimestepPower->setValue(val);

        if (p1.compare("Spatially Dynamic Timestep")==0)                            checkBox_SpatiallyDynamicTimestep->setChecked(check);
        if (p1.compare("Flow Minimum Timestep")==0)                                 E_SFMinimumDT->setValue(val);
        if (p1.compare("Surface Flow Courant Factor")==0)                           E_SFCourant->setValue(val);

        if (p1.compare("Drag Power Law Coefficient")==0)                            E_UFDragPower->setValue(val);
        if (p1.compare("Viscosity Alpha")==0)                                       E_UFViscosityAlpha->setValue(val);
        if (p1.compare("Viscosity Beta")==0)                                        E_UFViscosityBeta->setValue(val);

        if (p1.compare("Minimal Flood Water Depth")==0)                             E_FloodMinDepth->setValue(val);
        if (p1.compare("Minimum Debris Flow Volumetric Sediment Fraction")==0)      E_DebrisFlowMinSedConcentration->setValue(val);

        if (p1.compare("BL method")==0)                     E_BLMethod->setValue(val);
        if (p1.compare("SS method")==0)                     E_SSMethod->setValue(val);
        if (p1.compare("Sigma diffusion")==0)               E_SigmaDiffusion->setValue(val);

        //if (p1.compare("Include levees")==0)                 checkLevees->setChecked(check);

    //    if (p1.compare("OF method")==0)                       E_OFMethod->setValue(val);
        if (p1.compare("Advanced sediment")==0)               checkAdvancedSediment->setChecked(check);
        if (p1.compare("Advanced sediment configuration")==0)
        {
            checkBox_SedSingleSingle->setChecked(false);
            checkBox_SedMultiSingle->setChecked(false);
            checkBox_SedMultiMulti->setChecked(false);
            //tabWidgetOptions->setTabEnabled(6,true);
            if(val == 2)
                checkBox_SedMultiMulti->setChecked(true);
            else
                if(val == 1)
                    checkBox_SedMultiSingle->setChecked(true);
                else
                {
                    checkBox_SedSingleSingle->setChecked(true);
                    //tabWidgetOptions->setTabEnabled(6,false);
                }
        }
        if (p1.compare("Detachment efficiency")==0)          E_EfficiencyDET->setValue(val);
      //  if (p1.compare("Detachment stoniness")==0)           checkStoninessDET->setChecked(check);

        if (p1.compare("Estimate grain size distribution")==0)checkEstimateGrainSizeDistribution->setChecked(check);
        if (p1.compare("Read grain distribution maps")==0)    checkReadGrainSizeDistribution->setChecked(check);

        if (p1.compare("Number of grain size classes (simulated)")==0)  E_NumberClasses->setValue(val);
        //if (p1.compare("Grain size distribution type")==0)    E_GrainSizeDistributionType->setValue(val);

//        if (p1.compare("Number of grain size classes (maps)")==0)  E_NumberClassesMaps->setValue(val);
        if (p1.compare("Grain size class maps")==0)   E_GrainSizes->setText(p);
        if (p1.compare("Use material depth")==0)             checkMaterialDepth->setChecked(check);


        if (p1.compare("Include Rainfall")==0)               dummyrain = check;//checkRainfall->setChecked(check);
        if (p1.compare("Include Snowmelt")==0)               dummysnow = check;//checkSnowmelt->setChecked(check);
        // if (p1.compare("Alternative flow detachment")==0)    checkAltErosion->setChecked(check);
        // if (p1.compare("Simple depression storage")==0)      checkSimpleDepression->setChecked(check);
        if (p1.compare("Hard Surfaces")==0)                  checkHardsurface->setChecked(check);
        // if (p1.compare("Limit Deposition TC")==0)            checkLimitDepTC->setChecked(check);
        if (p1.compare("Include wheeltracks")==0)            checkInfilCompact->setChecked(check);
        if (p1.compare("Include compacted")==0)            checkInfilCompact->setChecked(check);
        if (p1.compare("Include grass strips")==0)           checkInfilGrass->setChecked(check);
        if (p1.compare("Grassstrip Mannings n")==0)           E_GrassStripN->setText(p);
        if (p1.compare("Include crusts")==0)                 checkInfilCrust->setChecked(check);
        if (p1.compare("Impermeable sublayer")==0)           checkImpermeable->setChecked(check);
        if (p1.compare("Geometric mean Ksat")==0)            checkGeometric->setChecked(check);
        if (p1.compare("Include percolation")==0)            checkPercolation->setChecked(check);

        //	  if (p1.compare("Matric head files")==0)              checkDumphead->setChecked(check);


        if (p1.compare("Timeseries as PCRaster")==0)         checkWritePCRnames->setChecked(check);
        if (p1.compare("plot as PCRaster")==0)           checkWritePCRaster->setChecked(check);
        if (p1.compare("Timeplot as CSV")==0)                checkWriteCommaDelimited->setChecked(check);

        //OBSOLETE
        //    if (p1.compare("Runoff maps in l/s/m")==0)           checkRunoffPerM->setChecked(check);
        //    if (p1.compare("Regular runoff output")==0)          checkOutputTimeStep->setChecked(check);
        //    if (p1.compare("User defined output")==0)            checkOutputTimeUser->setChecked(check);
        //    if (p1.compare("No erosion at outlet")==0)           checkNoErosionOutlet->setChecked(check);
        //    if (p1.compare("Subsoil drainage")==0)               checkDrainage->setChecked(check);
        //    if (p1.compare("Gully infiltration")==0)             checkGullyInfil->setChecked(check);
        //    if (p1.compare("Use initial gully dimensions")==0)   checkGullyInit->setChecked(check);

        if (p1.compare("Report point output separate")==0)   checkSeparateOutput->setChecked(check);
        //if (p1.compare("Report point output for SOBEK")==0)  checkWriteSOBEK->setChecked(check);
        //if (p1.compare("SOBEK date string")==0)              SOBEKdatestring->setText(p);
        if (p1.compare("Use canopy storage map")==0)          radioButton_9->setChecked(check);

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

        if (p1.compare("Infil Method")==0 || p1.compare("Method")==0) //<= very old runfile
        {
            namelist[j].name = QString("Infil Method");
            uiInfilMethod = 0;
            dummy2layerinfil = false;
            switch(iii) {
            case INFIL_SWATRE : uiInfilMethod = 1; break;
            case INFIL_GREENAMPT : uiInfilMethod = 2; break;
            case INFIL_GREENAMPT2 : uiInfilMethod = 2; dummy2layerinfil = true; break;
            case INFIL_SMITH : uiInfilMethod = 3; break;
            case INFIL_SMITH2 : uiInfilMethod = 3;  dummy2layerinfil = true; break;
            case INFIL_KSAT : uiInfilMethod = 4; break;
            }
            E_InfiltrationMethod->setCurrentIndex(uiInfilMethod);
            checkInfil2layer->setChecked(dummy2layerinfil);
        }

        //VJ 110705 KE equations
        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            radioButtonKE1->setChecked(param[0].toInt() == 1);
            spinKEparameterA1->setValue(param[1].toDouble());
            spinKEparameterB1->setValue(param[2].toDouble());
            spinKEparameterC1->setValue(param[3].toDouble());
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            radioButtonKE2->setChecked(param[0].toInt() == 1);
            spinKEparameterA2->setValue(param[1].toDouble());
            spinKEparameterB2->setValue(param[2].toDouble());
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
            param = p.split(",",QString::SkipEmptyParts);
            radioButtonKE3->setChecked(param[0].toInt() == 1);
            spinKEparameterA3->setValue(param[1].toDouble());
            spinKEparameterB3->setValue(param[2].toDouble());
        }
        if (p1.compare("KE time based")==0)      checkKETimebased->setChecked(check);

        if (p1.compare("Ksat calibration")==0)
        {
            if (oldRunfile)
            {
                val/=100;
                QMessageBox::warning(this,"openLISEM",QString("Old runfile detected: calibration value Ksat changed from % to fraction:\n"
                                                              "Ksat calibration divided by 100, check 'Calibration' options in main menu."));
            }
            E_CalibrateKsat->setValue(val);
        }
        if (p1.compare("N calibration")==0)            E_CalibrateN->setValue(val);
        if (p1.compare("Theta calibration")==0)        E_CalibrateTheta->setValue(val);
        if (p1.compare("Psi calibration")==0)          E_CalibratePsi->setValue(val);
        if (p1.compare("Channel Ksat calibration")==0) E_CalibrateChKsat->setValue(val);
        if (p1.compare("Channel N calibration")==0)    E_CalibrateChN->setValue(val);

        if (p1.compare("Erosive Power Calibration")==0)         E_CalibrateEP->setValue(val);
        if (p1.compare("Transport Capacity Calibration")==0)    E_CalibrateTC->setValue(val);
        if (p1.compare("Settling Velocity Calibration")==0)     E_CalibrateSV->setValue(val);
        if (p1.compare("Internal Friction Angle")==0)          E_CalibrateYS->setValue(val);
        if (p1.compare("Dynamic Viscosity Calibration")==0)     E_CalibrateDV->setValue(val);
        if (p1.compare("Drag Force Calibration")==0)            E_CalibrateDF->setValue(val);
        if (p1.compare("Solid Phase Friction Calibration")==0)  E_CalibrateSPF->setValue(val);
        if (p1.compare("Deposition Criteria Calibration")==0)   E_CalibrateDC->setValue(val);
        if (p1.compare("Limit Failure")==0)                     E_LimitFailure->setChecked(check);

        if (p1.compare("Soil Cohesion Calibration")==0)                 E_CalibrateSoilCohesion->setValue(val);
        if (p1.compare("Soil Internal Friction Angle Calibration")==0)  E_CalibrateInternalFrictionAngle->setValue(val);
        if (p1.compare("Soil Depth Calibration")==0)                    E_CalibrateSoilDepth->setValue(val);
        if (p1.compare("Create Stable Initial Safety Factor")==0)       E_IncreaseSlopeStabilityCohesion->setChecked(check);
        if (p1.compare("Minimum Safety Factor Calibration")==0)         E_CalibrateMinimumSafetyFactor->setValue(val);



        if (p1.compare("Splash Delivery Ratio")==0)    E_SplashDelibery->setValue(val);
        if (p1.compare("Particle Cohesion of Deposited Layer")==0) E_DepositedCohesion->setValue(val);
        if (p1.compare("Stemflow fraction")==0)        E_StemflowFraction->setValue(val);
        if (p1.compare("Canopy Openess")==0)           E_CanopyOpeness->setValue(val);
        // VJ 110209 canopy openess, factor Aston as user input
        //   if (p1.compare("Max flood level")==0)          E_maxFloodLevel->setValue(val);
        //   if (p1.compare("Min flood dt")==0)             E_minFloodDt->setValue(val);

        //3d display settings

        if(p1.compare("Light_Ambient_R") == 0)
        {
            qDebug() << namelist[j].value << " set  " << val;
            this->GL_Light_Ambient_R->setValue(val);
        }
        if(p1.compare("Light_Ambient_G") == 0)  this->GL_Light_Ambient_G->setValue(val);
        if(p1.compare("Light_Ambient_B") == 0)  this->GL_Light_Ambient_B->setValue(val);
        if(p1.compare("Light_Ambient_A") == 0)  this->GL_Light_Ambient_A->setValue(val);
        if(p1.compare("Light_Directional_R") == 0)  this->GL_Light_Directional_R->setValue(val);
        if(p1.compare("Light_Directional_G") == 0)  this->GL_Light_Directional_G->setValue(val);
        if(p1.compare("Light_Directional_B") == 0)  this->GL_Light_Directional_B->setValue(val);
        if(p1.compare("Light_Directional_A") == 0)  this->GL_Light_Directional_A->setValue(val);
        if(p1.compare("Light_Directional_X") == 0)  this->GL_Light_Directional_X->setValue(val);
        if(p1.compare("Light_Directional_Y") == 0)  this->GL_Light_Directional_Y->setValue(val);
        if(p1.compare("Light_Directional_Z") == 0)  this->GL_Light_Directional_Z->setValue(val);

        /*if(p1.compare("Surface_Draw") == 0)  this->GL_Surface_Draw->setChecked(check);
        if(p1.compare("Surface_Micro_Elevation_Scale") == 0)  this->GL_Surface_Micro_Elevation_Scale->setValue(val);
        if(p1.compare("Surface_Mipmap_Distance_1") == 0)  this->GL_Surface_Mipmap_Distance_1->setValue(val);
        if(p1.compare("Surface_Mipmap_Distance_2") == 0)  this->GL_Surface_Mipmap_Distance_2->setValue(val);
        if(p1.compare("Surface_Vegetated_Small_Color_R") == 0)  this->GL_Surface_Vegetated_Small_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Small_Color_G") == 0)  this->GL_Surface_Vegetated_Small_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Small_Color_B") == 0)  this->GL_Surface_Vegetated_Small_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Large_Color_R") == 0)  this->GL_Surface_Vegetated_Large_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Large_Color_G") == 0)  this->GL_Surface_Vegetated_Large_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Large_Color_B") == 0)  this->GL_Surface_Vegetated_Large_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Bare_Color_R") == 0)  this->GL_Surface_Bare_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Bare_Color_G") == 0)  this->GL_Surface_Bare_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Bare_Color_B") == 0)  this->GL_Surface_Bare_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Roads_Color_R") == 0)  this->GL_Surface_Road_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Roads_Color_G") == 0)  this->GL_Surface_Road_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Roads_Color_B") == 0)  this->GL_Surface_Road_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Buildings_Color_R") == 0)  this->GL_Surface_Building_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Buildings_Color_G") == 0)  this->GL_Surface_Building_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Buildings_Color_B") == 0)  this->GL_Surface_Building_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Erosion_Color_R") == 0)  this->GL_Surface_Erosion_Color_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Erosion_Color_G") == 0)  this->GL_Surface_Erosion_Color_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Erosion_Color_B") == 0)  this->GL_Surface_Erosion_Color_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Erosion_Color_A") == 0)  this->GL_Surface_Erosion_Color_A->setValue(val);
        if(p1.compare("Surface_Vegetated_Deposition_Color_R") == 0)  this->GL_Surface_Deposition_Color_R->setValue(val);
        if(p1.compare("Surface_Vegetated_Deposition_Color_G") == 0)  this->GL_Surface_Deposition_Color_G->setValue(val);
        if(p1.compare("Surface_Vegetated_Deposition_Color_B") == 0)  this->GL_Surface_Deposition_Color_B->setValue(val);
        if(p1.compare("Surface_Vegetated_Deposition_Color_A") == 0)  this->GL_Surface_Deposition_Color_A->setValue(val);

        if(p1.compare("Water_Draw") == 0)  this->GL_Water_Draw->setChecked(check);
        if(p1.compare("Water_Reflectivity") == 0)  this->GL_Water_Reflectivity->setValue(val);
        if(p1.compare("Water_Refractivity") == 0)  this->GL_Water_Refractivity->setValue(val);
        if(p1.compare("Water_Velocity_Scale") == 0)  this->GL_Water_Velocity_Scale->setValue(val);
        if(p1.compare("Water_Micro_Elevation_Scale") == 0)  this->GL_Water_Micro_Elevation_Scale->setValue(val);
        if(p1.compare("Water_Transparancy") == 0)  this->GL_Water_Transparancy->setValue(val);
        if(p1.compare("Water_Deep_Color_R") == 0)  this->GL_Water_Deep_R->setValue(val);
        if(p1.compare("Water_Deep_Color_G") == 0)  this->GL_Water_Deep_G->setValue(val);
        if(p1.compare("Water_Deep_Color_B") == 0)  this->GL_Water_Deep_B->setValue(val);
        if(p1.compare("Water_Deep_Color_A") == 0)  this->GL_Water_Deep_A->setValue(val);
        if(p1.compare("Water_Shallow_Color_R") == 0)  this->GL_Water_Shallow_R->setValue(val);
        if(p1.compare("Water_Shallow_Color_G") == 0)  this->GL_Water_Shallow_G->setValue(val);
        if(p1.compare("Water_Shallow_Color_B") == 0)  this->GL_Water_Shallow_B->setValue(val);
        if(p1.compare("Water_Shallow_Color_A") == 0)  this->GL_Water_Shallow_A->setValue(val);
        if(p1.compare("Water_Sediment_Color_R") == 0)  this->GL_Water_Sediment_R->setValue(val);
        if(p1.compare("Water_Sediment_Color_G") == 0)  this->GL_Water_Sediment_G->setValue(val);
        if(p1.compare("Water_Sediment_Color_B") == 0)  this->GL_Water_Sediment_B->setValue(val);
        if(p1.compare("Water_Sediment_Color_A") == 0)  this->GL_Water_Sediment_A->setValue(val);*/

        if(p1.compare("Clouds_Draw") == 0)  this->GL_Objects_Clouds_Draw->setChecked(check);
        if(p1.compare("Rain_Draw") == 0)  this->GL_Objects_Rain_Draw->setChecked(check);
        if(p1.compare("Roads_Draw") == 0)  this->GL_Objects_Roads_Draw->setChecked(check);
        //if(p1.compare("Roads_Distance") == 0)  this->GL_Objects_Roads_Distance->setValue(val);

        if(p1.compare("Buildings_Draw") == 0)  this->GL_Objects_Buildings_Draw->setChecked(check);
        if(p1.compare("Buildings_Distance") == 0)  this->GL_Objects_Buildings_Distance->setValue(val);

        if(p1.compare("Trees_Draw") == 0)  this->GL_Objects_Trees_Draw->setChecked(check);
        if(p1.compare("Trees_Distance") == 0)  this->GL_Objects_Trees_Distance->setValue(val);
        //if(p1.compare("Trees_Instances") == 0)  this->GL_Objects_Trees_Instances->setValue(val);
        if(p1.compare("Trees_Increment") == 0)  this->GL_Objects_Trees_Increment->setValue(val);

        if(p1.compare("Grass_Draw") == 0)  this->GL_Objects_Grass_Draw->setChecked(check);
        if(p1.compare("Grass_Distance") == 0)  this->GL_Objects_Grass_Distance->setValue(val);
        //if(p1.compare("Grass_Instances") == 0)  this->GL_Objects_Grass_Instances->setValue(val);
        if(p1.compare("Grass_Increment") == 0)  this->GL_Objects_Grass_Increment->setValue(val);
        //if(p1.compare("Grass_Vertical_Scale") == 0) this->GL_Objects_Grass_Vertical_Scale->setValue(val);



        //VJ 111120 water repellency
        if (p1.compare("Use Water Repellency")==0)      checkWaterRepellency->setChecked(check);
        if (p1.compare("Water Repellency A")==0)        E_waterRep_a->setValue(val);
        if (p1.compare("Water Repellency B")==0)        E_waterRep_b->setValue(val);
        if (p1.compare("Water Repellency C")==0)        E_waterRep_c->setValue(val);
        if (p1.compare("Water Repellency D")==0)        E_waterRep_d->setValue(val);


        if (p1.compare("Output interval")==0)   printinterval->setValue(std::max(1,iii));

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

            if (outputcheck.count() < 26)
               for (int k = outputcheck.count();k < 26; k++)
                   outputcheck << "0";

            checkBox_OutRunoff->setChecked(bool(outputcheck.at(0).toInt() == 1));
            checkBox_OutWH->setChecked(bool(outputcheck.at(2).toInt() == 1));
            //OBSOLETE checkBox_OutWHC->setChecked(bool(outputcheck.at(3).toInt() == 1));
            checkBox_OutInf->setChecked(bool(outputcheck.at(8).toInt() == 1));
            checkBox_OutV->setChecked(bool(outputcheck.at(7).toInt() == 1));
            checkBox_OutDet->setChecked(bool(outputcheck.at(5).toInt() == 1));
            checkBox_OutDep->setChecked(bool(outputcheck.at(6).toInt() == 1));
            checkBox_OutConc->setChecked(bool(outputcheck.at(1).toInt() == 1));
            checkBox_OutTC->setChecked(bool(outputcheck.at(4).toInt() == 1));
            checkBox_OutSurfStor->setChecked(bool(outputcheck.at(9).toInt() == 1));
            checkBox_OutChanVol->setChecked(bool(outputcheck.at(10).toInt() == 1));
            checkBox_OutTiledrain->setChecked(bool(outputcheck.at(11).toInt() == 1));
            checkBox_OutHmx->setChecked(bool(outputcheck.at(12).toInt() == 1));
            checkBox_OutQf->setChecked(bool(outputcheck.at(13).toInt() == 1));
            checkBox_OutVf->setChecked(bool(outputcheck.at(14).toInt() == 1));
            checkBox_OutHmxWH->setChecked(bool(outputcheck.at(15).toInt() == 1));
            checkBox_OutSL->setChecked(bool(outputcheck.at(16).toInt() == 1));
            checkBox_OutSed->setChecked(bool(outputcheck.at(17).toInt() == 1));

            checkBox_OutSafetyFactor->setChecked(bool(outputcheck.at(18).toInt() == 1));
            checkBox_OutSlopeFailure->setChecked(bool(outputcheck.at(19).toInt() == 1));
            checkBox_OutDebrisFlowHeight->setChecked(bool(outputcheck.at(20).toInt() == 1));
            checkBox_OutDebrisFlowVelocity->setChecked(bool(outputcheck.at(21).toInt() == 1));
            checkBox_OutFluidPhaseHeight->setChecked(bool(outputcheck.at(22).toInt() == 1));
            checkBox_OutSolidPhaseHeight->setChecked(bool(outputcheck.at(23).toInt() == 1));
            checkBox_OutEntrainment->setChecked(bool(outputcheck.at(24).toInt() == 1));
            checkBox_OutTimestep->setChecked(bool(outputcheck.at(25).toInt() == 1));









            // TODO replace these numbers with defines for clarity in model.h and lisemqt.h

            // checkboxes normal output map series, numbering according to original LISEM

            if (p1.compare("WH max level map")==0) E_WHmaxMap->setText(p);

        }
    }

    checkDoErosion->setChecked(dummyErosion);
    setErosionTab(dummyErosion);

    if (!dummyrain && !dummysnow)
        QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");

    checkRainfall->setChecked(dummyrain);
    checkSnowmelt->setChecked(dummysnow);

    // qDebug() << dummyFloodExplicit << dummyFloodSWOF1 << dummyFloodSWOF2 << E_floodSolution->value();
    // get directory and file names

    // first guess
    E_WorkDir = QFileInfo(E_runFileList->currentText()).dir().absolutePath();
    QDir dir(E_WorkDir);
    if (dir.cdUp())
        E_WorkDir = dir.absolutePath()+"/";

    for (j = 0; j < nrnamelist; j++)
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
        if (p1.compare("Work Directory")==0)
        {
            QString S = E_WorkDir;
            E_WorkDir = CheckDir(p);
            if (!QFileInfo(E_WorkDir).exists())
                E_WorkDir = S;
        }

        if (p1.compare("Map Directory")==0)
        {
            E_MapDir->setText(CheckDir(p));

            if (!p.isEmpty() && E_WorkDir.isEmpty())
            {
                E_WorkDir = E_MapDir->text();
                QDir dir(E_WorkDir);
                if (dir.cdUp())
                    E_WorkDir = dir.absolutePath()+"/";

            }

            if (E_MapDir->text().isEmpty() && !E_WorkDir.isEmpty())
            {
                E_MapDir->setText(E_WorkDir+"maps/");
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
            if (!QFileInfo(E_ResultDir->text()).exists())
                E_ResultDir->setText(E_WorkDir/*->text()*/+"res/");

        }
        if (p1.compare("Main results file")==0) E_MainTotals->setText(p);
        if (p1.compare("Filename point output")==0) E_PointResults->setText(p);
        if (p1.compare("Filename landunit output")==0) E_LandunitResults->setText(p);
        // resultDir is added in report operation

        if (p1.compare("Rainfall Directory")==0) RainFileDir = CheckDir(p);
        if (p1.compare("Rainfall file")==0)
        {
            E_RainfallName->setText(RainFileDir + p);
            RainFileName = p;///*rainFileDir + */E_RainfallName->text();
            if (!QFileInfo(E_RainfallName->text()).exists())
            {
                RainFileDir = QString(E_WorkDir/*->text()*/ + "rain/");
                E_RainfallName->setText(RainFileDir + p);
            }
        }

        if (p1.compare("Rainfall map")==0) E_RainfallMap->setText(p);
        if (p1.compare("Interception map")==0) E_InterceptionMap->setText(p);
        if (p1.compare("Infiltration map")==0) E_InfiltrationMap->setText(p);
        if (p1.compare("Runoff map")==0) E_RunoffMap->setText(p);
        if (p1.compare("Runoff fraction map")==0) E_RunoffFractionMap->setText(p);
        if (p1.compare("Channel discharge map")==0) E_ChannelQtotm3Map->setText(p);

        if (p1.compare("Erosion map")==0) E_DetachmentMap->setText(p);
        if (p1.compare("Deposition map")==0) E_DepositionMap->setText(p);
        if (p1.compare("Soilloss map")==0) E_SoillossMap->setText(p);
        if (p1.compare("Channel detachment map")==0) E_ChanDetachmentMap->setText(p);
        if (p1.compare("Channel deposition map")==0) E_ChanDepositionMap->setText(p);

        if (p1.compare("Maximum Debris Flow Height Map")==0) E_MaxDebrisFlowHeight->setText(p);
        if (p1.compare("Maximum Debris Flow Velocity Map")==0) E_MaxDebrisFlowVelocity->setText(p);
        if (p1.compare("Debris Flow Start Map")==0) E_DebrisFlowStart->setText(p);
        if (p1.compare("Entrainment Map")==0) E_Entrainment->setText(p);
        if (p1.compare("Slope Failure Map")==0) E_SlopeFailure->setText(p);
        if (p1.compare("Minimum Safety Factor Map")==0) E_MinimumSafetyFactor->setText(p);

        if (p1.compare("Flood level map")==0) E_FloodlevelMap->setText(p);
        if (p1.compare("Flood time map")==0) E_FloodTimeMap->setText(p);
        if (p1.compare("Flood start time")==0) E_FloodFEW->setText(p);
        if (p1.compare("Flood Max V")==0) E_FloodmaxVMap->setText(p);
        if (p1.compare("Channel Max Q")==0) E_ChannelMaxQ->setText(p);
        if (p1.compare("Channel Max WH")==0) E_ChannelMaxWH->setText(p);
        if (p1.compare("Flood stats")==0) E_FloodStats->setText(p);

        if (p1.compare("Snowmelt Directory")==0) SnowmeltFileDir = CheckDir(p);
        if (p1.compare("Snowmelt file")==0)
        {
            E_SnowmeltName->setText(SnowmeltFileDir + p);
            SnowmeltFileName = p;// /*SnowmeltFileDir + */E_SnowmeltName->text();
        }

        if (uiInfilMethod == 1 && p1.compare("Table Directory")==0)
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


    //checkBox_OutHmx->setEnabled(checkChannelFlood->isChecked());
    //checkBox_OutQf->setEnabled(checkChannelFlood->isChecked());
    //checkBox_OutVf->setEnabled(checkChannelFlood->isChecked());
    //checkBox_OutHmxWH->setEnabled(checkChannelFlood->isChecked());

    on_checkIncludeChannel_clicked();
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
    // the maplist contains the default names already
    // this is to get the correct names for the model run
    fillNamelistMapnames(false);
    for (int k = 0; k < nrmaplist; k++)
        mapList[k].dir = E_MapDir->text();
    // dir necessary?

    //RunAllChecks();
    //obsolete: is done in void lisemqt::on_E_runFileList_currentIndexChanged(int)
}
//---------------------------------------------------------------------------
QString lisemqt::CheckDir(QString p, bool makeit)
{
    /* TODO mulitplatform: fromNativeSeparators etc*/
    QString path;

    if (p.isEmpty())
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
            qDebug() << "NOTE: Result dir created !";
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
    if(!checkRainfall->isChecked() && !checkSnowmelt->isChecked())
        QMessageBox::warning(this,"openLISEM","No rainfall or snowmelt, running on empty!");


    for (int j = 0; j < nrnamelist; j++)
    {
        QString p1 = namelist[j].name;
        QString p;
        // erosion
         if (p1.compare("No Erosion simulation")==0)           namelist[j].value.setNum((int)(!checkDoErosion->isChecked()));
        // obsolete
        if (p1.compare("Include Erosion simulation")==0)      namelist[j].value.setNum((int)checkDoErosion->isChecked());

        //channels
        if (p1.compare("Include main channels")==0)          namelist[j].value.setNum((int)checkIncludeChannel->isChecked());
        if (p1.compare("Include channel infil")==0)          namelist[j].value.setNum((int)checkChannelInfil->isChecked());
        if (p1.compare("Include channel baseflow")==0)       namelist[j].value.setNum((int)checkChannelBaseflow->isChecked());

        //flooding
   //     if (p1.compare("Include channel flooding")==0)       namelist[j].value.setNum((int)checkChannelFlood->isChecked());
   //     if (p1.compare("Include rainfall flooding")==0)      namelist[j].value.setNum((int)checkRainfallFlood->isChecked());
   //     if (p1.compare("Rainfall flooding gradient")==0)     namelist[j].value = E_RainFloodGradient->text();
        if (p1.compare("Include road system")==0)            namelist[j].value.setNum((int)checkRoadsystem->isChecked());

        if (p1.compare("Include litter interception")==0)    namelist[j].value.setNum((int)checkIncludeLitter->isChecked());
        if (p1.compare("Litter interception storage")==0)    namelist[j].value = E_LitterSmax->text();

        //Surface Flow
        if (p1.compare("Enable Solid Phase")==0)                                    namelist[j].value.setNum((int)checkSolidPhase->isChecked());
        if (p1.compare("Enable Entrainment")==0)                                    namelist[j].value.setNum((int)checkEntrainment->isChecked());
        if (p1.compare("Enable Deposition")==0)                                     namelist[j].value.setNum((int)checkBox_Deposition->isChecked());
        if (p1.compare("Enable Compaction")==0)                                     namelist[j].value.setNum((int)checkBox_Compaction->isChecked());
        if (p1.compare("Enable Seismic") == 0)                                      namelist[j].value.setNum((int)checkSeismic->isChecked());
        if (p1.compare("Enable Slope Stability")==0)                                namelist[j].value.setNum((int)checkSlopeStability->isChecked());
        if (p1.compare("Enable Upslope Forcing") == 0)                              namelist[j].value.setNum((int)checkUpslopeForcing->isChecked());
        if (p1.compare("Enable Downslope Forcing") == 0)                              namelist[j].value.setNum((int)checkDownslopeForcing->isChecked());
        if (p1.compare("Include Bedrock Layer") == 0)                              namelist[j].value.setNum((int)checkBedRockLayer->isChecked());
        if (p1.compare("Enable Slope Failure")==0)                                  namelist[j].value.setNum((int)checkSlopeFailure->isChecked());
        if (p1.compare("Minimum Safety Factor")==0)                                 namelist[j].value = E_PRESF->text();
        if (p1.compare("Maximum Safety Factor")==0)                                 namelist[j].value = E_POSTSF->text();
        if (p1.compare("Maximum safety factor for display")==0)                     namelist[j].value = E_SFMaxDisplay->text();
        if (p1.compare("Entrainment Coefficient")==0)                               namelist[j].value = E_EntrainmentCoefficient->text();
        if (p1.compare("Minimum Entrainment Height")==0)                            namelist[j].value = E_EntrainmentMinimumWaterHeight->text();
        if (p1.compare("Minimum Failure Height")==0)                                namelist[j].value = E_MinimumFailureHeight->text();
        //if (p1.compare("Enable Levees")==0)                                         namelist[j].value.setNum((int)checkLevees->isChecked());
        if (p1.compare("Enable Barriers")==0)                                       namelist[j].value.setNum((int)checkBox_Barriers->isChecked());
        if (p1.compare("Enable Flow Barriers")==0)                                  namelist[j].value.setNum((int)checkBox_FlowBarriers->isChecked());
        if (p1.compare("Flow barrier table filename")==0)                           namelist[j].value = lineEdit_FlowBarriers->text();
        if (p1.compare("Enable Inflow")==0)                                   namelist[j].value.setNum((int)checkBox_Inflow->isChecked());
        if (p1.compare("Inflow table filename")==0)                           namelist[j].value = lineEdit_Inflow->text();
        if (p1.compare("Include Initial FluidSolid Mixture")==0)                    namelist[j].value.setNum((int)checkBox_UFInitial->isChecked());
        if (p1.compare("Include Forced FluidSolid Mixture")==0)                     namelist[j].value.setNum((int)checkBox_UFForced->isChecked());
        if (p1.compare("Incldue Maximum ChannelVolume")==0)                         namelist[j].value.setNum((int)checkBox_ChannelMaxVolume->isChecked());
        if (p1.compare("Incldue Maximum Volume")==0)                                namelist[j].value.setNum((int)checkBox_MaxVolume->isChecked());

        if (p1.compare("Suspended Viscosity")==0)                                    namelist[j].value.setNum((int)checkBox_SuspendedViscosity->isChecked());
        if (p1.compare("Lax Multiplier")==0)                                        namelist[j].value = E_LaxMultiplier->text();
        if (p1.compare("Friction force correction")==0)                             namelist[j].value = E_FrictionCorrection->text();
        if (p1.compare("Erosion Cohesion Calibration")==0)                          namelist[j].value = E_CalibrateESC->text();
        if (p1.compare("Erosion Grain Size Calibration")==0)                        namelist[j].value = E_CalibrateEGS->text();

        if (p1.compare("Solid-Fluid Drag Coefficient")==0)                          namelist[j].value = E_UFDrag->text();
        if (p1.compare("Use HLL2")==0)                                              namelist[j].value.setNum((int)checkBox_HLL2->isChecked());
        if (p1.compare("Spatially Dynamic Timestep")==0)                            namelist[j].value.setNum((int)checkBox_SpatiallyDynamicTimestep->isChecked());
        if (p1.compare("Flow Minimum Timestep")==0)                                 namelist[j].value = E_SFMinimumDT->text();
        if (p1.compare("Surface Flow Courant Factor")==0)                           namelist[j].value = E_SFCourant->text();

        if (p1.compare("Drag Power Law Coefficient")==0)                            namelist[j].value = E_UFDragPower->text();
        if (p1.compare("Viscosity Alpha")==0)                                       namelist[j].value = E_UFViscosityAlpha->text();
        if (p1.compare("Viscosity Beta")==0)                                        namelist[j].value = E_UFViscosityBeta->text();

        if (p1.compare("Limit Cores")==0)                                           namelist[j].value.setNum((int)checkBox_CoreLimit->isChecked());
        if (p1.compare("Core Limit")==0)                                            namelist[j].value = E_CoreLimit->text();

        if (p1.compare("Kinematic Timestep Power")==0)                              namelist[j].value = E_KinematicTimestepPower->text();

        if (p1.compare("Minimal Flood Water Depth")==0)                             namelist[j].value = E_FloodMinDepth->text();
        if (p1.compare("Minimum Debris Flow Volumetric Sediment Fraction")==0)      namelist[j].value = E_DebrisFlowMinSedConcentration->text();



        if (p1.compare("Sigma diffusion")==0)                namelist[j].value = E_SigmaDiffusion->text();


        if (p1.compare("Advanced sediment")==0)               namelist[j].value.setNum((int)checkAdvancedSediment->isChecked());

        if (p1.compare("Detachment efficiency")==0)          namelist[j].value = E_EfficiencyDET->text();
   //     if (p1.compare("Detachment stoniness")==0)          namelist[j].value.setNum((int)checkStoninessDET->isChecked());

        if (p1.compare("Advanced sediment configuration")==0)
        {
            if(checkBox_SedMultiMulti->isChecked())
            {
                namelist[j].value.setNum((int)2);
            }else
                if(checkBox_SedMultiSingle->isChecked())
                {
                    namelist[j].value.setNum((int)1);
                }else
                {
                    namelist[j].value.setNum((int)0);
                }
        }

        if (p1.compare("BL method")==0)                 namelist[j].value = E_BLMethod->text();
        if (p1.compare("SS method")==0)                 namelist[j].value = E_SSMethod->text();

        if (p1.compare("Estimate grain size distribution")==0)namelist[j].value.setNum((int)checkEstimateGrainSizeDistribution->isChecked());

        if (p1.compare("Read grain distribution maps")==0)    namelist[j].value.setNum((int)checkReadGrainSizeDistribution->isChecked());
        if (p1.compare("Number of grain size classes (simulated)")==0)  namelist[j].value = E_NumberClasses->text();
       // if (p1.compare("Grain size distribution type")==0)    namelist[j].value = E_GrainSizeDistributionType->text();

        //if (p1.compare("Number of grain size classes (maps)")==0)  namelist[j].value = E_NumberClassesMaps->text();
        if (p1.compare("Grain size class maps")==0)   namelist[j].value = E_GrainSizes->text();
        if (p1.compare("Use material depth")==0)             namelist[j].value.setNum((int)checkMaterialDepth->isChecked());

        //tile drains
        if (p1.compare("Include tile drains")==0)            namelist[j].value.setNum((int)checkIncludeTiledrains->isChecked());

        //houses
        if (p1.compare("Include house storage")==0)          namelist[j].value.setNum((int)checkHouses->isChecked());
        if (p1.compare("Include raindrum storage")==0)       namelist[j].value.setNum((int)checkRaindrum->isChecked());

        if (p1.compare("Include Rainfall")==0)               namelist[j].value.setNum((int)checkRainfall->isChecked());
        if (p1.compare("Include Snowmelt")==0)               namelist[j].value.setNum((int)checkSnowmelt->isChecked());
        if (p1.compare("Hard Surfaces")==0)                  namelist[j].value.setNum((int)checkHardsurface->isChecked());
        if (p1.compare("Include compacted")==0)            namelist[j].value.setNum((int)checkInfilCompact->isChecked());
        if (p1.compare("Include grass strips")==0)           namelist[j].value.setNum((int)checkInfilGrass->isChecked());
        if (p1.compare("Grassstrip Mannings n")==0)          namelist[j].value = E_GrassStripN->text();

        if (p1.compare("Include crusts")==0)                 namelist[j].value.setNum((int)checkInfilCrust->isChecked());
        if (p1.compare("Impermeable sublayer")==0)           namelist[j].value.setNum((int)checkImpermeable->isChecked());
        //if (p1.compare("Matric head files")==0)              namelist[j].value.setNum((int)checkDumphead->isChecked());
        if (p1.compare("Include percolation")==0)                 namelist[j].value.setNum((int)checkPercolation->isChecked());
        if (p1.compare("Geometric mean Ksat")==0)            namelist[j].value.setNum((int)checkGeometric->isChecked());
        if (p1.compare("Timeseries as PCRaster")==0)         namelist[j].value.setNum((int)checkWritePCRnames->isChecked());
        if (p1.compare("Timeseries as CSV")==0)              namelist[j].value.setNum((int)checkWriteCommaDelimited->isChecked());
        if (p1.compare("Timeplot as PCRaster")==0)           namelist[j].value.setNum((int)checkWritePCRaster->isChecked());
        //if (p1.compare("Regular runoff output")==0)          namelist[j].value.setNum((int)checkOutputTimeStep->isChecked());
        //if (p1.compare("User defined output")==0)            namelist[j].value.setNum((int)checkOutputTimeUser->isChecked());
        //if (p1.compare("No erosion at outlet")==0)           namelist[j].value.setNum((int)checkNoErosionOutlet->isChecked());
        if (p1.compare("Report point output separate")==0)   namelist[j].value.setNum((int)checkSeparateOutput->isChecked());
        //if (p1.compare("Report point output for SOBEK")==0)  namelist[j].value.setNum((int)checkWriteSOBEK->isChecked());
        //if (p1.compare("SOBEK date string")==0)              namelist[j].value = SOBEKdatestring->text();
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

        //VJ 110705 KE equations
        if (p1.compare("KE parameters EQ1")==0)
        {
            QStringList param;
            param << (radioButtonKE1->isChecked()?"1":"0") << spinKEparameterA1->text() << spinKEparameterB1->text() << spinKEparameterC1->text();
            namelist[j].value = param.join(",");
        }
        if (p1.compare("KE parameters EQ2")==0)
        {
            QStringList param;
            param << (radioButtonKE2->isChecked()?"1":"0") << spinKEparameterA2->text() << spinKEparameterB2->text();
            namelist[j].value = param.join(",");
        }
        if (p1.compare("KE parameters EQ3")==0)
        {
            QStringList param;
            param << (radioButtonKE3->isChecked()?"1":"0") << spinKEparameterA3->text() << spinKEparameterB3->text();
            namelist[j].value = param.join(",");
        }
        if (p1.compare("KE time based")==0)      namelist[j].value.setNum((int)checkKETimebased->isChecked());

        if (p1.compare("Begin time")==0) namelist[j].value = E_BeginTime->text();
        if (p1.compare("End time")==0)   namelist[j].value = E_EndTime->text();
        if (p1.compare("Timestep")==0)   namelist[j].value = E_Timestep->text();
        if (p1.compare("Work Directory")==0)    namelist[j].value = E_WorkDir;//->text();
        if (p1.compare("Map Directory")==0)    namelist[j].value = E_MapDir->text();
        if (p1.compare("Result Directory")==0) namelist[j].value = E_ResultDir->text();
        if (p1.compare("Main results file")==0) namelist[j].value = E_MainTotals->text();
        if (p1.compare("Filename point output")==0) namelist[j].value = E_PointResults->text();
        if (p1.compare("Filename landunit output")==0) namelist[j].value = E_LandunitResults->text();
        if (p1.compare("Rainfall Directory")==0) namelist[j].value = RainFileDir;
        if (p1.compare("Rainfall file")==0) namelist[j].value = RainFileName; //E_RainfallName->text();
        if (p1.compare("Snowmelt Directory")==0) namelist[j].value = SnowmeltFileDir;
        if (p1.compare("Snowmelt file")==0) namelist[j].value = SnowmeltFileName;//E_SnowmeltName->text();

        if (p1.compare("Rainfall map")==0) namelist[j].value = E_RainfallMap->text();
        if (p1.compare("Interception map")==0) namelist[j].value = E_InterceptionMap->text();
        if (p1.compare("Infiltration map")==0) namelist[j].value = E_InfiltrationMap->text();
        if (p1.compare("Runoff map")==0) namelist[j].value = E_RunoffMap->text();
        if (p1.compare("Runoff fraction map")==0) namelist[j].value = E_RunoffFractionMap->text();
        if (p1.compare("Channel discharge map")==0) namelist[j].value = E_ChannelQtotm3Map->text();
        if (p1.compare("WH max level map")==0) namelist[j].value = E_WHmaxMap->text();

        if (p1.compare("Flood level map")==0) namelist[j].value = E_FloodlevelMap->text();
        if (p1.compare("Flood time map")==0) namelist[j].value = E_FloodTimeMap->text();
        if (p1.compare("Flood start time")==0) namelist[j].value = E_FloodFEW->text();
        if (p1.compare("Flood Max V")==0) namelist[j].value = E_FloodmaxVMap->text();
        if (p1.compare("Channel Max Q")==0) namelist[j].value = E_ChannelMaxQ->text();
        if (p1.compare("Channel Max WH")==0) namelist[j].value = E_ChannelMaxWH->text();
        if (p1.compare("Flood stats")==0) namelist[j].value = E_FloodStats->text();

        if (p1.compare("Erosion map")==0) namelist[j].value = E_DetachmentMap->text();
        if (p1.compare("Deposition map")==0) namelist[j].value = E_DepositionMap->text();
        if (p1.compare("Soilloss map")==0) namelist[j].value = E_SoillossMap->text();
        if (p1.compare("Channel detachment map")==0) namelist[j].value = E_ChanDetachmentMap->text();
        if (p1.compare("Channel deposition map")==0) namelist[j].value = E_ChanDepositionMap->text();

        if (p1.compare("Maximum Debris Flow Height Map")==0) namelist[j].value = E_MaxDebrisFlowHeight->text();
        if (p1.compare("Maximum Debris Flow Velocity Map")==0) namelist[j].value = E_MaxDebrisFlowVelocity->text();
        if (p1.compare("Debris Flow Start Map")==0) namelist[j].value = E_DebrisFlowStart->text();
        if (p1.compare("Entrainment Map")==0) namelist[j].value = E_Entrainment->text();
        if (p1.compare("Slope Failure Map")==0) namelist[j].value = E_SlopeFailure->text();
        if (p1.compare("Minimum Safety Factor Map")==0) namelist[j].value = E_MinimumSafetyFactor->text();

        if (p1.compare("Ksat calibration")==0) namelist[j].value = E_CalibrateKsat->text();
        if (p1.compare("N calibration")==0) namelist[j].value = E_CalibrateN->text();
        if (p1.compare("Theta calibration")==0) namelist[j].value = E_CalibrateTheta->text();
        if (p1.compare("Psi calibration")==0) namelist[j].value = E_CalibratePsi->text();
        if (p1.compare("Channel Ksat calibration")==0) namelist[j].value = E_CalibrateChKsat->text();
        if (p1.compare("Channel N calibration")==0) namelist[j].value = E_CalibrateChN->text();

        if (p1.compare("Erosive Power Calibration")==0)         namelist[j].value = E_CalibrateEP->text();
        if (p1.compare("Transport Capacity Calibration")==0)    namelist[j].value = E_CalibrateTC->text();
        if (p1.compare("Settling Velocity Calibration")==0)     namelist[j].value = E_CalibrateSV->text();
        if (p1.compare("Internal Friction Angle")==0)          namelist[j].value = E_CalibrateYS->text();
        if (p1.compare("Dynamic Viscosity Calibration")==0)     namelist[j].value = E_CalibrateDV->text();
        if (p1.compare("Drag Force Calibration")==0)            namelist[j].value = E_CalibrateDF->text();
        if (p1.compare("Solid Phase Friction Calibration")==0)  namelist[j].value = E_CalibrateSPF->text();
        if (p1.compare("Deposition Criteria Calibration")==0)   namelist[j].value = E_CalibrateDC->text();
        if (p1.compare("Limit Failure")==0)                     namelist[j].value.setNum((int)E_LimitFailure->isChecked());

        if (p1.compare("Soil Cohesion Calibration")==0)                 namelist[j].value = E_CalibrateSoilCohesion->text();
        if (p1.compare("Soil Internal Friction Angle Calibration")==0)  namelist[j].value = E_CalibrateInternalFrictionAngle->text();
        if (p1.compare("Soil Depth Calibration")==0)                    namelist[j].value = E_CalibrateSoilDepth->text();
        if (p1.compare("Create Stable Initial Safety Factor")==0)       namelist[j].value.setNum((int)E_IncreaseSlopeStabilityCohesion->isChecked());
        if (p1.compare("Minimum Safety Factor Calibration")==0)         namelist[j].value = E_CalibrateMinimumSafetyFactor->text();


        if (p1.compare("Splash Delivery Ratio")==0) namelist[j].value = E_SplashDelibery->text();
        if (p1.compare("Particle Cohesion of Deposited Layer")==0) namelist[j].value = E_DepositedCohesion->text();
        if (p1.compare("Stemflow fraction")==0) namelist[j].value = E_StemflowFraction->text();
        if (p1.compare("Canopy Openess")==0) namelist[j].value = E_CanopyOpeness->text();
        // VJ 110209 canopy openess, factor Aston as user input

        //3d display settings
        if(p1.compare("Light_Ambient_R") == 0)
        {
            qDebug() << namelist[j].value << "  " << this->GL_Light_Ambient_R->text();
            namelist[j].value = this->GL_Light_Ambient_R->text();
        }
        if(p1.compare("Light_Ambient_G") == 0)  namelist[j].value = this->GL_Light_Ambient_G->text();
        if(p1.compare("Light_Ambient_B") == 0)  namelist[j].value = this->GL_Light_Ambient_B->text();
        if(p1.compare("Light_Ambient_A") == 0)  namelist[j].value = this->GL_Light_Ambient_A->text();
        if(p1.compare("Light_Directional_R") == 0)  namelist[j].value = this->GL_Light_Directional_R->text();
        if(p1.compare("Light_Directional_G") == 0)  namelist[j].value = this->GL_Light_Directional_G->text();
        if(p1.compare("Light_Directional_B") == 0)  namelist[j].value = this->GL_Light_Directional_B->text();
        if(p1.compare("Light_Directional_A") == 0)  namelist[j].value = this->GL_Light_Directional_A->text();
        if(p1.compare("Light_Directional_X") == 0)  namelist[j].value = this->GL_Light_Directional_X->text();
        if(p1.compare("Light_Directional_Y") == 0)  namelist[j].value = this->GL_Light_Directional_Y->text();
        if(p1.compare("Light_Directional_Z") == 0)  namelist[j].value = this->GL_Light_Directional_Z->text();

        /*if(p1.compare("Surface_Draw") == 0)  namelist[j].value.setNum(this->GL_Surface_Draw->isChecked());
        if(p1.compare("Surface_Micro_Elevation_Scale") == 0)  namelist[j].value = QString::number(this->GL_Surface_Micro_Elevation_Scale->value());
        if(p1.compare("Surface_Mipmap_Distance_1") == 0)  namelist[j].value = QString::number(this->GL_Surface_Mipmap_Distance_1->value());
        if(p1.compare("Surface_Mipmap_Distance_2") == 0)  namelist[j].value = QString::number(this->GL_Surface_Mipmap_Distance_2->value());
        if(p1.compare("Surface_Vegetated_Small_Color_R") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Small_R->text();
        if(p1.compare("Surface_Vegetated_Small_Color_G") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Small_G->text();
        if(p1.compare("Surface_Vegetated_Small_Color_B") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Small_B->text();
        if(p1.compare("Surface_Vegetated_Large_Color_R") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Large_R->text();
        if(p1.compare("Surface_Vegetated_Large_Color_G") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Large_G->text();
        if(p1.compare("Surface_Vegetated_Large_Color_B") == 0)  namelist[j].value = this->GL_Surface_Vegetated_Large_B->text();
        if(p1.compare("Surface_Vegetated_Bare_Color_R") == 0)  namelist[j].value = this->GL_Surface_Bare_R->text();
        if(p1.compare("Surface_Vegetated_Bare_Color_G") == 0)  namelist[j].value = this->GL_Surface_Bare_G->text();
        if(p1.compare("Surface_Vegetated_Bare_Color_B") == 0)  namelist[j].value = this->GL_Surface_Bare_B->text();
        if(p1.compare("Surface_Vegetated_Roads_Color_R") == 0)  namelist[j].value = this->GL_Surface_Road_R->text();
        if(p1.compare("Surface_Vegetated_Roads_Color_G") == 0)  namelist[j].value = this->GL_Surface_Road_G->text();
        if(p1.compare("Surface_Vegetated_Roads_Color_B") == 0)  namelist[j].value = this->GL_Surface_Road_B->text();
        if(p1.compare("Surface_Vegetated_Buildings_Color_R") == 0)  namelist[j].value = this->GL_Surface_Building_R->text();
        if(p1.compare("Surface_Vegetated_Buildings_Color_G") == 0)  namelist[j].value = this->GL_Surface_Building_G->text();
        if(p1.compare("Surface_Vegetated_Buildings_Color_B") == 0)  namelist[j].value = this->GL_Surface_Building_B->text();
        if(p1.compare("Surface_Vegetated_Erosion_Color_R") == 0)  namelist[j].value = this->GL_Surface_Erosion_Color_R->text();
        if(p1.compare("Surface_Vegetated_Erosion_Color_G") == 0)  namelist[j].value = this->GL_Surface_Erosion_Color_G->text();
        if(p1.compare("Surface_Vegetated_Erosion_Color_B") == 0)  namelist[j].value = this->GL_Surface_Erosion_Color_B->text();
        if(p1.compare("Surface_Vegetated_Erosion_Color_A") == 0)  namelist[j].value = this->GL_Surface_Erosion_Color_A->text();
        if(p1.compare("Surface_Vegetated_Deposition_Color_R") == 0)  namelist[j].value = this->GL_Surface_Deposition_Color_R->text();
        if(p1.compare("Surface_Vegetated_Deposition_Color_G") == 0)  namelist[j].value = this->GL_Surface_Deposition_Color_G->text();
        if(p1.compare("Surface_Vegetated_Deposition_Color_B") == 0)  namelist[j].value = this->GL_Surface_Deposition_Color_B->text();
        if(p1.compare("Surface_Vegetated_Deposition_Color_A") == 0)  namelist[j].value = this->GL_Surface_Deposition_Color_A->text();

        if(p1.compare("Water_Draw") == 0)  namelist[j].value.setNum(this->GL_Water_Draw->isChecked());
        if(p1.compare("Water_Reflectivity") == 0)  namelist[j].value = QString::number(this->GL_Water_Reflectivity->value());
        if(p1.compare("Water_Refractivity") == 0)  namelist[j].value = QString::number(this->GL_Water_Refractivity->value());
        if(p1.compare("Water_Velocity_Scale") == 0)  namelist[j].value = QString::number(this->GL_Water_Velocity_Scale->value());
        if(p1.compare("Water_Micro_Elevation_Scale") == 0)  namelist[j].value = QString::number(this->GL_Water_Micro_Elevation_Scale->value());
        if(p1.compare("Water_Transparancy") == 0)  namelist[j].value = QString::number(this->GL_Water_Transparancy->value());
        if(p1.compare("Water_Deep_Color_R") == 0)  namelist[j].value = this->GL_Water_Deep_R->text();
        if(p1.compare("Water_Deep_Color_G") == 0)  namelist[j].value = this->GL_Water_Deep_G->text();
        if(p1.compare("Water_Deep_Color_B") == 0)  namelist[j].value = this->GL_Water_Deep_B->text();
        if(p1.compare("Water_Deep_Color_A") == 0)  namelist[j].value = this->GL_Water_Deep_A->text();
        if(p1.compare("Water_Shallow_Color_R") == 0)  namelist[j].value = this->GL_Water_Shallow_R->text();
        if(p1.compare("Water_Shallow_Color_G") == 0)  namelist[j].value = this->GL_Water_Shallow_G->text();
        if(p1.compare("Water_Shallow_Color_B") == 0)  namelist[j].value = this->GL_Water_Shallow_B->text();
        if(p1.compare("Water_Shallow_Color_A") == 0)  namelist[j].value = this->GL_Water_Shallow_A->text();
        if(p1.compare("Water_Sediment_Color_R") == 0)  namelist[j].value = this->GL_Water_Sediment_R->text();
        if(p1.compare("Water_Sediment_Color_G") == 0)  namelist[j].value = this->GL_Water_Sediment_G->text();
        if(p1.compare("Water_Sediment_Color_B") == 0)  namelist[j].value = this->GL_Water_Sediment_B->text();
        if(p1.compare("Water_Sediment_Color_A") == 0)  namelist[j].value = this->GL_Water_Sediment_A->text();*/

        if(p1.compare("Clouds_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Clouds_Draw->isChecked());
        if(p1.compare("Rain_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Rain_Draw->isChecked());
        if(p1.compare("Roads_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Roads_Draw->isChecked());
        //if(p1.compare("Roads_Distance") == 0)  namelist[j].value = QString::number(this->GL_Objects_Roads_Distance->value());

        if(p1.compare("Buildings_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Buildings_Draw->isChecked());
        if(p1.compare("Buildings_Distance") == 0)  namelist[j].value = QString::number(this->GL_Objects_Buildings_Distance->value());

        if(p1.compare("Trees_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Trees_Draw->isChecked());
        if(p1.compare("Trees_Distance") == 0)  namelist[j].value = QString::number(this->GL_Objects_Trees_Distance->value());
        //if(p1.compare("Trees_Instances") == 0)  namelist[j].value = QString::number(this->GL_Objects_Trees_Instances->value());
        if(p1.compare("Trees_Increment") == 0)  namelist[j].value = QString::number(this->GL_Objects_Trees_Increment->value());

        if(p1.compare("Grass_Draw") == 0)  namelist[j].value.setNum(this->GL_Objects_Grass_Draw->isChecked());
        if(p1.compare("Grass_Distance") == 0)  namelist[j].value = QString::number(this->GL_Objects_Grass_Distance->value());
        //if(p1.compare("Grass_Instances") == 0)  namelist[j].value = QString::number(this->GL_Objects_Grass_Instances->value());
        if(p1.compare("Grass_Increment") == 0)  namelist[j].value = QString::number(this->GL_Objects_Grass_Increment->value());
        //if(p1.compare("Grass_Vertical_Scale") == 0)  namelist[j].value = QString::number(this->GL_Objects_Grass_Vertical_Scale->value());


        //   if (p1.compare("Max flood level")==0) namelist[j].value = E_maxFloodLevel->text();
        //   if (p1.compare("Min flood dt")==0) namelist[j].value = E_minFloodDt->text();

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

        //VJ 111120 water repellency
        if (p1.compare("Use Water Repellency")==0)      namelist[j].value.setNum((int)checkWaterRepellency->isChecked());
        if (p1.compare("Water Repellency A")==0)        namelist[j].value.setNum(E_waterRep_a->value(),'g',6);
        if (p1.compare("Water Repellency B")==0)        namelist[j].value.setNum(E_waterRep_b->value(),'g',6);
        if (p1.compare("Water Repellency C")==0)        namelist[j].value.setNum(E_waterRep_c->value(),'g',6);
        if (p1.compare("Water Repellency D")==0)        namelist[j].value.setNum(E_waterRep_d->value(),'g',6);

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
        if (p1.compare("CheckOutputMaps")==0)
        {
            outputcheck.clear();
            if ( 	checkBox_OutRunoff->isChecked()) outputcheck << "1"; else outputcheck << "0"; //0
            if ( 	  checkBox_OutConc->isChecked()) outputcheck << "1"; else outputcheck << "0"; //1
            if ( 	    checkBox_OutWH->isChecked()) outputcheck << "1"; else outputcheck << "0"; //2
            /* OBSOLETE if ( 	   checkBox_OutWHC->isChecked()) outputcheck << "1"; else */      //3
            outputcheck << "0";
            if (       checkBox_OutTC->isChecked())  outputcheck << "1"; else outputcheck << "0"; //4
            if ( 	   checkBox_OutDet->isChecked()) outputcheck << "1"; else outputcheck << "0"; //5
            if ( 	   checkBox_OutDep->isChecked()) outputcheck << "1"; else outputcheck << "0"; //6
            if ( 	     checkBox_OutV->isChecked()) outputcheck << "1"; else outputcheck << "0"; //7
            if ( 	   checkBox_OutInf->isChecked()) outputcheck << "1"; else outputcheck << "0"; //8
            if (checkBox_OutSurfStor->isChecked())   outputcheck << "1"; else outputcheck << "0"; //9
            if (checkBox_OutChanVol->isChecked())    outputcheck << "1"; else outputcheck << "0"; //10
            if (checkBox_OutTiledrain->isChecked())  outputcheck << "1"; else outputcheck << "0"; //11

            if (checkBox_OutHmx->isChecked())        outputcheck << "1"; else outputcheck << "0"; //12
            if (checkBox_OutQf->isChecked())         outputcheck << "1"; else outputcheck << "0"; //13
            if (checkBox_OutVf->isChecked())         outputcheck << "1"; else outputcheck << "0"; //14
            if (checkBox_OutHmxWH->isChecked())      outputcheck << "1"; else outputcheck << "0"; //15
            if (checkBox_OutSL->isChecked())         outputcheck << "1"; else outputcheck << "0"; //16
            if (checkBox_OutSed->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17

            if (checkBox_OutSafetyFactor->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutSlopeFailure->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutDebrisFlowHeight->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutDebrisFlowVelocity->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutFluidPhaseHeight->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutSolidPhaseHeight->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutEntrainment->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17
            if (checkBox_OutTimestep->isChecked())         outputcheck << "1"; else outputcheck << "0"; //17


            outputcheck << "0";

            // twenty-six places for now
            namelist[j].value = outputcheck.join(",");
        }
    }

    //get all actual mapnames from the mapList structure
    fillNamelistMapnames(true);

    currentDir = E_WorkDir;//->text();
    QDir::setCurrent(currentDir);
}
