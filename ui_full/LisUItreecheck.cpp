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
  \file LisUItreecheck.cpp
  \brief map tree interaction

 * function that determine reactions of the map tree structure when
 * the user checks main options in the interface: disable or enbale braches
 * e.g. channel maps become enabled when this option is chosen
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"

/*
#define RAINFALLMAPS 0
#define CATCHMENTMAPS 1
#define LANDUSEMAPS 2
#define SURFACEMAPS 3
#define EROSIONMAPS 4
#define SLOPESTABILITYMAPS 5
#define ENTRAINMENTMAPS 6
#define INFILTRATIONMAPS 7
#define CHANNELMAPS 8
#define SURFACEFLOWMAPS 9
#define SNOWMELTMAPS 10
#define HOUSESMAPS 11

#define CHANNELPROPERTIES 1
#define CHANNELINFIL 2
#define CHANNELBASEFLOW 3

#define TOPLAYERSLOPESTABILITYMAPS 1
#define BOTTOMLAYERSLOPESTABILITYMAPS 2
#define SEISMICSLOPESTABILITYMAPS 2

#define SWATREMAPS 1
#define GREENANDAMPTLAYER1MAPS 2
#define GREENANDAMPTLAYER2MAPS 3
#define KSATSUBSTRACTIONMAPS 4

#define FLOWBARRIERS 1
#define FLOWLIMITING 2
#define FLOWINITIAL 3
#define FLOWFORCED 4
#define FLOWCELLBARRIERS 5
*/

void lisemqt::on_checkDoInitialFlow_clicked()
{
    checkMapNameModel(SURFACEFLOWMAPS, 10 +FLOWINITIAL, checkBox_UFInitial->isChecked());
}
void lisemqt::on_checkDoForcedFlow_clicked()
{
    checkMapNameModel(SURFACEFLOWMAPS, 10 +FLOWFORCED, checkBox_UFForced->isChecked());
}
void lisemqt::on_checkDoFlowLimiting_clicked()
{
    //checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());
}
void lisemqt::on_checkDoFlowBarriers_clicked()
{
    //checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());
}
void lisemqt::on_checkDoFlowCellBarriers_clicked()
{
    //checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());
}
void lisemqt::on_checkDoSlopeStabilityTop_clicked()
{
    checkMapNameModel(SLOPESTABILITYMAPS, 10 +TOPLAYERSLOPESTABILITYMAPS, checkSlopeStability->isChecked());
}
void lisemqt::on_checkDoSlopeStabilityBottom_clicked()
{
    checkMapNameModel(SLOPESTABILITYMAPS, 10 +BOTTOMLAYERSLOPESTABILITYMAPS, checkBedRockLayer->isChecked());
}
void lisemqt::on_checkDoSlopeStabilitySeismic_clicked()
{
    checkMapNameModel(SLOPESTABILITYMAPS, 10 +SEISMICSLOPESTABILITYMAPS, checkSeismic->isChecked());
}
void lisemqt::on_checkDoEntrainment_clicked()
{
    checkMapNameModel(ENTRAINMENTMAPS, 0, checkSlopeStability->isChecked());
}

//--------------------------------------------------------------------
void lisemqt::on_checkDoErosion_clicked()
{
    checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
    if (checkIncludeChannel->isChecked())
    {
        checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
        checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
        checkMapNameModel(CHANNELMAPS, 13, checkChannelFlood->isChecked());
    }

    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());


    checkChannelInfil->setEnabled(checkIncludeChannel->isChecked());
    checkChannelBaseflow->setEnabled(checkIncludeChannel->isChecked());
    checkChannelFlood->setEnabled(checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeTiledrains_clicked()
{

}
//--------------------------------------------------------------------
void lisemqt::on_checkHouses_clicked()
{
    checkMapNameModel(HOUSESMAPS, 0, checkHouses->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{

    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());

    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelFlood_clicked()
{

}
//--------------------------------------------------------------------
//2nd number is number of rows at a level. e.g. green and ampt starts at
// after swatre, swatre has 11 rows (maps), starting at 0, so G&A starts at 11
//
void lisemqt::on_E_InfiltrationMethod_currentIndexChanged(int inr)
{
    int nr = std::max(0, inr);
    checkInfil2layer->setEnabled(bool(nr == 2 || nr == 3));
    groupBox_SwatreOptions->setEnabled(nr == 1);

    uiInfilMethod = nr;
    // set runfile var to infil nr

    checkMapNameModel(INFILTRATIONMAPS, 0, true);
    checkMapNameModel(INFILTRATIONMAPS, 10, false);//SW
    checkMapNameModel(INFILTRATIONMAPS, 11, false);//GA1
    checkMapNameModel(INFILTRATIONMAPS, 12, false);//GA2
    checkMapNameModel(INFILTRATIONMAPS, 13, false);//KS
    checkMapNameModel(INFILTRATIONMAPS, 14, false);//SP

    if (nr == 0)
    {
        checkMapNameModel(INFILTRATIONMAPS, 0, false);
    }
    else
    {
        checkMapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                          || checkInfilCompact->isChecked()
                          || checkInfilGrass->isChecked()
                          );

        checkMapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());

        if (nr == 1) checkMapNameModel(INFILTRATIONMAPS, 10, true);
        else
            if (nr == 2 || nr == 3) checkMapNameModel(INFILTRATIONMAPS, 11, true);
            else
                if (nr == 4)
                    checkMapNameModel(INFILTRATIONMAPS, 13, true);
    }
    groupBox_InfilOptions->setDisabled(bool(nr == 0));
    checkInfil2layer->setEnabled(bool(nr == 2 || nr == 3));

    groupBox_SwatreOptions->setEnabled(bool(nr == 1));
    checkPercolation->setDisabled(bool(nr == 1));
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfil2layer_clicked()
{
    if (E_InfiltrationMethod->currentIndex() == 2 ||
            E_InfiltrationMethod->currentIndex() == 3)
        checkMapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilCompact_clicked()
{
    if (E_InfiltrationMethod->currentIndex() > 0)
        checkMapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                          || checkInfilCompact->isChecked()
                          || checkInfilGrass->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilCrust_clicked()
{
    on_checkInfilCompact_clicked();
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilGrass_clicked()
{
    on_checkInfilCompact_clicked();
}
//--------------------------------------------------------------------
void lisemqt::on_checkBuffers_clicked()
{
    /*
    checkMapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked()||checkSedtrap->isChecked());
    buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
    label_33->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
    E_BulkDens->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
    if (checkSedtrap->isChecked())
        checkSedtrap->setChecked(false);
        */

}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{

}
//--------------------------------------------------------------------
void lisemqt::on_checkSnowmelt_clicked()
{
    checkMapNameModel(SNOWMELTMAPS, 0, checkSnowmelt->isChecked());

    E_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
    label_5->setEnabled(checkSnowmelt->isChecked());
    toolButton_SnowmeltShow->setEnabled(checkSnowmelt->isChecked());
    toolButton_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::doCheckSnowmelt(bool check)
{
    checkMapNameModel(SNOWMELTMAPS, 0, check);

    E_SnowmeltName->setEnabled(check);
    label_5->setEnabled(check);
    toolButton_SnowmeltShow->setEnabled(check);
    toolButton_SnowmeltName->setEnabled(check);
}
//--------------------------------------------------------------------
void lisemqt::doCheckRainfall(bool check)
{
    checkMapNameModel(RAINFALLMAPS, 0, check);

    E_RainfallName->setEnabled(check);
    label_4->setEnabled(check);
    toolButton_RainfallShow->setEnabled(check);
    toolButton_RainfallName->setEnabled(check);
}
//--------------------------------------------------------------------
void lisemqt::doCheckPesticides(bool check)
{

}
//--------------------------------------------------------------------
void lisemqt::on_checkExpandActive_clicked()
{
    if (!checkExpandActive->isChecked())
        treeView->collapseAll();
    else
        for (int i = 0; i < MapNameModel->rowCount(); i++)
        {
            if (MapNameModel->getFlag(i))
                treeView->expand(MapNameModel->index(i,0));
        }
}
//--------------------------------------------------------------------
void lisemqt::RunAllChecks()
{
    for (int i = 0; i < 12; i++)
        checkMapNameModel(i, 0, false);
    // set all to false

    // PROCESS IN REVERSE ORDER

    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
    checkMapNameModel(SNOWMELTMAPS, 0, checkSnowmelt->isChecked());

    int nr = E_InfiltrationMethod->currentIndex();
    if (nr == 0)
        checkMapNameModel(INFILTRATIONMAPS, 0, false);
    else
    {
        checkMapNameModel(INFILTRATIONMAPS, 0, true); //all

        checkMapNameModel(INFILTRATIONMAPS, 10, false);//SW
        checkMapNameModel(INFILTRATIONMAPS, 11, false);//GA1
        checkMapNameModel(INFILTRATIONMAPS, 12, false);//GA2
        checkMapNameModel(INFILTRATIONMAPS, 13, false);//KS
        checkMapNameModel(INFILTRATIONMAPS, 14, false);//SP

        checkMapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                          || checkInfilCompact->isChecked()
                          || checkInfilGrass->isChecked()
                          );

        checkMapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());

        if (nr == 1) checkMapNameModel(INFILTRATIONMAPS, 10, true);
        else
            if (nr == 2 || nr == 3) checkMapNameModel(INFILTRATIONMAPS, 11, true);
            else
                checkMapNameModel(INFILTRATIONMAPS, 12, false);
    }
    checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());
    checkMapNameModel(SURFACEMAPS, 0, true);
    checkMapNameModel(LANDUSEMAPS, 0, true);
    checkMapNameModel(CATCHMENTMAPS, 0, true);
    checkMapNameModel(SURFACEFLOWMAPS, 0, true);
    checkMapNameModel(RAINFALLMAPS, 0, checkRainfall->isChecked());

    //houses
    checkMapNameModel(HOUSESMAPS, 0, checkHouses->isChecked());
    checkMapNameModel(ENTRAINMENTMAPS, 0, checkSlopeStability->isChecked());
    checkMapNameModel(SURFACEFLOWMAPS, 10+FLOWINITIAL, checkBox_UFInitial->isChecked());
    checkMapNameModel(SURFACEFLOWMAPS, 10+FLOWFORCED, checkBox_UFForced->isChecked());
    checkMapNameModel(SLOPESTABILITYMAPS, 10 +SEISMICSLOPESTABILITYMAPS, checkSeismic->isChecked());
    checkExpandActive->setChecked(false);
    treeView->collapseAll();
}
