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
//#include "model.h"
#include "global.h"
/*
RAINFALLMAPS
CATCHMENTMAPS
LANDUSEMAPS
SURFACEMAPS
EROSIONMAPS
INFILTRATIONMAPS
CHANNELMAPS
BUFFERSMAPS
SNOWMELTMAPS

NUTRIENTSMAPS

*/


//--------------------------------------------------------------------
void lisemqt::on_checkFlowBarriers_clicked()
{
    checkMapNameModel(BARRIERMAPS, 0, checkFlowBarriers->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkDoErosion_clicked()
{
    checkMapNameModel(EROSIONMAPS, 0, checkDoErosion->isChecked());

    setErosionMapOutput(checkDoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeTiledrains_clicked()
{
    checkMapNameModel(TILEDRAINMAPS, 0, checkIncludeTiledrains->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkHouses_clicked()
{
    checkMapNameModel(HOUSESMAPS, 0, checkHouses->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
//    if (checkIncludeChannel->isChecked())
//    {
//        checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
//        checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
//        checkMapNameModel(CHANNELMAPS, 13, checkChannelFlood->isChecked());
//    }

    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());

    checkChannelInfil->setEnabled(checkIncludeChannel->isChecked());
    checkChannelBaseflow->setEnabled(checkIncludeChannel->isChecked());
    checkChannelCulverts->setEnabled(checkIncludeChannel->isChecked());
    //checkChannelFlood->setEnabled(checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{

//    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
//    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
//    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
    on_checkIncludeChannel_clicked();

    if (checkChannelBaseflow->isChecked())
        checkChannelBaseflow->setChecked(false);
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
//    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
//    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
//    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
    on_checkIncludeChannel_clicked();

    label_195->setEnabled(checkChannelBaseflow->isChecked());
    label_baseflowtot->setEnabled(checkChannelBaseflow->isChecked());

    if (checkChannelInfil->isChecked())
        checkChannelInfil->setChecked(false);
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
    checkBox_OutTheta->setEnabled(bool(nr == 2 || nr == 3));

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
    //checkPercolation->setDisabled(bool(nr == 1));
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
                          );
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilCrust_clicked()
{
    on_checkInfilCompact_clicked();
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfilGrass_clicked()
{
    checkMapNameModel(CONSERVATIONMAPS, 0, checkSedtrap->isChecked()||checkInfilGrass->isChecked());
    E_GrassStripN->setEnabled(checkInfilGrass->isChecked());
    label_15->setEnabled(checkInfilGrass->isChecked());
    label_98->setEnabled(checkInfilGrass->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{
    bool yes = checkSedtrap->isChecked();
    checkMapNameModel(CONSERVATIONMAPS, 0, checkSedtrap->isChecked()||checkInfilGrass->isChecked());
    E_SedTrapN->setEnabled(yes);
    E_BulkDens->setEnabled(yes);
     label_189->setEnabled(yes);
     label_191->setEnabled(yes);
     label_104->setEnabled(yes);
      label_33->setEnabled(yes);

    // groupBoxConservation->setEnabled(checkSedtrap->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkMaterialDepth_clicked()
{
    bool yes = checkMaterialDepth->isChecked();
    E_DepositedCohesion->setEnabled(yes);
    E_BulkDens2->setEnabled(yes);
    label_130->setEnabled(yes);
    //label_190->setEnabled(yes);
    //matdepthparams->setEnabled(yes);
    label_105->setEnabled(yes);
}
//--------------------------------------------------------------------
//void lisemqt::on_E_BulkDens2_editingFinished()
//{
//     E_BulkDens->setText(E_BulkDens2->text());
//     qDebug() << "dens";

//}
////--------------------------------------------------------------------
//void lisemqt::on_E_BulkDens_editingFinished()
//{
//     E_BulkDens2->setText(E_BulkDens->text());
//     qDebug() << "dens2";
//}

//--------------------------------------------------------------------
//void lisemqt::on_checkSnowmelt_clicked()
//{
//    checkMapNameModel(RAINFALLMAPS, 0, checkSnowmelt->isChecked());

//    E_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
//    label_5->setEnabled(checkSnowmelt->isChecked());
//    toolButton_SnowmeltShow->setEnabled(checkSnowmelt->isChecked());
//    toolButton_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
//}
////--------------------------------------------------------------------
//void lisemqt::doCheckSnowmelt(bool check)
//{
//    checkMapNameModel(RAINFALLMAPS, 0, check);

//    E_SnowmeltName->setEnabled(check);
//    label_5->setEnabled(check);
//    toolButton_SnowmeltShow->setEnabled(check);
//    toolButton_SnowmeltName->setEnabled(check);
//}
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
   // checkMapNameModel(NUTRIENTSMAPS, 11, check);
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

 //   checkMapNameModel(CHANNELFLOODMAPS, 0, checkChannelFlood->isChecked());
    checkMapNameModel(CHANNELMAPS, 12, checkChannelBaseflow->isChecked());
    checkMapNameModel(CHANNELMAPS, 11, checkChannelInfil->isChecked());
    checkMapNameModel(CHANNELMAPS, 10, checkIncludeChannel->isChecked());
   // checkMapNameModel(RAINFALLMAPS, 0, checkSnowmelt->isChecked());

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
    checkMapNameModel(RAINFALLMAPS, 0, true);

    checkMapNameModel(TILEDRAINMAPS, 0, checkIncludeTiledrains->isChecked());
    //houses
    checkMapNameModel(HOUSESMAPS, 0, checkHouses->isChecked());

   // checkMapNameModel(NUTRIENTSMAPS, 10, checkPesticides->isChecked());

    checkExpandActive->setChecked(false);
    treeView->collapseAll();
}
