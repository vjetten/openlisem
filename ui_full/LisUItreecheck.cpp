/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/



/*
 * function that determine reactions of the map tree structure when
 * the user checks options in the interface
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"
/*
RAINFALLMAPS
CATCHMENTMAPS 
LANDUSEMAPS   
SURFACEMAPS   
EROSIONMAPS   
INFILTRATIONMAPS 
CHANNELSMAPS 
BUFFERSMAPS  
SNOWMELTMAPS 
WHEELTRACKSMAPS
TEXTUREMAPS   
NUTRIENTSMAPS 
GULLIESMAPS  
*/

//--------------------------------------------------------------------
void lisemqt::on_checkNoErosion_clicked()
{
   checkMapNameModel(EROSIONMAPS, 0, !checkNoErosion->isChecked());
	sedgroup->setEnabled(!checkNoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
	if (checkIncludeChannel->isChecked())
	{
      checkMapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
      checkMapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
	}
   checkMapNameModel(CHANNELSMAPS, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeTiledrains_clicked()
{
   checkMapNameModel(TILEDRAINMAPS, 0, checkIncludeTiledrains->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{
	if (checkChannelBaseflow->isChecked())
		checkChannelBaseflow->setChecked(false);
   checkMapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
   checkMapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
	if (checkChannelInfil->isChecked())
		checkChannelInfil->setChecked(false);
   checkMapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
   checkMapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
   QMessageBox::warning(this,"openLISEM",QString("NOT IMPLEMENTED YET"));
   checkChannelBaseflow->setChecked(false);
}
//--------------------------------------------------------------------
//2nd number is number of rows at a level. e.g. green and ampt starts at
// after swatre, swatre has 11 rows (maps), starting at 0, so G&A starts at 11
//
void lisemqt::on_E_InfiltrationMethod_currentIndexChanged(int inr)
{
	int nr = inr;//E_InfiltrationMethod->currentIndex();
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
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfil2layer_clicked()
{
	if (E_InfiltrationMethod->currentIndex() == 2 ||
        E_InfiltrationMethod->currentIndex() == 3)
      checkMapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked());
    //	else
    //	{
    //		checkInfil2layer->setChecked(false);
    //	}
}//--------------------------------------------------------------------
void lisemqt::on_checkInfilCompact_clicked()
{
	if (E_InfiltrationMethod->currentIndex() > 0)
      checkMapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                            || checkInfilCompact->isChecked()
                            || checkInfilGrass->isChecked());
    //	else
    //	{
    //		checkInfilCompact->setChecked(false);
    //	}
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
   checkMapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked()||checkSedtrap->isChecked());
	buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	label_33->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	E_BulkDens->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	if (checkSedtrap->isChecked())
		checkSedtrap->setChecked(false);
    
}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{
   checkMapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked()||checkSedtrap->isChecked());
	buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	label_33->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	E_BulkDens->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	if(checkBuffers->isChecked())
		checkBuffers->setChecked(false);
    
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
//    if (!check && !checkRainfall->isChecked())
//    {
//      QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");
//        checkSnowmelt->setChecked(true);
//        check = true;
//    }

   checkMapNameModel(SNOWMELTMAPS, 0, check);

    E_SnowmeltName->setEnabled(check);
    label_5->setEnabled(check);
    toolButton_SnowmeltShow->setEnabled(check);
    toolButton_SnowmeltName->setEnabled(check);
}
//--------------------------------------------------------------------
void lisemqt::doCheckRainfall(bool check)
{
//    if (!check && !checkSnowmelt->isChecked())
//    {
//		QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");
//        checkRainfall->setChecked(true);
//        check = true;
//    }

   checkMapNameModel(RAINFALLMAPS, 0, check);

    E_RainfallName->setEnabled(check);
    label_4->setEnabled(check);
    toolButton_RainfallShow->setEnabled(check);
    toolButton_RainfallName->setEnabled(check);
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
    
   checkMapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
   checkMapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
   checkMapNameModel(CHANNELSMAPS, 10, checkIncludeChannel->isChecked());
   checkMapNameModel(SNOWMELTMAPS, 0, checkSnowmelt->isChecked());
   checkMapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked());
   checkMapNameModel(BUFFERSMAPS, 0, checkSedtrap->isChecked());
    
	int nr = E_InfiltrationMethod->currentIndex();
	checkInfil2layer->setEnabled(bool(nr ==2 || nr == 3));
    
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
   checkMapNameModel(EROSIONMAPS, 0, !checkNoErosion->isChecked());
   checkMapNameModel(SURFACEMAPS, 0, true);
   checkMapNameModel(LANDUSEMAPS, 0, true);
   checkMapNameModel(CATCHMENTMAPS, 0, true);
   checkMapNameModel(RAINFALLMAPS, 0, checkRainfall->isChecked());

   checkChannelInfil->setEnabled(checkIncludeChannel->isChecked());
   checkChannelBaseflow->setEnabled(checkIncludeChannel->isChecked());
   // VJ 110110

	buffergroup->setEnabled(checkBuffers->isChecked());
	buffergroup->setEnabled(checkSedtrap->isChecked());
	sedgroup->setEnabled(!checkNoErosion->isChecked());
    
	groupBox_SwatreOptions->setEnabled(E_InfiltrationMethod->currentIndex() == 1);
	checkExpandActive->setChecked(false);
	treeView->collapseAll();
}
//--------------------------------------------------------------------
