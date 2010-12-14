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
	change_MapNameModel(EROSIONMAPS, 0, !checkNoErosion->isChecked());
	sedgroup->setEnabled(!checkNoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
	if (checkIncludeChannel->isChecked())
	{
		change_MapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
		change_MapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
	}
	change_MapNameModel(CHANNELSMAPS, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{
	if (checkChannelBaseflow->isChecked())
		checkChannelBaseflow->setChecked(false);
	change_MapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
	change_MapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
	if (checkChannelInfil->isChecked())
		checkChannelInfil->setChecked(false);
	change_MapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
	change_MapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
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
    
	change_MapNameModel(INFILTRATIONMAPS, 0, true);
	change_MapNameModel(INFILTRATIONMAPS, 10, false);//SW
	change_MapNameModel(INFILTRATIONMAPS, 11, false);//GA1
	change_MapNameModel(INFILTRATIONMAPS, 12, false);//GA2
	change_MapNameModel(INFILTRATIONMAPS, 13, false);//KS
	change_MapNameModel(INFILTRATIONMAPS, 14, false);//SP
    
	if (nr == 0)
	{
		change_MapNameModel(INFILTRATIONMAPS, 0, false);
	}
	else
	{
		change_MapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                            || checkInfilCompact->isChecked()
                            || checkInfilGrass->isChecked()
                            );
        
		change_MapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());
        
		if (nr == 1) change_MapNameModel(INFILTRATIONMAPS, 10, true);
		else
			if (nr == 2 || nr == 3) change_MapNameModel(INFILTRATIONMAPS, 11, true);
        else
            if (nr == 4)
                change_MapNameModel(INFILTRATIONMAPS, 13, true);
	}
}
//--------------------------------------------------------------------
void lisemqt::on_checkInfil2layer_clicked()
{
	if (E_InfiltrationMethod->currentIndex() == 2 ||
        E_InfiltrationMethod->currentIndex() == 3)
		change_MapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked());
    //	else
    //	{
    //		checkInfil2layer->setChecked(false);
    //	}
}//--------------------------------------------------------------------
void lisemqt::on_checkInfilCompact_clicked()
{
	if (E_InfiltrationMethod->currentIndex() > 0)
		change_MapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
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
	change_MapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked()||checkSedtrap->isChecked());
	buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	label_33->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	E_BulkDens->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	if (checkSedtrap->isChecked())
		checkSedtrap->setChecked(false);
    
}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{
	change_MapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked()||checkSedtrap->isChecked());
	buffergroup->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	label_33->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	E_BulkDens->setEnabled(checkBuffers->isChecked()||checkSedtrap->isChecked());
	if(checkBuffers->isChecked())
		checkBuffers->setChecked(false);
    
}
//--------------------------------------------------------------------
void lisemqt::on_checkSnowmelt_clicked()
{
	change_MapNameModel(SNOWMELTMAPS, 0, checkSnowmelt->isChecked());

    E_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
    label_5->setEnabled(checkSnowmelt->isChecked());
    toolButton_SnowmeltShow->setEnabled(checkSnowmelt->isChecked());
    toolButton_SnowmeltName->setEnabled(checkSnowmelt->isChecked());
}

void lisemqt::doCheckSnowmelt(bool check)
{
    if (!check && !checkRainfall->isChecked())
    {
		QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");
        checkSnowmelt->setChecked(true);
        check = true;
    }   

	change_MapNameModel(SNOWMELTMAPS, 0, check);

    E_SnowmeltName->setEnabled(check);
    label_5->setEnabled(check);
    toolButton_SnowmeltShow->setEnabled(check);
    toolButton_SnowmeltName->setEnabled(check);
}
//--------------------------------------------------------------------
void lisemqt::doCheckRainfall(bool check)
{
    if (!check && !checkSnowmelt->isChecked())
    {
		QMessageBox::warning(this,"openLISEM","Must have rainfall, snowmelt or both");
        checkRainfall->setChecked(true);
        check = true;
    }   

	change_MapNameModel(RAINFALLMAPS, 0, check);

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
		change_MapNameModel(i, 0, false);
    // set all to false
    
    // PROCESS IN REVERSE ORDER
    
	change_MapNameModel(CHANNELSMAPS, 11, checkChannelInfil->isChecked());
	change_MapNameModel(CHANNELSMAPS, 12, checkChannelBaseflow->isChecked());
	change_MapNameModel(CHANNELSMAPS, 10, checkIncludeChannel->isChecked());
	change_MapNameModel(SNOWMELTMAPS, 0, checkSnowmelt->isChecked());
	change_MapNameModel(BUFFERSMAPS, 0, checkBuffers->isChecked());
	change_MapNameModel(BUFFERSMAPS, 0, checkSedtrap->isChecked());
    
	int nr = E_InfiltrationMethod->currentIndex();
	checkInfil2layer->setEnabled(bool(nr ==2 || nr == 3));
    
	if (nr == 0)
		change_MapNameModel(INFILTRATIONMAPS, 0, false);
	else
	{
		change_MapNameModel(INFILTRATIONMAPS, 0, true); //all
        
		change_MapNameModel(INFILTRATIONMAPS, 10, false);//SW
		change_MapNameModel(INFILTRATIONMAPS, 11, false);//GA1
		change_MapNameModel(INFILTRATIONMAPS, 12, false);//GA2
		change_MapNameModel(INFILTRATIONMAPS, 13, false);//KS
		change_MapNameModel(INFILTRATIONMAPS, 14, false);//SP
        
		change_MapNameModel(INFILTRATIONMAPS, 14, checkInfilCrust->isChecked()
                            || checkInfilCompact->isChecked()
                            || checkInfilGrass->isChecked()
                            );
        
		change_MapNameModel(INFILTRATIONMAPS, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());
        
		if (nr == 1) change_MapNameModel(INFILTRATIONMAPS, 10, true);
		else
			if (nr == 2 || nr == 3) change_MapNameModel(INFILTRATIONMAPS, 11, true);
        else
            change_MapNameModel(INFILTRATIONMAPS, 12, false);
	}
	change_MapNameModel(EROSIONMAPS, 0, !checkNoErosion->isChecked());
	change_MapNameModel(SURFACEMAPS, 0, true);
	change_MapNameModel(LANDUSEMAPS, 0, true);
	change_MapNameModel(CATCHMENTMAPS, 0, true);
	change_MapNameModel(RAINFALLMAPS, 0, checkRainfall->isChecked());
    
	buffergroup->setEnabled(checkBuffers->isChecked());
	buffergroup->setEnabled(checkSedtrap->isChecked());
	sedgroup->setEnabled(!checkNoErosion->isChecked());
    
	groupBox_SwatreOptions->setEnabled(E_InfiltrationMethod->currentIndex() == 1);
	checkExpandActive->setChecked(false);
	treeView->collapseAll();
}
//--------------------------------------------------------------------
