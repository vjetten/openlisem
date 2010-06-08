/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/



/*
 * function that determine reactions of the map tree structure when
 * the user checks options in the interface
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"


//--------------------------------------------------------------------
void lisemqt::on_checkNoErosion_clicked()
{
	change_MapNameModel(3, 0, !checkNoErosion->isChecked());
	sedgroup->setEnabled(!checkNoErosion->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkIncludeChannel_clicked()
{
	if (checkIncludeChannel->isChecked())
	{
		change_MapNameModel(5, 11, checkChannelInfil->isChecked());
		change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
	}
	change_MapNameModel(5, 10, checkIncludeChannel->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelInfil_clicked()
{
	if (checkChannelBaseflow->isChecked())
		checkChannelBaseflow->setChecked(false);
	change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
	change_MapNameModel(5, 11, checkChannelInfil->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkChannelBaseflow_clicked()
{
	if (checkChannelInfil->isChecked())
		checkChannelInfil->setChecked(false);
	change_MapNameModel(5, 11, checkChannelInfil->isChecked());
	change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
}
//--------------------------------------------------------------------
//2nd number is number of rows at a level. e.g. green and ampt starts at
// after swatre, swatre has 11 rows (maps), starting at 0, so G&A starts at 11
//
void lisemqt::on_E_InfiltrationMethod_currentIndexChanged(int)
{
	int nr = E_InfiltrationMethod->currentIndex();
	checkInfil2layer->setEnabled(bool(nr ==2 || nr == 3));
	groupBox_SwatreOptions->setEnabled(nr == 1);

	uiInfilMethod = nr;

	change_MapNameModel(4, 0, true);
	change_MapNameModel(4, 10, false);//SW
	change_MapNameModel(4, 11, false);//GA1
	change_MapNameModel(4, 12, false);//GA2
	change_MapNameModel(4, 13, false);//KS
	change_MapNameModel(4, 14, false);//SP

	if (nr == 0)
	{
		change_MapNameModel(4, 0, false);
	}
	else
	{
		change_MapNameModel(4, 14, checkInfilCrust->isChecked()
				|| checkInfilCompact->isChecked()
				|| checkInfilGrass->isChecked()
		);

		change_MapNameModel(4, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());

		if (nr == 1) change_MapNameModel(4, 10, true);
		else
			if (nr == 2 || nr == 3) change_MapNameModel(4, 11, true);
			else
				if (nr == 4)
					change_MapNameModel(4, 13, true);
	}
}
//--------------------------------------------------------------------
void lisemqt::on_toolButton_SwatreTable_clicked()
{
	QString path;
	path = QFileDialog::getExistingDirectory(
			this,
			tr("Select directory with SWATRE tables"),
			QString::null,
			QFileDialog::ShowDirsOnly);

	E_SWATRETableDir->setText( path );
}

//--------------------------------------------------------------------
void lisemqt::on_checkInfil2layer_clicked()
{
	if (E_InfiltrationMethod->currentIndex() == 2 ||
			E_InfiltrationMethod->currentIndex() == 3)
		change_MapNameModel(4, 12, checkInfil2layer->isChecked());
//	else
//	{
//		checkInfil2layer->setChecked(false);
//	}
}//--------------------------------------------------------------------
void lisemqt::on_checkInfilCompact_clicked()
{
	if (E_InfiltrationMethod->currentIndex() > 0)
		change_MapNameModel(4, 14, checkInfilCrust->isChecked()
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
	change_MapNameModel(6, 0, checkBuffers->isChecked());
	buffergroup->setEnabled(checkBuffers->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkSedtrap_clicked()
{
	change_MapNameModel(6, 0, checkSedtrap->isChecked());
	buffergroup->setEnabled(checkSedtrap->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkSnowmelt_clicked()
{
	change_MapNameModel(7, 0, checkSnowmelt->isChecked());
}
//--------------------------------------------------------------------
void lisemqt::on_checkExpandActive_clicked()
{
	if (!checkExpandActive->isChecked())
		treeView->collapseAll();
	else
		for (int i = 0; i < MapNameModel->rowCount(); i++)
		{
			if (MapNameModel->getflag(i))
				treeView->expand(MapNameModel->index(i,0));
		}
}
//--------------------------------------------------------------------
void lisemqt::RunAllChecks()
{
	for (int i = 0; i < 12; i++)
		change_MapNameModel(i, 0, false);

	change_MapNameModel(5, 11, checkChannelInfil->isChecked());
	change_MapNameModel(5, 12, checkChannelBaseflow->isChecked());
	change_MapNameModel(5, 10, checkIncludeChannel->isChecked());
	change_MapNameModel(7, 0, checkSnowmelt->isChecked());
	change_MapNameModel(6, 0, checkBuffers->isChecked());

	int nr = E_InfiltrationMethod->currentIndex();
	checkInfil2layer->setEnabled(bool(nr ==2 || nr == 3));

	if (nr == 0)
		change_MapNameModel(4, 0, false);
	else
	{
		change_MapNameModel(4, 0, true); //all

		change_MapNameModel(4, 10, false);//SW
		change_MapNameModel(4, 11, false);//GA1
		change_MapNameModel(4, 12, false);//GA2
		change_MapNameModel(4, 13, false);//KS
		change_MapNameModel(4, 14, false);//SP

		change_MapNameModel(4, 14, checkInfilCrust->isChecked()
				|| checkInfilCompact->isChecked()
				|| checkInfilGrass->isChecked()
		);

		change_MapNameModel(4, 12, checkInfil2layer->isChecked() && checkInfil2layer->isEnabled());

		if (nr == 1) change_MapNameModel(4, 10, true);
		else
			if (nr == 2 || nr == 3) change_MapNameModel(4, 11, true);
			else
				change_MapNameModel(4, 12, false);
	}
	change_MapNameModel(3, 0, !checkNoErosion->isChecked());
	change_MapNameModel(2, 0, true);
	change_MapNameModel(1, 0, true);
	change_MapNameModel(0, 0, true);

	groupBox_SwatreOptions->setEnabled(E_InfiltrationMethod->currentIndex() == 1);
	checkExpandActive->setChecked(false);
	treeView->collapseAll();
}
//--------------------------------------------------------------------
