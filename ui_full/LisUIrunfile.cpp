/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
 * functions to read and parse the runfile
 */


#include "lisemqt.h"
#include "model.h"
#include "global.h"



//---------------------------------------------------------------------------
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

	for (int i = 0; i < NUMNAMES; i++)
	{
		namelist[i].name.clear();
		namelist[i].value.clear();
	}
	nrnamelist = 0;

	while (!fin.atEnd())
	{
		QString S = fin.readLine();
		if (!S.trimmed().isEmpty())
		{
			if (S.contains("="))
			{
				QStringList SL = S.split(QRegExp("="));
				namelist[nrnamelist].name = SL[0].trimmed();
				namelist[nrnamelist].value = SL[1].trimmed();
				nrnamelist++;
			}
			else
			{
				namelist[nrnamelist].name = S;
				namelist[nrnamelist].value = "";
				nrnamelist++;
			}
		}
	}
}
//---------------------------------------------------------------------------
void lisemqt::ParseInputData()
{
	int j=0;
	for (j = 0; j < nrnamelist; j++)
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
		//if (p1.compare("All water and sediment to outlet")==0) checkAllinChannel->setChecked(check);

		if (p1.compare("Include snowmelt")==0)               checkSnowmelt->setChecked(check);
		if (p1.compare("Alternative flow detachment")==0)    checkAltErosion->setChecked(check);
		if (p1.compare("Simple depression storage")==0)      checkSimpleDepression->setChecked(check);
		if (p1.compare("Hard Surfaces")==0)                  checkHardsurface->setChecked(check);
		if (p1.compare("Include buffers")==0)                checkBuffers->setChecked(check);
		if (p1.compare("Include Sediment traps")==0)         checkSedtrap->setChecked(check);
		if (p1.compare("Include wheeltracks")==0)            checkInfilCompact->setChecked(check);
		if (p1.compare("Include grass strips")==0)           checkInfilGrass->setChecked(check);
		if (p1.compare("Include crusts")==0)                 checkInfilCrust->setChecked(check);
		if (p1.compare("Impermeable sublayer")==0)           checkImpermeable->setChecked(check);
		if (p1.compare("Matric head files")==0)              checkDumphead->setChecked(check);
		if (p1.compare("Geometric mean Ksat")==0)            checkGeometricMean->setChecked(check);
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
		if (p1.compare("Use canopy storage map")==0)   	     checkInterceptionLAI->setChecked(!check);
		if (p1.compare("Sediment bulk density")==0)          E_BulkDens->setText(p);

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
			}
		}
		if (p1.compare("Infil Method")==0)
		{
			uiInfilMethod = 0;
			switch(iii){
				case INFIL_SWATRE : uiInfilMethod = 1; break;
				case INFIL_GREENAMPT : uiInfilMethod = 2; break;
				case INFIL_SMITH : uiInfilMethod = 3; break;
				case INFIL_GREENAMPT2 : uiInfilMethod = 2; checkInfil2layer->setChecked(true); break;
				case INFIL_KSAT : uiInfilMethod = 4; break;
			}
			E_InfiltrationMethod->setCurrentIndex(uiInfilMethod);
		}

		if (p1.compare("Ksat calibration")==0)         E_CalibrateKsat->setValue(namelist[j].value.toDouble());
		if (p1.compare("N calibration")==0)            E_CalibrateN->setValue(namelist[j].value.toDouble());
		if (p1.compare("Channel Ksat calibration")==0) E_CalibrateChKsat->setValue(namelist[j].value.toDouble());
		if (p1.compare("Channel N calibration")==0)    E_CalibrateChN->setValue(namelist[j].value.toDouble());
		if (p1.compare("Splash Delivery Ratio")==0)    E_SplashDelibery->setValue(namelist[j].value.toDouble());
		if (p1.compare("Stemflow fraction")==0)        E_StemflowFraction->setValue(namelist[j].value.toDouble());


		if (p1.compare("Output interval")==0)   printinterval->setValue(iii);

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

	// get main directories
	for (j = 0; j < nrnamelist; j++)
	{
		QString p1 = namelist[j].name;
		QString p = namelist[j].value;

		if (p1.compare("Begin time")==0) E_BeginTime->setText(p);
		if (p1.compare("End time")==0) E_EndTime->setText(p);
		if (p1.compare("Timestep")==0) E_Timestep->setText(p);

		// input ourput dirs and file names
		if (p1.compare("Map Directory")==0) E_MapDir->setText(CheckDir(p, p1));
		if (p1.compare("Result Directory")==0) E_ResultDir->setText(CheckDir(p, p1));

		//   if (p1.compare("Table Directory")==0) tableDir = CheckDir(p, p1);
		// move to swatre later when infiltration method is known!
		if (p1.compare("Main results file")==0) E_MainTotals->setText(p);
		if (p1.compare("Filename point output")==0) E_PointResults->setText(p);
		// resultDir is added in report operation

		if (p1.compare("Rainfall Directory")==0) rainFileDir = CheckDir(p, p1);
		if (p1.compare("Rainfall file")==0)
		{
			E_RainfallName->setText(p);
			RainFileName = rainFileDir + E_RainfallName->text();
		}

		if (p1.compare("Erosion map")==0) E_DetachmentMap->setText(p);
		if (p1.compare("Deposition map")==0) E_DepositionMap->setText(p);
		if (p1.compare("Soilloss map")==0) E_SoillossMap->setText(p);
		// resultDir is added in report operation
		if (checkSnowmelt->isChecked())
		{
			if (p1.compare("Snowmelt Directory")==0) snowmeltFileDir = CheckDir(p, p1);
			if (p1.compare("Snowmelt file")==0)
			{
				E_SnowmeltName->setText(p);
				SnowmeltFileName = snowmeltFileDir + E_SnowmeltName->text();
			}
		}
	}

	// get input names
	for (j = 0; j < nrnamelist; j++)
	{
		for (int i = 0; i < DEFmaps.size(); i++)
		{
			QStringList S = DEFmaps.at(i).split(";");
			if (S.contains(namelist[j].name))
			{
				QFileInfo fil(namelist[j].value);
				S.replace(2, fil.fileName() );
				DEFmaps.replace(i, S.join(";") );
			}
		}
	}

	RunAllChecks();
	//buffergroup->setEnabled(checkBuffers->isChecked());
	//buffergroup->setEnabled(checkSedtrap->isChecked());
	//sedgroup->setEnabled(!checkNoErosion->isChecked());

}
//---------------------------------------------------------------------------
QString lisemqt::CheckDir(QString p, QString p1)
{
	p.replace("/","\\");
	//	if (!p.endsWith("/"))
	if (!p.endsWith("\\"))
		p = p + "\\";
	if (QDir(p).exists())
		return p;
	else
		return "";
}
//---------------------------------------------------------------------------
// change runfile strings with current interface options
void lisemqt::UpdateModelData()
{
	for (int j = 0; j < nrnamelist; j++)
	{
		QString p1 = namelist[j].name;
		QString p = namelist[j].value;
		if (p1.compare("No Erosion simulation")==0) 			  namelist[j].value.setNum((int)checkNoErosion->isChecked());
		if (p1.compare("Include main channels")==0) 			  namelist[j].value.setNum((int)checkIncludeChannel->isChecked());
		if (p1.compare("Include channel infil")==0)          namelist[j].value.setNum((int)checkChannelInfil->isChecked());
		if (p1.compare("Include channel baseflow")==0)       namelist[j].value.setNum((int)checkChannelBaseflow->isChecked());
		if (p1.compare("Include snowmelt")==0)               namelist[j].value.setNum((int)checkSnowmelt->isChecked());
		if (p1.compare("Alternative flow detachment")==0)    namelist[j].value.setNum((int)checkAltErosion->isChecked());
		if (p1.compare("Simple depression storage")==0)      namelist[j].value.setNum((int)checkSimpleDepression->isChecked());
		if (p1.compare("Hard Surfaces")==0)                  namelist[j].value.setNum((int)checkHardsurface->isChecked());
		if (p1.compare("Include buffers")==0)                namelist[j].value.setNum((int)checkBuffers->isChecked());
		if (p1.compare("Include Sediment traps")==0)         namelist[j].value.setNum((int)checkSedtrap->isChecked());
		if (p1.compare("Include wheeltracks")==0)            namelist[j].value.setNum((int)checkInfilCompact->isChecked());
		if (p1.compare("Include grass strips")==0)           namelist[j].value.setNum((int)checkInfilGrass->isChecked());
		if (p1.compare("Include crusts")==0)                 namelist[j].value.setNum((int)checkInfilCrust->isChecked());
		if (p1.compare("Impermeable sublayer")==0)           namelist[j].value.setNum((int)checkImpermeable->isChecked());
		if (p1.compare("Matric head files")==0)              namelist[j].value.setNum((int)checkDumphead->isChecked());
		if (p1.compare("Geometric mean Ksat")==0)            namelist[j].value.setNum((int)checkGeometricMean->isChecked());
		if (p1.compare("Timeseries as PCRaster")==0)         namelist[j].value.setNum((int)checkWritePCRnames->isChecked());
		if (p1.compare("Timeplot as PCRaster")==0)           namelist[j].value.setNum((int)checkWritePCRtimeplot->isChecked());
		if (p1.compare("Regular runoff output")==0)          namelist[j].value.setNum((int)checkOutputTimeStep->isChecked());
		if (p1.compare("User defined output")==0)            namelist[j].value.setNum((int)checkOutputTimeUser->isChecked());
		if (p1.compare("No erosion at outlet")==0)           namelist[j].value.setNum((int)checkNoErosionOutlet->isChecked());
		if (p1.compare("Report point output separate")==0)   namelist[j].value.setNum((int)checkSeparateOutput->isChecked());
		if (p1.compare("Report point output for SOBEK")==0)  namelist[j].value.setNum((int)checkSOBEKOutput->isChecked());
		if (p1.compare("SOBEK date string")==0)              namelist[j].value = SOBEKdatestring->text();
		if (p1.compare("Sediment bulk density")==0)          namelist[j].value = E_BulkDens->text();
		if (p1.compare("Use canopy storage map")==0)   	     namelist[j].value.setNum((int)!checkInterceptionLAI->isChecked());
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
			namelist[j].value.setNum(i);
		}

		if (p1.compare("Begin time")==0) namelist[j].value = E_BeginTime->text();
		if (p1.compare("End time")==0)   namelist[j].value = E_EndTime->text();
		if (p1.compare("Timestep")==0)   namelist[j].value = E_Timestep->text();
		if (p1.compare("Map Directory")==0)    namelist[j].value = E_MapDir->text();
		if (p1.compare("Result Directory")==0) namelist[j].value = E_ResultDir->text();
		if (p1.compare("Main results file")==0) namelist[j].value = E_MainTotals->text();
		if (p1.compare("Filename point output")==0) namelist[j].value = E_PointResults->text();
		if (p1.compare("Rainfall Directory")==0) namelist[j].value = rainFileDir;
		if (p1.compare("Rainfall file")==0) namelist[j].value = E_RainfallName->text();
		if (p1.compare("Erosion map")==0) namelist[j].value = E_DetachmentMap->text();
		if (p1.compare("Deposition map")==0) namelist[j].value = E_DepositionMap->text();
		if (p1.compare("Soilloss map")==0) namelist[j].value = E_SoillossMap->text();
		if (p1.compare("Snowmelt Directory")==0) namelist[j].value = snowmeltFileDir;
		if (p1.compare("Snowmelt file")==0) namelist[j].value = E_SnowmeltName->text();
		if (p1.compare("Ksat calibration")==0) namelist[j].value = E_CalibrateKsat->text();
		if (p1.compare("N calibration")==0) namelist[j].value = E_CalibrateN->text();
		if (p1.compare("Channel Ksat calibration")==0) namelist[j].value = E_CalibrateChKsat->text();
		if (p1.compare("Channel N calibration")==0) namelist[j].value = E_CalibrateChN->text();
		if (p1.compare("Splash Delivery Ratio")==0) namelist[j].value = E_SplashDelibery->text();
		if (p1.compare("Stemflow fraction")==0) namelist[j].value = E_StemflowFraction->text();

		if (p1.compare("Infil Method")==0)
		{
			switch(uiInfilMethod)
			{
				case 0 : namelist[j].value.setNum(INFIL_NONE); break;
				case 1 : namelist[j].value.setNum(INFIL_SWATRE);break;
				case 2 : if(checkInfil2layer->isChecked()) namelist[j].value.setNum(INFIL_GREENAMPT2);
							else namelist[j].value.setNum(INFIL_GREENAMPT);break;
				case 3: namelist[j].value.setNum(INFIL_SMITH); break;
				case 4: namelist[j].value.setNum(INFIL_KSAT); break;
			}
		}
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
	}
}
//---------------------------------------------------------------------------