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
		//		if (p1.compare("Matric head files")==0)              checkDumphead->setChecked(check);
		if (p1.compare("Geometric mean Ksat")==0)            checkGeometric->setChecked(check);
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

		if (p1.compare("Table Directory")==0)
		{
			SwatreTableDir = CheckDir(p, p1);
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

	// get all map names
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

}
//---------------------------------------------------------------------------
void lisemqt::InsertVariable(QString q, QString p, QString p1)
{
	int j, pos = 0;
	for (j = 0; j < nrnamelist; j++)
	{
		if(namelist[j].name.compare(q)==0)
			break;
	}
	pos = j;
	for (j = nrnamelist; j > pos; j--)
	{
		namelist[j].name = namelist[j-1].name;
		namelist[j].value = namelist[j-1].value;
	}
	namelist[pos].name = p;
	namelist[pos].value = p1;
	nrnamelist++;
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
	DefaultRunFile();
	/*
	// add new variables here
	bool check = false;
	for (int j = 0; j < nrnamelist; j++)
	{
		if (namelist[j].name.compare("Table File")==0)
			check = true;
	}
	if (!check)
		InsertVariable(QString("Table Directory"), QString("Table file"), SwatreTableName);
*/

	for (int j = 0; j < nrdefnamelist; j++)
	{
		QString p1 = defnamelist[j].name;
		QString p;
		if (p1.compare("No Erosion simulation")==0) 			  defnamelist[j].value.setNum((int)checkNoErosion->isChecked());
		if (p1.compare("Include main channels")==0) 			  defnamelist[j].value.setNum((int)checkIncludeChannel->isChecked());
		if (p1.compare("Include channel infil")==0)          defnamelist[j].value.setNum((int)checkChannelInfil->isChecked());
		if (p1.compare("Include channel baseflow")==0)       defnamelist[j].value.setNum((int)checkChannelBaseflow->isChecked());
		if (p1.compare("Include snowmelt")==0)               defnamelist[j].value.setNum((int)checkSnowmelt->isChecked());
		if (p1.compare("Alternative flow detachment")==0)    defnamelist[j].value.setNum((int)checkAltErosion->isChecked());
		if (p1.compare("Simple depression storage")==0)      defnamelist[j].value.setNum((int)checkSimpleDepression->isChecked());
		if (p1.compare("Hard Surfaces")==0)                  defnamelist[j].value.setNum((int)checkHardsurface->isChecked());
		if (p1.compare("Include buffers")==0)                defnamelist[j].value.setNum((int)checkBuffers->isChecked());
		if (p1.compare("Include Sediment traps")==0)         defnamelist[j].value.setNum((int)checkSedtrap->isChecked());
		if (p1.compare("Include wheeltracks")==0)            defnamelist[j].value.setNum((int)checkInfilCompact->isChecked());
		if (p1.compare("Include grass strips")==0)           defnamelist[j].value.setNum((int)checkInfilGrass->isChecked());
		if (p1.compare("Include crusts")==0)                 defnamelist[j].value.setNum((int)checkInfilCrust->isChecked());
		if (p1.compare("Impermeable sublayer")==0)           defnamelist[j].value.setNum((int)checkImpermeable->isChecked());
		//if (p1.compare("Matric head files")==0)              defnamelist[j].value.setNum((int)checkDumphead->isChecked());
		if (p1.compare("Geometric mean Ksat")==0)            defnamelist[j].value.setNum((int)checkGeometric->isChecked());
		if (p1.compare("Timeseries as PCRaster")==0)         defnamelist[j].value.setNum((int)checkWritePCRnames->isChecked());
		if (p1.compare("Timeplot as PCRaster")==0)           defnamelist[j].value.setNum((int)checkWritePCRtimeplot->isChecked());
		if (p1.compare("Regular runoff output")==0)          defnamelist[j].value.setNum((int)checkOutputTimeStep->isChecked());
		if (p1.compare("User defined output")==0)            defnamelist[j].value.setNum((int)checkOutputTimeUser->isChecked());
		if (p1.compare("No erosion at outlet")==0)           defnamelist[j].value.setNum((int)checkNoErosionOutlet->isChecked());
		if (p1.compare("Report point output separate")==0)   defnamelist[j].value.setNum((int)checkSeparateOutput->isChecked());
		if (p1.compare("Report point output for SOBEK")==0)  defnamelist[j].value.setNum((int)checkSOBEKOutput->isChecked());
		if (p1.compare("SOBEK date string")==0)              defnamelist[j].value = SOBEKdatestring->text();
		if (p1.compare("Sediment bulk density")==0)          defnamelist[j].value = E_BulkDens->text();
		if (p1.compare("Use canopy storage map")==0)   	     defnamelist[j].value.setNum((int)!checkInterceptionLAI->isChecked());
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
			defnamelist[j].value.setNum(i);
		}

		if (p1.compare("Begin time")==0) defnamelist[j].value = E_BeginTime->text();
		if (p1.compare("End time")==0)   defnamelist[j].value = E_EndTime->text();
		if (p1.compare("Timestep")==0)   defnamelist[j].value = E_Timestep->text();
		if (p1.compare("Map Directory")==0)    defnamelist[j].value = E_MapDir->text();
		if (p1.compare("Result Directory")==0) defnamelist[j].value = E_ResultDir->text();
		if (p1.compare("Main results file")==0) defnamelist[j].value = E_MainTotals->text();
		if (p1.compare("Filename point output")==0) defnamelist[j].value = E_PointResults->text();
		if (p1.compare("Rainfall Directory")==0) defnamelist[j].value = rainFileDir;
		if (p1.compare("Rainfall file")==0) defnamelist[j].value = E_RainfallName->text();
		if (p1.compare("Erosion map")==0) defnamelist[j].value = E_DetachmentMap->text();
		if (p1.compare("Deposition map")==0) defnamelist[j].value = E_DepositionMap->text();
		if (p1.compare("Soilloss map")==0) defnamelist[j].value = E_SoillossMap->text();
		if (p1.compare("Snowmelt Directory")==0) defnamelist[j].value = snowmeltFileDir;
		if (p1.compare("Snowmelt file")==0) defnamelist[j].value = E_SnowmeltName->text();
		if (p1.compare("Ksat calibration")==0) defnamelist[j].value = E_CalibrateKsat->text();
		if (p1.compare("N calibration")==0) defnamelist[j].value = E_CalibrateN->text();
		if (p1.compare("Channel Ksat calibration")==0) defnamelist[j].value = E_CalibrateChKsat->text();
		if (p1.compare("Channel N calibration")==0) defnamelist[j].value = E_CalibrateChN->text();
		if (p1.compare("Splash Delivery Ratio")==0) defnamelist[j].value = E_SplashDelibery->text();
		if (p1.compare("Stemflow fraction")==0) defnamelist[j].value = E_StemflowFraction->text();

		if (p1.compare("Table Directory")==0) defnamelist[j].value = SwatreTableDir;
		if (p1.compare("Table File")==0) defnamelist[j].value = SwatreTableName;

		if (p1.compare("Infil Method")==0)
		{
			switch(uiInfilMethod)
			{
			case 0 : defnamelist[j].value.setNum(INFIL_NONE); break;
			case 1 : defnamelist[j].value.setNum(INFIL_SWATRE);break;
			case 2 : if(checkInfil2layer->isChecked()) defnamelist[j].value.setNum(INFIL_GREENAMPT2);
				else defnamelist[j].value.setNum(INFIL_GREENAMPT);break;
			case 3 : if(checkInfil2layer->isChecked()) defnamelist[j].value.setNum(INFIL_SMITH2);
				else defnamelist[j].value.setNum(INFIL_SMITH); break;
			case 4: defnamelist[j].value.setNum(INFIL_KSAT); break;
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
			defnamelist[j].value = outputcheck.join(",");
		}
		//defnamelist[j].value = p;
	}
	// get all map names
	for (int j = 0; j < nrdefnamelist; j++)
		for (int i = 0; i < DEFmaps.size(); i++)
		{
			QStringList S = DEFmaps.at(i).split(";");
			if (S.contains(defnamelist[j].name))
				defnamelist[j].value = S.at(2);
		}

}
//---------------------------------------------------------------------------
// change runfile strings with current interface options
void lisemqt::DefaultRunFile()
{
	int i = 0, outputpos;
	defnamelist[i++].name = QString("[openLISEM runfile version 4]");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[LISEM main type]");
	defnamelist[i++].name = QString("LISEM Type");
	defnamelist[i++].name = QString("");
//	defnamelist[i++].name = QString("[Work Directory]");
//	defnamelist[i++].name = QString("WorkDir");
//	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Input]");
	defnamelist[i++].name = QString("Map Directory");
	defnamelist[i++].name = QString("Rainfall Directory");
	defnamelist[i++].name = QString("Rainfall file");
	defnamelist[i++].name = QString("Incude Snowmelt");
	defnamelist[i++].name = QString("Snowmelt Directory");
	defnamelist[i++].name = QString("Snowmelt file");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Output]");
	defnamelist[i++].name = QString("Result Directory");
	defnamelist[i++].name = QString("Main results file");
	defnamelist[i++].name = QString("Filename point output");
	defnamelist[i++].name = QString("Report point output separate");
	defnamelist[i++].name = QString("Report point output for SOBEK");
	defnamelist[i++].name = QString("SOBEK date string");
	defnamelist[i++].name = QString("Erosion map");
	defnamelist[i++].name = QString("Deposition map");
	defnamelist[i++].name = QString("Soilloss map");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Simulation times]");
	defnamelist[i++].name = QString("Begin time");
	defnamelist[i++].name = QString("End time");
	defnamelist[i++].name = QString("Timestep");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[General options]");
	defnamelist[i++].name = QString("No Erosion simulation");
	defnamelist[i++].name = QString("Include main channels");
	defnamelist[i++].name = QString("Include channel infil");
	defnamelist[i++].name = QString("Include channel baseflow");
	defnamelist[i++].name = QString("All water and sediment to outlet");
	defnamelist[i++].name = QString("Include snowmelt");
	defnamelist[i++].name = QString("No erosion at outlet");
	defnamelist[i++].name = QString("Alternative flow detachment");
	defnamelist[i++].name = QString("Simple depression storage");
	defnamelist[i++].name = QString("Hard Surfaces");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Interception]");
	defnamelist[i++].name = QString("Use canopy storage map");
	defnamelist[i++].name = QString("Canopy storage equation");
	defnamelist[i++].name = QString("Stemflow fraction");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Conservation]");
	defnamelist[i++].name = QString("Include grass strips");
	defnamelist[i++].name = QString("Grassstrip Mannings n");
	defnamelist[i++].name = QString("Include buffers");
	defnamelist[i++].name = QString("Sediment bulk density");
	defnamelist[i++].name = QString("Include Sediment traps");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Calibration]");
	defnamelist[i++].name = QString("Ksat calibration");
	defnamelist[i++].name = QString("N calibration");
	defnamelist[i++].name = QString("Channel Ksat calibration");
	defnamelist[i++].name = QString("Channel N calibration");
	defnamelist[i++].name = QString("Splash Delivery Ratio");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Gully options]");
	defnamelist[i++].name = QString("Fcrit relation");
	defnamelist[i++].name = QString("Threshold gradient");
	defnamelist[i++].name = QString("QW relation");
	defnamelist[i++].name = QString("QW param A");
	defnamelist[i++].name = QString("QW param B");
	defnamelist[i++].name = QString("Gully infiltration");
	defnamelist[i++].name = QString("Use initial gully dimensions");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Infiltration]");
	defnamelist[i++].name = QString("Infil Method");
	defnamelist[i++].name = QString("Include wheeltracks");
	defnamelist[i++].name = QString("Include crusts");
	defnamelist[i++].name = QString("Impermeable sublayer");
	defnamelist[i++].name = QString("Subsoil drainage");
	defnamelist[i++].name = QString("Table Directory");
	defnamelist[i++].name = QString("Table File");
	defnamelist[i++].name = QString("SWATRE internal minimum timestep");
	defnamelist[i++].name = QString("Matric head files");
	defnamelist[i++].name = QString("Geometric mean Ksat");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Output maps]");
	defnamelist[i++].name = QString("Runoff maps in l/s/m");
	defnamelist[i++].name = QString("Timeseries as PCRaster");
	defnamelist[i++].name = QString("Timeplot as PCRaster");
	defnamelist[i++].name = QString("Regular runoff output");
	defnamelist[i++].name = QString("Erosion map units (0/1/2)");
	defnamelist[i++].name = QString("Output interval");
	defnamelist[i++].name = QString("User defined output");
	defnamelist[i++].name = QString("Output times");
	defnamelist[i++].name = QString("CheckOutputMaps");
	defnamelist[i++].name = QString("CheckOutputMapsNUT");
	defnamelist[i++].name = QString("CheckOutputMapsMC");
	defnamelist[i++].name = QString("CheckOutputMapsGUL");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Texture classes]");
	defnamelist[i++].name = QString("ClassMu");
	defnamelist[i++].name = QString("[map names]");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputBASIC]");
	outputpos = i;
	defnamelist[i++].name = QString("OUTRUNOFF");
	defnamelist[i++].name = QString("OUTCONC");
	defnamelist[i++].name = QString("OUTWH");
	defnamelist[i++].name = QString("OUTRWH");
	defnamelist[i++].name = QString("OUTTC");
	defnamelist[i++].name = QString("OUTEROS");
	defnamelist[i++].name = QString("OUTDEPO");
	defnamelist[i++].name = QString("OUTVELO");
	defnamelist[i++].name = QString("OUTINF");
	defnamelist[i++].name = QString("OUTSS");
	defnamelist[i++].name = QString("OUTCHVOL");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputMC]");
	defnamelist[i++].name = QString("OUTMU0");
	defnamelist[i++].name = QString("OUTMU1");
	defnamelist[i++].name = QString("OUTMU2");
	defnamelist[i++].name = QString("OUTMU3");
	defnamelist[i++].name = QString("OUTMU4");
	defnamelist[i++].name = QString("OUTMU5");
	defnamelist[i++].name = QString("OUTD50SUSP");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputNut]");
	defnamelist[i++].name = QString("OUTPSOLUT");
	defnamelist[i++].name = QString("OUTPSUS");
	defnamelist[i++].name = QString("OUTPINF");
	defnamelist[i++].name = QString("OUTNH4SOLUT");
	defnamelist[i++].name = QString("OUTNH4SUS");
	defnamelist[i++].name = QString("OUTNH4INF");
	defnamelist[i++].name = QString("OUTNO3SOLUT");
	defnamelist[i++].name = QString("OUTNO3SUS");
	defnamelist[i++].name = QString("OUTNO3INF");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[OutputNutErosDep]");
	defnamelist[i++].name = QString("OUTPDEP");
	defnamelist[i++].name = QString("OUTNH4DEP");
	defnamelist[i++].name = QString("OUTNO3DEP");
	defnamelist[i++].name = QString("OUTPDET");
	defnamelist[i++].name = QString("OUTNH4DET");
	defnamelist[i++].name = QString("OUTNO3DET");
	defnamelist[i++].name = QString("");
   defnamelist[i++].name = QString("[OutputGul]");
	defnamelist[i++].name = QString("OUTGULD");
	defnamelist[i++].name = QString("OUTGULW");
	defnamelist[i++].name = QString("OUTGULA");
	defnamelist[i++].name = QString("OUTGULF");
	defnamelist[i++].name = QString("OUTGULDEM");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Catchment]");
	defnamelist[i++].name = QString("grad");
	defnamelist[i++].name = QString("ldd");
	defnamelist[i++].name = QString("outlet");
	defnamelist[i++].name = QString("ID");
	defnamelist[i++].name = QString("outpoint");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Landuse]");
	defnamelist[i++].name = QString("cover");
	defnamelist[i++].name = QString("lai");
	defnamelist[i++].name = QString("ch");
	defnamelist[i++].name = QString("smax");
	defnamelist[i++].name = QString("road");
	defnamelist[i++].name = QString("grasswidth");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Buffers]");
	defnamelist[i++].name = QString("bufferID");
	defnamelist[i++].name = QString("bufferVolume");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Snowmelt]");
	defnamelist[i++].name = QString("SnowID");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Erosion]");
	defnamelist[i++].name = QString("coh");
	defnamelist[i++].name = QString("cohadd");
	defnamelist[i++].name = QString("aggrstab");
	defnamelist[i++].name = QString("d50");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Surface]");
	defnamelist[i++].name = QString("rr");
	defnamelist[i++].name = QString("manning");
	defnamelist[i++].name = QString("crustfrc");
	defnamelist[i++].name = QString("compfrc");
	defnamelist[i++].name = QString("stonefrc");
	defnamelist[i++].name = QString("hardsurf");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilSwatre]");
	defnamelist[i++].name = QString("profmap");
	defnamelist[i++].name = QString("profcrst");
	defnamelist[i++].name = QString("profwltr");
	defnamelist[i++].name = QString("profgras");
	defnamelist[i++].name = QString("inithead");
	defnamelist[i++].name = QString("headout");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilExtra]");
	defnamelist[i++].name = QString("ksatcrst");
	defnamelist[i++].name = QString("ksatcomp");
	defnamelist[i++].name = QString("ksatgras");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilGA1]");
	defnamelist[i++].name = QString("ksat1");
	defnamelist[i++].name = QString("psi1");
	defnamelist[i++].name = QString("thetas1");
	defnamelist[i++].name = QString("thetai1");
	defnamelist[i++].name = QString("soildep1");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[InfilGA2]");
	defnamelist[i++].name = QString("ksat2");
	defnamelist[i++].name = QString("psi2");
	defnamelist[i++].name = QString("thetas2");
	defnamelist[i++].name = QString("thetai2");
	defnamelist[i++].name = QString("soildep2");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Channelinfil]");
	defnamelist[i++].name = QString("chanksat");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Channels]");
	defnamelist[i++].name = QString("lddchan");
	defnamelist[i++].name = QString("chanwidth");
	defnamelist[i++].name = QString("chanside");
	defnamelist[i++].name = QString("changrad");
	defnamelist[i++].name = QString("chanman");
	defnamelist[i++].name = QString("chancoh");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[ChannelBaseflow]");
	defnamelist[i++].name = QString("chanbaseflux");
	defnamelist[i++].name = QString("chanincrease");
	defnamelist[i++].name = QString("chanvolini");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Wheeltrack]");
	defnamelist[i++].name = QString("lddwheel");
	defnamelist[i++].name = QString("wheelnbr");
	defnamelist[i++].name = QString("wheelwidth");
	defnamelist[i++].name = QString("wheeldepth");
	defnamelist[i++].name = QString("wheelgradient");
	defnamelist[i++].name = QString("wheelman");
	defnamelist[i++].name = QString("wheelcohesion");
	defnamelist[i++].name = QString("ksatwt");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Texture]");
	defnamelist[i++].name = QString("fractionmu0");
	defnamelist[i++].name = QString("fractionmu1");
	defnamelist[i++].name = QString("fractionmu2");
	defnamelist[i++].name = QString("fractionmu3");
	defnamelist[i++].name = QString("fractionmu4");
	defnamelist[i++].name = QString("fractionmu5");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsP]");
	defnamelist[i++].name = QString("pcont");
	defnamelist[i++].name = QString("psolute");
	defnamelist[i++].name = QString("pefficiency");
	defnamelist[i++].name = QString("psorp");
	defnamelist[i++].name = QString("pconv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsNO3]");
	defnamelist[i++].name = QString("no3cont");
	defnamelist[i++].name = QString("no3solute");
	defnamelist[i++].name = QString("no3efficiency");
	defnamelist[i++].name = QString("no3sorp");
	defnamelist[i++].name = QString("no3conv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsNH4]");
	defnamelist[i++].name = QString("nh4cont");
	defnamelist[i++].name = QString("nh4solute");
	defnamelist[i++].name = QString("nh4efficiency");
	defnamelist[i++].name = QString("nh4sorp");
	defnamelist[i++].name = QString("nh4conv");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[NutsBD]");
	defnamelist[i++].name = QString("bulk");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[Gully]");
	defnamelist[i++].name = QString("dem");
	defnamelist[i++].name = QString("gullyn");
	defnamelist[i++].name = QString("bulkdens1");
	defnamelist[i++].name = QString("gulksat1");
	defnamelist[i++].name = QString("gullydep");
	defnamelist[i++].name = QString("gullycoh");
	defnamelist[i++].name = QString("bulkdens2");
	defnamelist[i++].name = QString("gulksat2");
	defnamelist[i++].name = QString("nonfcrit");
	defnamelist[i++].name = QString("");
	defnamelist[i++].name = QString("[GullyInit]");
	defnamelist[i++].name = QString("gulwinit");
	defnamelist[i++].name = QString("guldinit");
	nrdefnamelist = i;
	for (i = 0; i < nrdefnamelist; i++)
		defnamelist[i].value.clear();
   i = outputpos;
	defnamelist[i++].value=QString("ro");
	defnamelist[i++].value=QString("conc");
	defnamelist[i++].value=QString("wh");
	defnamelist[i++].value=QString("whc");
	defnamelist[i++].value=QString("tc");
	defnamelist[i++].value=QString("det");
	defnamelist[i++].value=QString("depo");
	defnamelist[i++].value=QString("velo");
	defnamelist[i++].value=QString("inf");
	defnamelist[i++].value=QString("sstor");
	defnamelist[i++].value=QString("chanvol");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("smu0");
	defnamelist[i++].value=QString("smu1");
	defnamelist[i++].value=QString("smu2");
	defnamelist[i++].value=QString("smu3");
	defnamelist[i++].value=QString("smu4");
	defnamelist[i++].value=QString("smu5");
	defnamelist[i++].value=QString("D50s");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("NPsol");
	defnamelist[i++].value=QString("NPsus");
	defnamelist[i++].value=QString("NPinf");
	defnamelist[i++].value=QString("NNH4sol");
	defnamelist[i++].value=QString("NNH4sus");
	defnamelist[i++].value=QString("NNH4inf");
	defnamelist[i++].value=QString("NNO3sol");
	defnamelist[i++].value=QString("NNO3sus");
	defnamelist[i++].value=QString("NNO3inf");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("NPdep.map");
	defnamelist[i++].value=QString("NNH4dep.map");
	defnamelist[i++].value=QString("NNO3dep.map");
	defnamelist[i++].value=QString("NPdet.map");
	defnamelist[i++].value=QString("NNH4det.map");
	defnamelist[i++].value=QString("NNO3det.map");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("");
	defnamelist[i++].value=QString("guld");
	defnamelist[i++].value=QString("gulw");
	defnamelist[i++].value=QString("gula");
	defnamelist[i++].value=QString("gulf");
	defnamelist[i++].value=QString("");
}
//---------------------------------------------------------------------------
