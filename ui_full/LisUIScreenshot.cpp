/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2022  Victor Jetten
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

#include "lisemqt.h"
#include "global.h"

//--------------------------------------------------------------------
void lisemqt::shootMScreen()
{
    doShootScreens = shootMscreenAct->isChecked();
}

//--------------------------------------------------------------------
void lisemqt::shootScreen()
{
    shootSingleScreen(0);
/*
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this, "openLISEM",QString("Select a run file first"));
        return;
    }

    QPixmap originalPixmap;
    QString format = "png";
    QFileInfo fi(op.runfilename);

    QString fileName = screenShotDir + fi.baseName();
    QString number = QString("-%1").arg(op.runstep,5,'d',0,'0');
    QString name;

    // is multiple screenshots while running
    if (doShootScreens)
    {
        if (op.runstep % printinterval->value() > 0)
            return;
        // continue each printinterval else return

        // always take screenshot of the hydrograph
        tabWidget_out->setCurrentIndex(0);
        originalPixmap = tabWidget->widget(2)->grab();
        fileName = screenShotDir + fi.baseName()+ "_Q" + number + ".png";
        originalPixmap.save(fileName, format.toLatin1());

        // grab the map screen, whatever map is present
        tabWidget_out->setCurrentIndex(1);
        originalPixmap = tabWidget->widget(2)->grab();
        // get the proper name of the variable shown
        QString name = "";
        if (checkBoxComboMaps->isChecked()) {
            int index = DisplayComboBox->currentIndex();
            if( index > -1 && index < NameList.length())
                name = NameList.at(index);
        } else if (checkBoxComboMaps2->isChecked()) {
            int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
         //   qDebug() << index;
            if( index > -1 && index < NameList.length())
                name = NameList.at(index);
        }

        fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
        originalPixmap.save(fileName, format.toLatin1());
    }
    else // only one screenshot
    {
        if (tabWidget->currentIndex() == 0) {
            originalPixmap = tabWidgetOptions->currentWidget()->grab();
            name = QString("input_%1").arg(tabWidgetOptions->currentIndex());
            number = "";
        } else
            originalPixmap = tabWidget->currentWidget()->grab();

        if (tabWidget->currentIndex() == 1)
            name = "inputmaps";
        if (tabWidget->currentIndex() == 2) // output
        {
            if (tabWidget_out->currentIndex() == 0) {
                name = "_Q";
            }
            if (tabWidget_out->currentIndex() == 1) {
                if (checkBoxComboMaps->isChecked()) {
                    int index = DisplayComboBox->currentIndex();
                    if( index > -1 && index < NameList.length())
                        name = NameList.at(index);
                } else
                    if (checkBoxComboMaps2->isChecked()) {
                    int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
                    if( index > -1 && index < NameList.length())
                        name = NameList.at(index);
                }
            }
        }

        fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
        originalPixmap.save(fileName, format.toLatin1());
    }
*/
}
//--------------------------------------------------------------------
void lisemqt::shootMultipleScreens()
{
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this, "openLISEM",QString("Select a run file first"));
        return;
    }

    QPixmap originalPixmap;
    QString format = "png";
    QFileInfo fi(op.runfilename);

    QString fileName = screenShotDir + fi.baseName();
    QString number = QString("-%1").arg(op.runstep,5,'d',0,'0');


    if (op.runstep % printinterval->value() > 0)
        return;
    // continue each printinterval else return

    // always take screenshot of the hydrograph
    tabWidget_out->setCurrentIndex(0);
    originalPixmap = tabWidget->widget(2)->grab();
    fileName = screenShotDir + fi.baseName()+ "_Q" + number + ".png";
    originalPixmap.save(fileName, format.toLatin1());

    // grab the map screen, whatever map is present
    tabWidget_out->setCurrentIndex(1);
    originalPixmap = tabWidget->widget(2)->grab();
    // get the proper name of the variable shown
    QString name = "";
    if (checkBoxComboMaps->isChecked()) {
        int index = DisplayComboBox->currentIndex();
        if( index > -1 && index < NameList.length())
            name = NameList.at(index);
    } else if (checkBoxComboMaps2->isChecked()) {
        int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
        //   qDebug() << index;
        if( index > -1 && index < NameList.length())
            name = NameList.at(index);
    }

    fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
    originalPixmap.save(fileName, format.toLatin1());

}
//--------------------------------------------------------------------
void lisemqt::shootSingleScreen(int options)
{
    if (op.runfilename.isEmpty())
    {
        QMessageBox::warning(this, "openLISEM",QString("Select a run file first"));
        return;
    }

    QPixmap originalPixmap;
    QString format = "png";
    QFileInfo fi(op.runfilename);

    QString fileName = screenShotDir + fi.baseName();
    QString number = "";//QString("-%1").arg(op.runstep,5,'d',0,'0');
    QString name;

    if (options == 1) {
        tabWidget->setCurrentIndex(0);
        originalPixmap = tabWidgetOptions->currentWidget()->grab();
        name = QString("input_%1").arg(tabWidgetOptions->currentIndex());
        number = "";
    } else {
        originalPixmap = tabWidget->currentWidget()->grab();
        // grab the screen

        if (tabWidget->currentIndex() == 1) { // maplist
            name = "inputmaps";
        }
        if (tabWidget->currentIndex() == 2) // output
        {
            number = QString("-%1").arg(op.runstep,5,'d',0,'0');
            if (tabWidget_out->currentIndex() == 0) {
                name = "_Q";
            }
            if (tabWidget_out->currentIndex() == 1) {
                if (checkBoxComboMaps->isChecked()) {
                    int index = DisplayComboBox->currentIndex();
                    if( index > -1 && index < NameList.length())
                        name = NameList.at(index);
                } else
                    if (checkBoxComboMaps2->isChecked()) {
                        int index = DisplayComboBox2->currentIndex()+DisplayComboBox->count();
                        if( index > -1 && index < NameList.length())
                            name = NameList.at(index);
                    }
            }

        }
    }
    fileName = screenShotDir + fi.baseName()+ name + number  + ".png";
    originalPixmap.save(fileName, format.toLatin1());
}
//--------------------------------------------------------------------
void lisemqt::convertScreenshotsToVideo()
{
    lisMpeg->setMencoderDir(mencoderDir);
    lisMpeg->setWorkDir(E_ResultDir->text());

    lisMpeg->exec();

    mencoderDir = lisMpeg->mencoderDir;
    qDebug() << "new" << mencoderDir;
}
