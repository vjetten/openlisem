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

/*
 * UI basic interface: very simple ntest interface, depreciated
 */

#ifndef ifacebasic_H
#define ifacebasic_H

#include <QtGui>

#include "ui_ifacebasic.h"
#include "model.h"
#include "lisuioutput.h"

/// Very simple basic test interface, depreciated
class ifacebasic : public QWidget, public Ui::ifacebasicClass
{
    Q_OBJECT

public:
    ifacebasic(QWidget *parent = 0);
    ~ifacebasic();

private slots:
    void on_runButton_clicked();
    void on_toolButton_ShowRunfile_clicked();
    void Showit();
    void worldDone(const QString &results);
    void worldDebug(const QString &results);
   // void on_checkChannel_clicked();
   // void on_checkErosion_clicked();
    void on_toolButton_runfilename_clicked();
    void StorePath();
    void GetStorePath();
    void SetStyleUI();

private:
    TWorld *W;
};


#endif // ifacebasic_H
