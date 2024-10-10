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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
\file lutio.cpp
\brief  SWATRE read theta, h, k-table

functions:
- double* TWorld::ReadSoilTable(const char *fileName, int *nrRows); \n
- void TWorld::ReadCols(const char *fileName, double *inLut, const char *buf); \n
*/


//#include "swatresoillut.h"
#include "lerror.h"
#include "model.h"


// ///  number of elements added to the lut when initializing (malloc) or resizing (realloc)
// #define EXP_NR_COLS 3
// #define INC_STEP (45*3)

// /// error messages:
// #define OPEN_ERRORs "SWATRE: Can't open %1"
// #define READ_ERRORs "SWATRE: Read error on %1"
// #define COL_ERRORs  "SWATRE: Table %1, entry nr. %2 contains %3 than 3 columns"
// #define EOF_MESSs   "SWATRE: Unexpected end of file while reading lookup table"
// #define NO_ENTRIESs "SWATRE: Table %1 contains no entries"
// #define NR_COLSs    "SWATRE: Encountered line containing %1 while first row had %2 columns"
// #define NAN_MESSs   "SWATRE: Table %1 contains a non number symbol: %2"
// #define SMALLERs    "SWATRE: Table %1 column %2 on entry %3 has smaller value (%4) than previous element"

static const char *colName[3] = { "theta", "h", "k" };

//----------------------------------------------------------------------------------------
//	const char *fileName,
LUT *TWorld::ReadSoilTableNew(QString fileName)
{
    QChar subChar(26); //SUB
    // read the table in a stringlist
    QStringList list;
    QFile file(fileName);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);

        while (!in.atEnd()) {
            QString line = in.readLine();
           // line.remove(subChar);
            line.replace("\u001A"," ");
            // Skip empty or space-only lines
            if (!line.trimmed().isEmpty() && line != " ")
                list.append(line);           

        }
        file.close();
    }
    LUT *l = new LUT;
    l->nrRows = list.count();

    for (int i = 0; i < list.count(); i++) {
        QStringList SL = list[i].split(QRegularExpression("\\s+"),Qt::SkipEmptyParts);

        bool ok;
        SL[0].toDouble(&ok);
        if (!ok || SL.count() < 3) {
            qDebug() << "not ok" << SL;
            l->nrRows--;
            break; // sometimes table ends with some char code
        }
        l->hydro[THETA_COL].append(SL[THETA_COL].toDouble());
        l->hydro[H_COL].append(SL[H_COL].toDouble());
        l->hydro[K_COL].append(SL[K_COL].toDouble()/86400 * ksatCalibration); // cm/day to mm/h 0.41667!
        //NOTE: Ksat in cm/day needs to me cm/sec ??
    }

    for (int i = 0; i < l->nrRows-1; i++) {
        if (l->hydro[H_COL][i+1] <= l->hydro[H_COL][i])
            Error(QString("matrix head not increasing in table %1 at h = %2.").arg(fileName).arg(l->hydro[H_COL][i]));
        if (l->hydro[THETA_COL][i+1] <= l->hydro[THETA_COL][i])
            Error(QString("moisture content not increasing in table %1 at theta = %2.").arg(fileName).arg(l->hydro[THETA_COL][i]));
        if (l->hydro[K_COL][i+1] <= l->hydro[K_COL][i])
            Error(QString("Hydraulic conductivity not increasing in table %1 at K = %2.").arg(fileName).arg(l->hydro[K_COL][i]));
    }

    for (int i = 0; i < l->nrRows - 1; i++) {
        double v = 0.5*(l->hydro[H_COL][i] + l->hydro[H_COL][i+1]);
        l->hydro[DMCH_COL] << v;
        v = (l->hydro[THETA_COL][i+1] - l->hydro[THETA_COL][i])/(l->hydro[H_COL][i+1] - l->hydro[H_COL][i]);

        if (i > 0 && v < l->hydro[DMCC_COL][i-1]) {
            double dv = l->hydro[DMCC_COL][i-1] - l->hydro[DMCC_COL][i-2];
            v = l->hydro[DMCC_COL][i-1]+dv;
        }
        l->hydro[DMCC_COL] << v;
    }

    // fill l->nrRows-1
    l->hydro[DMCH_COL] << 0;
    l->hydro[DMCC_COL] << l->hydro[DMCC_COL][l->nrRows-2] + (l->hydro[DMCC_COL][l->nrRows-2] - l->hydro[DMCC_COL][l->nrRows-3]);

    // qDebug() << fileName;
    // for (int i = 0; i < l->nrRows; i++) {
    // qDebug() << l->hydro[0][i] <<  l->hydro[1][i] << l->hydro[2][i] << l->hydro[3][i] << l->hydro[4][i];
    // }

    // WORKS

    return(l);
}
