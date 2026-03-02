/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file swatinp.cpp
  \brief SWATRE: initialize and read profile data

  functions:
- void TWorld::ReadSwatreInputNew(void) \n
- ZONE * TWorld::ReadNodeDefinitionNew(void) \n
- PROFILE * TWorld::ReadProfileDefinitionNew(int pos,ZONE *z) \n
- HORIZON * TWorld::ReadHorizonNew(const char *tablePath, const char *tableName) \n
- void  TWorld::FreeSwatreInfo(void) \n

profile node setup:
    endComp is what is in the profile.inp file, the bottom of the layer
    dz = (endComp[i-1] - endComp[i]) is negative layer thickness
    z = 0.5*(dz[i-1]+dz[i]) is negative centre of compartment, nodes
    disnod = z[i]-z[i-1] is negative distance between centres, nodes
    NOTE: Independent of individuual soil profiles table depths!!!

     -------   surface    -       - z[0]-
        o                  |dz[0] -      | disnod[0]
     -------   endComp[0] -        |z[1]-
        o                  |dz[1] -      | disnod[1]
     -------   endcomp[1] -        |z[2]-
        o                  |dz[2] -      | disnod[2]
     -------   endcomp[2] -
    etc.
*/

#include "lerror.h"
#include "model.h"

#define LIST_INC	20

#define ROOTMAX 60  // rootzone depth

//----------------------------------------------------------------------------------------------
/// read and parse profile.inp
/// new version using swatreProfileDef QStringList
void TWorld::ReadSwatreInputNew(void)
{
//qDebug() << "ReadSwatreInputNew";

    // reset stuff
    int nrProfileList = 0;
    sizeProfileList = 0;
    nrHorizonList = 0;
    sizeHorizonList = 0;
    swatreProfileDef.clear();
    swatreProfileNr.clear();

    QFile file(SwatreTableName); // table name has full path (profile.inp)

    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);

        while (!in.atEnd()) {
            QString line = in.readLine();

            // Skip empty or space-only lines
            if (!line.trimmed().isEmpty()) {
                swatreProfileDef.append(line);
            } else {
                swatreProfileDef.append("###"); // replace empty lines "###" for parsing
            }
        }

        file.close();
    } else {
        Error(QString("SWATRE: Can't open profile definition file %1").arg(SwatreTableName));
        throw 1;
    }

    for (int i = swatreProfileDef.size() - 1; i >= 0; --i) {
        if (swatreProfileDef[i].trimmed().isEmpty()) {
            swatreProfileDef.removeAt(i);
        }
    }
   while(swatreProfileDef.last() == "###")
        swatreProfileDef.removeLast();

   // for(int i; i < swatreProfileDef.count(); i++)
   //     qDebug() << swatreProfileDef[i];

    // read and make nodes
    zone = new ZONE;

    bool ok;
    zone->nrNodes = swatreProfileDef[0].toInt(&ok, 10);
    if (!ok)
        Error(QString("SWATRE: Can't read number of nodes %1 from input file: %2").arg(zone->nrNodes).arg(SwatreTableName));
    if (zone->nrNodes < 3 )
        Error(QString("SWATRE: you need to define 3 nodes or more").arg(zone->nrNodes));
    if (zone->nrNodes > MAX_NODES)
        Error(QString("SWATRE: number of nodes %1 larger than %2").arg(zone->nrNodes).arg(MAX_NODES));

    for (int i=0; i <= zone->nrNodes; i++) {
        zone->dz.append(0.0);
        zone->z.append(0.0);
        zone->endComp.append(0.0);
        zone->disnod.append(0.0);
        zone->rootz.append(0.0);
    }

    int pos = 2;
    for (int i = 0; i < zone->nrNodes; i++) {
        bool ok;
        zone->endComp[i] = toDOUBLE(swatreProfileDef[i+pos],&ok);//swatreProfileDef[i+pos].toDouble(&ok);
        if (!ok)
            Error(QString("SWATRE: Can't read compartment end of node %1").arg(i+pos));
        if (zone->endComp[i] <= 0)
            Error(QString("SWATRE: compartment end of node %1 <= 0").arg(i+pos));
    }
    zone->dz[0]= -zone->endComp[0];
    zone->z[0]= zone->dz[0]*0.5;
    zone->disnod[0] = zone->z[0];
    zone->rootz[0] = 0;
    double rootmax = -ROOTMAX; // 60 cm fixed independent of crop
    double sum = 0;
    for (int i = 1; i < zone->nrNodes; i++) {
        zone->dz[i]= (zone->endComp[i-1]-zone->endComp[i]);
        zone->z[i]= zone->z[i-1] + 0.5*(zone->dz[i-1]+zone->dz[i]);
        zone->disnod[i] = zone->z[i] - zone->z[i-1];
        if (zone->z[i] > rootmax) {
            zone->rootz[i] = (rootmax - zone->z[i])/rootmax;
            sum = sum + zone->rootz[i];
        }
    }
    // tapering distribution of from layer 1 to < rootmax, sum is 1
    for (int i = 1; i < zone->nrNodes; i++) {
        if (zone->z[i] > rootmax) {
            zone->rootz[i] /= sum;
        }
    }
    // for (int i = 0; i < zone->nrNodes; i++) {
    //      qDebug() << zone->rootz[i];
    // }
    zone->disnod[zone->nrNodes] = 0.5 * zone->dz[zone->nrNodes-1];

 // for (int i = 0; i <= zone->nrNodes; i++)
   //    qDebug() << i << "dz" << zone->dz[i] << "z" << zone->z[i] << "dist" << zone->disnod[i] << "root" << zone->rootz[i];

    //  count and check valid profiles
    QStringList checkList; // temp list to check for double profile nrs
    for (int i = zone->nrNodes+1; i < swatreProfileDef.count(); i++) {
        if (swatreProfileDef[i].contains("###")) {
            checkList << swatreProfileDef[i+1];
            nrProfileList++;
        }
    }
    sizeProfileList = nrProfileList;

    if (nrProfileList == 0)
        Error(QString("SWATRE: no profiles read from %1").arg(SwatreTableName));

    // sort profile on increasing number
    swatreProfileNr.clear();
    for (int i = 0; i < checkList.count(); i++)
        swatreProfileNr << checkList[i].toInt();
   // std::sort(swatreProfileNr.begin(), swatreProfileNr.end());
    // DO NOT SORT! IT BREAKS THE LINK BETWEEN THIS LIST AND THE MAP

    for (int i = 0; i < swatreProfileNr.count()-1; i++)
    {
       // qDebug() << swatreProfileNr[i];
        if (swatreProfileNr[i] == swatreProfileNr[i+1])
            DEBUG(QString("Warning SWATRE: profile id %1 defined more than once").arg(swatreProfileNr[i+1]));
    }

    // profile list is a list of pointers to PROFILE
    profileList = new PROFILE*[nrProfileList + 1];

    nrProfileList = 0;
   // qDebug() << swatreProfileDef.count();
    for (int i = zone->nrNodes+1; i < swatreProfileDef.count(); i++) {
        if (swatreProfileDef[i].contains("###")) {
            // start of a new profilel
            i++;
        //    qDebug() << swatreProfileDef[i];
            profileList[nrProfileList] = ReadProfileDefinitionNew(i, zone);
            // creates a profile and gives the pointer to profilelist
            // i is the place in the StrinList where a profile starts
            nrProfileList++;
         //   qDebug() << "nrProfileList" << nrProfileList;
        }
    }

    for (int i = 0; i < ProfileIDList.count(); i++) {
        if (ProfileIDList[i]> 0 && !swatreProfileNr.contains(ProfileIDList[i])) {
            Error(QString("SWATRE: map profile nr %1 does not exist in profile table: %2").arg(ProfileIDList[i]).arg(SwatreTableName));
        }
    }



}
//----------------------------------------------------------------------------------------------
// for reference:
// typedef struct PROFILE {
//     int            profileId; 	/** number identifying this profile  >= 0 */
//     const ZONE     *zone; 		/** array with zone.nrNodes elements: containing z, dz, node distance etc*/
//     const HORIZON  **horizon; 	/** ptr to horizon information this node belongs to */
//     QVector <double> KsatCal;
// } PROFILE;

//note: all horizons and tables are read that are needed, not only those in profile.map
// but the horizons are in each profile with a pointer, not a full copy
PROFILE * TWorld::ReadProfileDefinitionNew(int pos, ZONE *z)
{
    QString tableName;
    double endHor = 0, endHorPrev = 0;
    PROFILE *p;
    HORIZON *hor;
    bool ok;

    p = new PROFILE;

    p->profileId = swatreProfileDef[pos].toInt(&ok, 10);

    if (!ok)
        Error(QString("SWATRE: read error: error in profile id %1 definition").arg(p->profileId));

    p->horizon = (const HORIZON **)malloc(sizeof(HORIZON *) * z->nrNodes); // array of pointers to horizon
    p->zone = z; // also pointer to zone info
    for (int i = 0; i < z->nrNodes; i++)
        p->KsatCal << 1.0; // create ksat cal 1,2,3 for each horizon, initialize to 1.0

    int i = 0;
    int hornr = 0;
    while (i != z->nrNodes) {
        pos++; // move one line to read the horizon table name

        tableName = swatreProfileDef[pos];
        if (!QFileInfo(SwatreTableDir + tableName).exists())
            Error(QString("SWATRE: Table %1 not found.").arg(tableName));
        //"Can't read the LUT for profile nr %1 node nr %2 and up").arg(p->profileId).arg(i+1));

        endHorPrev = endHor;
        pos++; // move one line to read the horizon depth endhor in cm
        hornr++; // increase horizon number

        endHor = swatreProfileDef[pos].toDouble(&ok);
        if (!ok)
            Error(QString("SWATRE: Can't read end of horizon for profile nr %1").arg(p->profileId));
        if (endHor <= endHorPrev)
            Error(QString("SWATRE: Error in profile definition nr %1, depth horizons does not increase").arg(p->profileId));

        // read the horizon and the luts for each node
        hor = ReadHorizonNew(SwatreTableDir, tableName);

        // copy horizon info to all nodes of this horizon
        // add the proper calibration factor (ksat1 cal for hor 1, ksat2cal for hor 2 adn the rest hor 3)
        while (i < z->nrNodes && z->endComp[i] <= endHor) {

            p->horizon[i] = hor;

            if (hornr == 1) p->KsatCal.replace(i, ksatCalibration);
            if (hornr == 2)  p->KsatCal.replace(i, ksat2Calibration);
            if (hornr > 2)  p->KsatCal.replace(i, ksat3Calibration);

         //   qDebug() << i << p->profileId << p->horizon[i]->name << z->endComp[i] << endHor;
                        //hornr << p->profileId << p->horizon[i]->name << z->endComp[i] << z->z[i] << z->disnod[i] << z->dz[i];
            i++;
        }

        // if (z->endComp[i-1] != endHor)
        //     Error(QString("SWATRE: Compartment does not end on depth '%1' (found in profile nr %2 for horizon %3)")
        //           .arg(endHor).arg(p->profileId).arg(tableName));
        //? what does this error mean exactly, horizons do not have to end exacvtly on nodes

    }

    return(p);
}


//----------------------------------------------------------------------------------------------
// copy horizon info to all nodes of this horizon
HORIZON * TWorld::ReadHorizonNew(QString tablePath, QString tableName)
{
    //qDebug() << "ReadHorizonNew" << nrHorizonList;
    // look if it's already loaded
    for(int i = 0; i < nrHorizonList; i++)
        if (tableName == horizonList[i]->name)
            return(horizonList[i]);

    // check for space in list, if not then add LIST_INC (20) pointers
    if (nrHorizonList == sizeHorizonList) {
        sizeHorizonList += LIST_INC; // add 20
        horizonList = (HORIZON **)realloc(horizonList, sizeof(HORIZON *)*sizeHorizonList);
    }

    // make a horizon and add it to the list, horizon is a name and a pointer to LUT
    HORIZON	*h = new HORIZON;
    horizonList[nrHorizonList++] = h;

    // read the lut with this horizon and link the pointer
    h->lut = ReadSoilTableNew(tablePath + tableName);
    h->name = tableName;

    return(h);
}
//----------------------------------------------------------------------------------------------
LUT *TWorld::ReadSoilTableNew(QString fileName)
{
    // read the table in a stringlist
    QStringList list;

    QRegularExpression regex("[a-df-zA-DF-Z]");

    QFile file(fileName);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);

        int j=1;
        while (!in.atEnd()) {
            QString line = in.readLine();
            if (regex.match(line).hasMatch()) {
                ErrorString = QString("SWATRE: Please check: Invalid characters found in file %1, line %2: [%3]").arg(fileName).arg(j).arg(line);
                throw 1;
            }
            j++;
            line.replace("\u001A"," "); // some limburg tables have the char "substitute" in them, this is just a hack
            // Skip empty or space-only lines
            if (!line.trimmed().isEmpty() && line != " ")
                list.append(line);

        }
        file.close();
    }
    LUT *l = new LUT;
    l->Rows = list.count();

    for (int i = 0; i < list.count(); i++) {
        QStringList SL = list[i].split(QRegularExpression("\\s+"),Qt::SkipEmptyParts);

        bool ok;
        double v = toDOUBLE(SL[0],&ok);
        if (!ok || SL.count() < 3) {
            l->Rows--;
            break; // sometimes table ends with a non empty line with some char code
        }
        l->hydro[THETA_COL].append(toDOUBLE(SL[THETA_COL], &ok));
        if (!ok)
            Error(QString("Cannot read theta in table %1.").arg(fileName));
        l->hydro[H_COL].append(toDOUBLE(SL[H_COL], &ok));
        if (!ok)
            Error(QString("Cannot read h in table %1.").arg(fileName));
        l->hydro[K_COL].append(toDOUBLE(SL[K_COL], &ok)/86400.0); // cm/day to cm/sec
        if (!ok)
            Error(QString("Cannot read K in table %1.").arg(fileName));
    }

    for (int i = 0; i < l->Rows-1; i++) {
        if (l->hydro[H_COL][i+1] <= l->hydro[H_COL][i])
            Error(QString("matrix head not increasing in table %1 at h = %2.").arg(fileName).arg(l->hydro[H_COL][i]));
        if (l->hydro[THETA_COL][i+1] <= l->hydro[THETA_COL][i])
            Error(QString("moisture content not increasing in table %1 at theta = %2.").arg(fileName).arg(l->hydro[THETA_COL][i]));
        if ((l->hydro[K_COL][i+1] < l->hydro[K_COL][i]) && (l->hydro[K_COL][i+1] > 0 && l->hydro[K_COL][i] > 0))
            Error(QString("Hydraulic conductivity not increasing in table %1 at K = %2.").arg(fileName).arg(l->hydro[K_COL][i]));
    }

    // dimoica
    for (int i = 0; i < l->Rows - 1; i++) {
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
    l->hydro[DMCC_COL] << l->hydro[DMCC_COL][l->Rows-2] + (l->hydro[DMCC_COL][l->Rows-2] - l->hydro[DMCC_COL][l->Rows-3]);

    return(l);
}
//----------------------------------------------------------------------------------------------
void TWorld::checkFileForInvalidLetters(const QString &filePath)
{
    QFile file(filePath);
    QTextStream in(&file);
    // Regex to match any letters except 'E' and 'e' used in scientific notation
    QRegularExpression regex("[a-df-zA-DF-Z]");
    int lineNumber = 1;

    while (!in.atEnd()) {
        QString line = in.readLine();
        if (regex.match(line).hasMatch()) {
            ErrorString = QString("SWATRE: Please check: Invalid characters found in file %1, line %2: %3").arg(filePath).arg(lineNumber).arg(line);
            throw 1;
        }
        lineNumber++;
    }
    file.close();
}
//----------------------------------------------------------------------------------------------
