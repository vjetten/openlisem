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
  \file swatinp.cpp
  \brief SWATRE: initialize and read profile data

  functions:
- void TWorld::InitializeProfile( void )
- void TWorld::ReadSwatreInputNew(void) \n
- ZONE * TWorld::ReadNodeDefinitionNew(void) \n
- PROFILE * TWorld::ReadProfileDefinitionNew(int pos,ZONE *z) \n
- HORIZON * TWorld::ReadHorizon(const char *tablePath, const char *tableName) \n
- void  TWorld::FreeSwatreInfo(void) \n

- obsolete !!!!!!!!:
- int TWorld::ReadSwatreInput(QString fileName, QString tablePath) \n
- ZONE * TWorld::ReadNodeDefinition(FILE *f) \n
- PROFILE * TWorld::ReadProfileDefinitionNew(FILE *f,ZONE *z,const char *tablePath) \n
- PROFILE * TWorld::ProfileNr(int profileNr) \n

profile node setup:
    endComp is what is in the profile.inp file, the bottom of the layer
    dz = (endComp[i-1] - endComp[i]) is negative layer thickness
    z = 0.5*(dz[i-1]+dz[i]) is negative centre of compartment, nodes
    disnod = z[i]-z[i-1] is negative distance between centres, nodes

     -------   surface    -       - z[0]-
        o                  |dz[0] -      | disnod[0]
     -------   endComp[0] -        |z[1]-
        o                  |dz[1] -      | disnod[1]
     -------   endcomp[1] -        |z[2]-
        o                  |dz[2] -      | disnod[2]
     -------   endcomp[2] -
    etc.
*/

#include <algorithm>
#include "lerror.h"
#include "model.h"

#define LIST_INC	20
#define LUT_COLS  5
#define IND(r,c)  ((r)*LUT_COLS+(c))

/// array of pointers to horizons, nullptr if not allocated
static HORIZON **horizonList = nullptr;
static int nrHorizonList=0, sizeHorizonList=0;

//----------------------------------------------------------------------------------------------
void TWorld::InitializeProfile( void )
{

    zone = nullptr;
    profileList = nullptr;
    //FreeSwatreInfo();

    //profileList = nullptr;

    // if (profileList != nullptr) {
    //     if (profileList[0] != nullptr) {
    //          for(int i=0; i < nrProfileList; i++)
    //             if (profileList[i] != nullptr)
    //                 free(profileList[i]);
    //     }
    //     free(profileList);
    //     profileList = nullptr;
    // }

    nrProfileList = 0;
    sizeProfileList = 0;

    // zone=nullptr;
    horizonList = nullptr;
    nrHorizonList = 0;
    sizeHorizonList = 0;

    //swatreProfileDef.clear();
    //swatreProfileNr.clear();
}
//----------------------------------------------------------------------------------------------
/// read and parse profile.inp
/// new version using swatreProfileDef QStringList
void TWorld::ReadSwatreInputNew(void)
{

    // get the profile inp file contents and put in stringlist
    // set profiles to nullptr
    InitializeProfile();

    QFile file(SwatreTableName);

    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);

        while (!in.atEnd()) {
            QString line = in.readLine();

            // Skip empty or space-only lines
            if (!line.trimmed().isEmpty()) {
                swatreProfileDef.append(line);
            } else {
                swatreProfileDef.append("###");
            }
        }

        file.close();
    } else {
        Error(QString("SWATRE: Can't open profile definition file %1").arg(SwatreTableName));
        throw 1;
    }
//qDebug() << swatreProfileDef;
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
    }

    int pos = 1;
    for (int i = 1; i < zone->nrNodes; i++) {
        pos++;

        bool ok;
        zone->endComp[i] = swatreProfileDef[pos].toDouble(&ok);

    //    qDebug() << "zone->endComp[i]" << pos << zone->endComp[pos] << swatreProfileDef[pos];

        if (!ok)
            Error(QString("SWATRE: Can't read compartment end of node %1").arg(pos));
        if (zone->endComp[i] <= 0)
            Error(QString("SWATRE: compartment end of node nr. %1 <= 0").arg(pos));

        /* compute dz and make negative */
        zone->dz[pos] = (zone->endComp[pos-1]-zone->endComp[pos]);
        zone->z[pos] = zone->z[pos-1] + 0.5*(zone->dz[pos-1]+zone->dz[pos]);
        zone->disnod[pos] = zone->z[pos] - zone->z[pos-1];
    }
    zone->dz[0] = -zone->endComp[1];
    zone->z[0] = zone->dz[1]*0.5;
    zone->disnod[0] = zone->z[1];

    zone->disnod[zone->nrNodes] = 0.5 * zone->dz[zone->nrNodes-1];

 //   for (int i = 0; i <= zone->nrNodes; i++)
    //    qDebug() << i << "dz" << zone->dz[i] << "z" << zone->z[i] << "dist" << zone->disnod[i];


    //  count and check valid profiles
    QStringList checkList; // temp list to check for double profile nrs
    for (int i = zone->nrNodes+1; i < swatreProfileDef.count(); i++) {
        if (swatreProfileDef[i].contains("###")) {
            checkList << swatreProfileDef[i+1];
            nrProfileList++;
        }
    }
    sizeProfileList = nrProfileList;
qDebug() << "nr profiles" << nrProfileList << checkList.count();

    if (nrProfileList == 0)
        Error(QString("SWATRE: no profiles read from %1").arg(SwatreTableName));

    checkList.sort();
    for (int i = 0; i < checkList.count()-1; i++)
    {
        if (checkList[i] == checkList[i+1])
            DEBUG(QString("Warning SWATRE: profile id %1 defined more than once").arg(checkList[i+1]));
    }

    // ordered int list with profile nrs
    // makes checklist redundant
    swatreProfileNr.clear();
    for (int i = 0; i < checkList.count(); i++)
        swatreProfileNr << checkList[i].toInt();

    std::sort(swatreProfileNr.begin(), swatreProfileNr.end());

    // qDebug() << swatreProfileNr;
    // qDebug() << checkList;

    profileList = (PROFILE **)realloc(profileList,sizeof(PROFILE *)*(nrProfileList+1));
    // why realloc? profile list is a list of pointers to PROFILE

    nrProfileList = 0;
    for (int i = zone->nrNodes+1; i < swatreProfileDef.count(); i++) {
        if (swatreProfileDef[i].contains("###")) {
            i++;
            profileList[nrProfileList] = ReadProfileDefinitionNew(i, zone);
            // creates a profile and gives the pointer to profilelist
            // i is the place in the StrinList where a profile starts
            nrProfileList++;
        }
    }

}
//----------------------------------------------------------------------------------------------
/// read a new profile from profile.inp, construct and return it
/**
* returns pointer to each new profile in file profile.inp\n
* for info:
\code
 typedef struct PROFILE {
   int            profileId; 	// number identifying this profile  >= 0
   const ZONE     *zone; 		// array with zone.nrNodes elements: dz, z etc
   const HORIZON  **horizon; 	// ptr to horizon information this node belongs to
 } PROFILE;
\endcode
*/
PROFILE * TWorld::ReadProfileDefinitionNew(int pos, ZONE *z)
{
    QString tableName;
    double endHor = 0, endHorPrev = 0;
    PROFILE *p;
    HORIZON *h;
    bool ok;
//qDebug() << pos << "readprofdefnew" << z->nrNodes;

    // profile has a pointer to LUT
    // allocate profile memory, PROFILE is defined in swatre_p.h
    p = new PROFILE;

    p->profileId = swatreProfileDef[pos].toInt(&ok, 10);
    if (!ok)
        Error(QString("SWATRE: read error: error in profile id %1 definition").arg(p->profileId));

    p->horizon = (const HORIZON **)malloc(sizeof(HORIZON *) * z->nrNodes);
    p->zone = z; // pointer

    int i = 0;
    while (i != z->nrNodes) {
        pos++; // profile nr, move one line

        tableName = swatreProfileDef[pos];

        if (!QFileInfo(SwatreTableDir + tableName).exists())
            Error(QString("SWATRE: Can't read the LUT for profile nr %1 node nr %2 and up").arg(p->profileId).arg(i+1));

        endHorPrev = endHor;
        pos++;
        endHor = swatreProfileDef[pos].toDouble(&ok);
        if (!ok)
            Error(QString("SWATRE: Can't read end of horizon for profile nr %1").arg(p->profileId));
        if (endHor <= endHorPrev)
            Error(QString("SWATRE: Error in profile definition nr %1, depth horizons do not increase").arg(p->profileId));

        h = ReadHorizonNew(SwatreTableDir, tableName);

        // copy horizon info to all nodes of this horizon
        while (i < z->nrNodes && z->endComp[i] <= endHor )
            p->horizon[i++] = h;

     //   if (z->endComp[i-1] != endHor)
       //    Error(QString("SWATRE: Compartment does not end on depth '%1' (found in profile nr %2 for horizon %3)")
         //        .arg(endHor).arg(p->profileId).arg(tableName));
        //????????????? what does this error mean exactly?
    }

    return(p); // return the profile
}
//----------------------------------------------------------------------------------------------
/// OBSOLETE, replaced by ReadSwatreInputNew
/** make profile list and read all profile data */
//int TWorld::ReadSwatreInput(QString fileName, QString tablePath)
//{
// // PROFILE.INP is opened here
// FILE *f = fopen(fileName.toLatin1().constData(), "r");
// if (f == nullptr)
// {
//    Error(QString("SWATRE: Can't open profile definition file %1").arg(fileName));
//    throw 1;
// }
// //All file name checking in main program

// // set profiles to nullptr
// InitializeProfile();

// // read node distances and make structure
// ZONE *z = ReadNodeDefinition(f);

// // check if list can hold new one
// do {

//    if (nrProfileList == sizeProfileList)
//    {
//       int i = sizeProfileList;
//       sizeProfileList += LIST_INC;
//       profileList = (PROFILE **)realloc(profileList,sizeof(PROFILE *)*sizeProfileList);
//       while (i < sizeProfileList)
//          profileList[i++] = nullptr;
//    }
//    profileList[nrProfileList] =
//          ReadProfileDefinition(f,z,tablePath.toLatin1().constData());

// } while (profileList[nrProfileList++] != nullptr);
// // correct for eof-marker and test if something is read
// if ( --nrProfileList == 0)
//    Error(QString("SWATRE: no profiles read from %1").arg(fileName));


// // make profileList index match the profileId's
// //WHY???? is a profile is numbered 100 I make 100 empty profiles???
// // becaus ethe map has the profile nrs, but that can be solved with a code list!
// int mmax = 0;
// for (int i = 0 ; i < nrProfileList; i++) {
//      mmax = std::max(mmax, profileList[i]->profileId);
//      mmax++;
// }

// PROFILE **tmpList = (PROFILE **)malloc(mmax*sizeof(PROFILE *));
// for (int i = 0 ; i < mmax; i++)
//    tmpList[i] = nullptr;

// for (int i = 0 ; i < nrProfileList; i++)
//    if (tmpList[profileList[i]->profileId] == nullptr)
//       tmpList[profileList[i]->profileId] = profileList[i];
//    else
//       Error(QString("SWATRE: profile with id '%!' declared more than once").arg(profileList[i]->profileId));

// free(profileList);

// profileList = tmpList;
// nrProfileList = mmax;
// sizeProfileList = mmax;

// /* PROFILE.INP is closed here */
// fclose(f);

// return (0);
//}
//----------------------------------------------------------------------------------------------
/// OBSOLETE, no longer used
/** return PROFILE or throw an error if not found, profileNr is the profile map value
    only used in Swatinp and no longer necessary
*/
PROFILE *TWorld::ProfileNr(int profileNr)
{
    if (profileNr < 0 || profileNr >= nrProfileList)
        return(nullptr);
    return(profileList[profileNr]);
}
//----------------------------------------------------------------------------------------------
void  TWorld::FreeSwatreInfo(void)
{
    int i;

    // if (zone == nullptr)
    //    return;

    // if (zone != nullptr) {
    //     zone->dz.clear();
    //     zone->z.clear();
    //     zone->endComp.clear();
    //     zone->disnod.clear();
    //     delete(zone);
    //     zone = nullptr;
    // }

    // // if (zone != nullptr) {
    // //     free(zone->dz);
    // //     free(zone->z);
    // //     free(zone->endComp);
    // //     free(zone->disnod);
    // //     free(zone);
    // //     zone = nullptr;
    // // }

    // if (profileList != nullptr) {
    //     if (profileList[0] != nullptr) {
    //         for(int i=0; i < nrProfileList; i++)
    //             if (profileList[i] != nullptr)
    //                 free(profileList[i]);
    //     }
    //     free(profileList);
    //     profileList = nullptr;
    // }
    // // for(i=0; i < nrProfileList; i++)
    // //    if (profileList[i] != nullptr)
    // //       free(profileList[i]);
    // // free(profileList);
    // // profileList = nullptr;
     nrProfileList = 0;
     sizeProfileList = 0;

    // if (horizonList != nullptr) {
    //     for(i=0; i < nrHorizonList; i++)
    //     {
    //         //free(horizonList[i]->name);
    //         //FreeLut(horizonList[i]->lut);
    //         horizonList[i]->lut = nullptr;
    //         free(horizonList[i]);
    //     }
    //     free(horizonList);
    //     horizonList = nullptr;
    // }

    nrHorizonList = 0;
    sizeHorizonList = 0;
}
//----------------------------------------------------------------------------------------------
/// copy horizon info to all nodes of this horizon
HORIZON * TWorld::ReadHorizonNew(QString tablePath, QString tableName)
{
    //qDebug() << "new hor";
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
/// copy horizon info to all nodes of this horizon
//HORIZON * TWorld::ReadHorizon(QString const char *tablePath,	const char *tableName)
HORIZON * TWorld::ReadHorizon(QString tablePath, QString tableName)
{
return nullptr;
}
//----------------------------------------------------------------------------------------------
