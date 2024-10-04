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
- obsolete:
- int TWorld::ReadSwatreInput(QString fileName, QString tablePath) \n
- ZONE * TWorld::ReadNodeDefinition(FILE *f) \n
- PROFILE * TWorld::ReadProfileDefinition(FILE *f,ZONE *z,const char *tablePath) \n
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
   profileList = nullptr;
   nrProfileList=0;
   sizeProfileList=0;
   zone=nullptr;

   horizonList = nullptr;
   nrHorizonList=0;
   sizeHorizonList=0;

   swatreProfileDef.clear();
   swatreProfileNr.clear();
}
//----------------------------------------------------------------------------------------------
/// read and parse profile.inp
/// new version using swatreProfileDef QStringList
void TWorld::ReadSwatreInputNew(void)
{
    // get the profile inp file contents and put in stringlist
    QFile fin(SwatreTableName);
    if (!fin.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Error(QString("SWATRE: Can't open profile definition file %1").arg(SwatreTableName));
        throw 1;
    }

    // set profiles to nullptr
    InitializeProfile();

    // read profile.inp and put in StringList, trimmed and without blank lines
    while (!fin.atEnd())
    {
        QString S = fin.readLine();
        S = S.trimmed();
        if (!S.isEmpty())
        {
            if (S.indexOf("#")< 0)
                swatreProfileDef << S.trimmed();
            else
            {
                int pos = S.indexOf("#");
                if (pos > 0)
                    swatreProfileDef << S.remove(pos,S.size()).trimmed();
            }
        }
    }
    fin.close();

    for (int i = 0; i < swatreProfileDef.count(); i++) {
        qDebug() << swatreProfileDef[i];
    }

    ZONE *z;
    int  i;
    QStringList checkList; // temp list to check for double profile nrs


    // read node distances and make structure
    z = ReadNodeDefinitionNew();


    //   mark count nr profiles
    for (int i = z->nrNodes+1; i < swatreProfileDef.count(); i++) {
        bool ok = false, oki = false;
        double dummy = swatreProfileDef[i-1].toDouble(&ok);
        int dummi = swatreProfileDef[i].toInt(&oki);
        // if(swatreProfileDef[i-1].toDouble() && swatreProfileDef[i].toInt())
        if (ok && oki)
        {
            // if two values follow each other this is a new profile
            nrProfileList++;
            checkList << swatreProfileDef[i];
        }
    }
    sizeProfileList = nrProfileList;

    if (nrProfileList == 0)
        Error(QString("SWATRE: no profiles read from %1").arg(SwatreTableName));

    checkList.sort();
    for (int i = 0; i < checkList.count()-1; i++)
    {
        if (checkList[i] == checkList[i+1])
            DEBUG(QString("Warning SWATRE: profile id %1 declared more than once").arg(checkList[i+1]));
    }
    // alternative duplicate check but with no info on duplicates
    //   if (checkList.removeDuplicates() > 0)
    //      Error(QString("SWATRE: there are duplicate profile ID's.");

    // ordered int list with profile nrs
    // makes checklist redundant
    swatreProfileNr.clear();
    for (int i = 0; i < checkList.count(); i++)
        swatreProfileNr << checkList[i].toInt();
    std::sort(swatreProfileNr.begin(), swatreProfileNr.end());

    profileList = (PROFILE **)realloc(profileList,sizeof(PROFILE *)*sizeProfileList+1);

    nrProfileList = 0;
    for (i = z->nrNodes+1; i < swatreProfileDef.count(); i++) {
        bool ok = false, oki = false;
        double dummy = swatreProfileDef[i-1].toDouble(&ok);
        int dummi = swatreProfileDef[i].toInt(&oki);
        if (ok && oki)
            //  if(swatreProfileDef[i-1].toDouble() && swatreProfileDef[i].toInt())
        {
            profileList[nrProfileList] = ReadProfileDefinitionNew(i, z);
            // read the profile tables for each defined prpfile,
            // i is the place in the StrinList where a profile starts
            nrProfileList++;
        }
    }

}
//----------------------------------------------------------------------------------------------
/** allocates ZONE structure,
 *   reads compartment ends:
 *   2.5 5 10 means dz[0] = dz[1] = 2.5, dz[2] = 5, etc.\n
 *   computes all parameters stored in ZONE -structure
 */

ZONE * TWorld::ReadNodeDefinitionNew(void)
{
   int  i;
   bool ok;
   int pos = 0;

   zone = (ZONE *)malloc(sizeof(ZONE));
   zone->nrNodes = swatreProfileDef[pos].toInt(&ok, 10);
   pos++;
   if (!ok)
      Error(QString("SWATRE: Can't read number of nodes %1 from input file").arg(zone->nrNodes));
   if (zone->nrNodes < 1 )
      Error(QString("SWATRE: number of nodes %1 smaller than 1").arg(zone->nrNodes));
   if (zone->nrNodes > MAX_NODES)
      Error(QString("SWATRE: number of nodes %1 larger than %2").arg(zone->nrNodes).arg(MAX_NODES));

   zone->dz     = (double *)malloc(sizeof(double)*zone->nrNodes);
   zone->z      = (double *)malloc(sizeof(double)*zone->nrNodes);
   zone->disnod = (double *)malloc(sizeof(double)*(zone->nrNodes+1));
   zone->endComp= (double *)malloc(sizeof(double)*zone->nrNodes);


   for (i=0; i < zone->nrNodes; i++)
   {
      zone->endComp[i] = swatreProfileDef[pos+i].toDouble(&ok);
      if (!ok)
         Error(QString("SWATRE: Can't read compartment end of node %1").arg(i+1));
      if (zone->endComp[i] <= 0)
         Error(QString("SWATRE: compartment end of node nr. %1 <= 0").arg(i+1));

      /* compute dz and make negative */
      zone->dz[i]= ( (i == 0) ? -zone->endComp[0] : (zone->endComp[i-1]-zone->endComp[i]));
      zone->z[i]= ( (i == 0) ? zone->dz[i]*0.5 : zone->z[i-1] + 0.5*(zone->dz[i-1]+zone->dz[i]));
      zone->disnod[i] = ( (i == 0) ? zone->z[i]: zone->z[i] - zone->z[i-1]);

      //qDebug() << i << "dz" << zone->dz[i] << "z" << zone->z[i] << "dist" << zone->disnod[i];
   }
   zone->disnod[zone->nrNodes] = 0.5 * zone->dz[zone->nrNodes-1];

   return(zone);
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
PROFILE * TWorld::ReadProfileDefinitionNew(
      int pos,
      ZONE *z)         /* zone division this profile */
{
   QString tableName;
   int  i;
   double endHor = 0, endHorPrev = 0;
   PROFILE *p;
   HORIZON *h;
   bool ok;


   /* profile has a pointer to LUT */
   // allocate profile memory, PROFILE is defined in swatre_p.h
   p = (PROFILE *)malloc(sizeof(PROFILE));
   p->profileId = swatreProfileDef[pos].toInt(&ok, 10);
   if (!ok)
      Error(QString("SWATRE: read error: error in profile id %1 definition").arg(p->profileId));

   p->horizon = (const HORIZON **)malloc(sizeof(HORIZON *)*z->nrNodes);
   p->zone = z;
   p->nrNodes = z->nrNodes;

   //pos++; // lastline of prev profile

   i = 0;
   while (i != z->nrNodes)
   {
      pos++; // profile nr

      // read tablename from swatreProfileDef (= profile.inp without blanks)
      tableName = swatreProfileDef[pos];
      if (!QFileInfo(SwatreTableDir + tableName).exists())
         Error(QString("SWATRE: Can't read a LUT for profile nr %1 node nr %2 and up").arg(p->profileId).arg(i+1));

      endHorPrev = endHor;
      pos++; // profile depth, endHor
      endHor = swatreProfileDef[pos].toDouble(&ok);//toInt(&ok, 10);
      if (!ok)
         Error(QString("SWATRE: Can't read end of horizon for profile nr %1").arg(p->profileId));

      if (endHor <= endHorPrev)
         Error(QString("SWATRE: Error in profile definition nr %1").arg(p->profileId));

      h = ReadHorizon(SwatreTableDir.toLatin1().constData(), tableName.toLatin1().constData());
      // copy horizon info to all nodes of this horizon

      while (i < z->nrNodes && z->endComp[i] <= endHor )
         p->horizon[i++] = h;

      if (z->endComp[i-1] != endHor)
         Error(QString("SWATRE: No compartment ends on depth '%1' (found in profile nr %2 for horizon %3)")
               .arg(endHor).arg(p->profileId).arg(tableName));
   }
   return(p);
}
//----------------------------------------------------------------------------------------------
/// OBSOLETE, replaced by ReadSwatreInputNew
/** make profile list and read all profile data */
int TWorld::ReadSwatreInput(QString fileName, QString tablePath)
{
   FILE *f;
   ZONE *z;
   int  i, mmax;
   PROFILE **tmpList;
   // PROFILE.INP is opened here
   f = fopen(fileName.toLatin1().constData(), "r");

   if (f == nullptr)
   {
      Error(QString("SWATRE: Can't open profile definition file %1").arg(fileName));
      throw 1;
   }
   //All file name checking in main program

   // set profiles to nullptr
   InitializeProfile();

   // read node distances and make structure
   z = ReadNodeDefinition(f);

   // check if list can hold new one
   do {

      if (nrProfileList == sizeProfileList)
      {
         int i = sizeProfileList;
         sizeProfileList += LIST_INC;
         profileList = (PROFILE **)realloc(profileList,sizeof(PROFILE *)*sizeProfileList);
         while (i < sizeProfileList)
            profileList[i++] = nullptr;
      }
      profileList[nrProfileList] =
            ReadProfileDefinition(f,z,tablePath.toLatin1().constData());

   } while (profileList[nrProfileList++] != nullptr);
   // correct for eof-marker and test if something is read
   if ( --nrProfileList == 0)
      Error(QString("SWATRE: no profiles read from %1").arg(fileName));

   /* make profileList index match the profileId's */
   mmax = 0;
   for (i = 0 ; i < nrProfileList; i++)
      mmax = std::max(mmax, profileList[i]->profileId);
   mmax++;

   tmpList = (PROFILE **)malloc(mmax*sizeof(PROFILE *));
   for (i = 0 ; i < mmax; i++)
      tmpList[i] = nullptr;

   for (i = 0 ; i < nrProfileList; i++)
      if (tmpList[profileList[i]->profileId] == nullptr)
         tmpList[profileList[i]->profileId] = profileList[i];
      else
         Error(QString("SWATRE: profile with id '%!' declared more than once").arg(profileList[i]->profileId));

   free(profileList);

   profileList = tmpList;
   nrProfileList = mmax;
   sizeProfileList = mmax;

   /* PROFILE.INP is closed here */
   fclose(f);

   return (0);
}
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

   /* currently, all profiles have the same zoning */
   if (zone == nullptr)
      return;

   free(zone->dz);
   free(zone->z);
   free(zone->endComp);
   free(zone->disnod);
   free(zone);

   for(i=0; i < nrProfileList; i++)
      if (profileList[i] != nullptr)
         free(profileList[i]);
   free(profileList);
   profileList = nullptr;
   nrProfileList=sizeProfileList=0;

   for(i=0; i < nrHorizonList; i++)
   {
      free(horizonList[i]->name);
      FreeLut(horizonList[i]->lut);
      horizonList[i]->lut = nullptr;
      free(horizonList[i]);
   }
   free(horizonList);
   horizonList = nullptr;
   nrHorizonList=sizeHorizonList=0;
}
//----------------------------------------------------------------------------------------------
/** allocates ZONE structure,
    reads compartment ends:
    2.5 5 10 means dz[0] = dz[1] = 2.5, dz[2] = 5, etc.\n
    computes all parameters stored in ZONE -structure
*/
ZONE * TWorld::ReadNodeDefinition(FILE *f)
{
   int  i;
   zone = (ZONE *)malloc(sizeof(ZONE));
   if ( fscanf(f,"%d",&(zone->nrNodes)) != 1 )
      Error(QString("SWATRE: Can't read number of nodes %1 from input file").arg(zone->nrNodes));
   if (zone->nrNodes < 1 )
      Error(QString("SWATRE: number of nodes %1 smaller than 1").arg(zone->nrNodes));
   if (zone->nrNodes > MAX_NODES)
      Error(QString("SWATRE: number of nodes %1 larger than %2").arg(zone->nrNodes).arg(MAX_NODES));

   zone->dz     = (double *)malloc(sizeof(double)*zone->nrNodes);
   zone->z      = (double *)malloc(sizeof(double)*zone->nrNodes);
   zone->disnod = (double *)malloc(sizeof(double)*(zone->nrNodes+1));
   zone->endComp= (double *)malloc(sizeof(double)*zone->nrNodes);

   for (i=0; i < zone->nrNodes; i++)
   {
      if ( fscanf(f,"%lf",&(zone->endComp[i])) != 1 )
         Error(QString("SWATRE: Can't read compartment end of node %1").arg(i+1));
      if (zone->endComp[i] <= 0)
         Error(QString("SWATRE: compartment end of node nr. %1 <= 0").arg(i+1));
      /* compute dz and make negative */
      zone->dz[i]= ( (i == 0) ? -zone->endComp[0] : (zone->endComp[i-1]-zone->endComp[i]));
      zone->z[i]= ( (i == 0) ? zone->dz[i]*0.5 : zone->z[i-1] + 0.5*(zone->dz[i-1]+zone->dz[i]));
      zone->disnod[i] = ( (i == 0) ? zone->z[i]: zone->z[i] - zone->z[i-1]);

      //qDebug() << i << "dz" << zone->dz[i] << "z" << zone->z[i] << "dist" << zone->disnod[i];
   }
   zone->disnod[zone->nrNodes] = 0.5 * zone->dz[zone->nrNodes-1];

   return(zone);
}
//----------------------------------------------------------------------------------------------
/// OBSOLETE, replaced by ReadProfileDefinitionNew
/**
   returns pointer to new profile or nullptr if eof is encountered
   while reading first token of profile definition
*/
PROFILE * TWorld::ReadProfileDefinition(
      FILE *f,
      ZONE *z,         /* zone division this profile */
      const char *tablePath) /* pathName ended with a '/' */
{
   char tableName[256];
   int  i;
   double endHor;
   PROFILE *p;
   HORIZON *h;
   /* profile has a pointer to LUT */
   p = (PROFILE *)malloc(sizeof(PROFILE));


   if ( fscanf(f,"%d",&(p->profileId)) != 1 )
   {
      if (feof(f))
      {
         free(p);
         return nullptr;
      }
      Error(QString("SWATRE: read error: can't read profile id %1").arg(p->profileId));
   }
   if (p->profileId < 0)
      Error(QString("SWATRE: profile id smaller that 0: %1").arg(p->profileId));

   p->horizon = (const HORIZON **)malloc(sizeof(HORIZON *)*z->nrNodes);
   p->zone = z;

   i = 0;
   while (i != z->nrNodes)
   {
      if ( fscanf(f,"%s",tableName) != 1 )
         Error(QString("SWATRE: Can't read a LUT for profile nr %1 node nr %2 and up").arg(p->profileId).arg(i+1));
      if ( fscanf(f,"%lf", &endHor) != 1 )
         Error(QString("SWATRE: Can't read end of horizon for profile nr %1").arg(p->profileId));

      h = ReadHorizon(tablePath, tableName);
      // copy horizon info to all nodes of this horizon

      while (i < z->nrNodes && z->endComp[i] <= endHor )
         p->horizon[i++] = h;
      if (z->endComp[i-1] != endHor)
         Error(QString("SWATRE: No compartment ends on depth '%1' (found in profile nr %2 for horizon %3)")
               .arg(endHor).arg(p->profileId).arg(tableName));
   }
   return(p);
}
//----------------------------------------------------------------------------------------------
/// copy horizon info to all nodes of this horizon
HORIZON * TWorld::ReadHorizon(const char *tablePath,	const char *tableName)
{
   HORIZON	*h;
   char fileName[256];
   double *t, *lutCont;
   int i, nrRowsa=0;

   // look if it's already loaded
   for( i= 0; i < nrHorizonList; i++)
      if (!strcmp(tableName, horizonList[i]->name))
         return(horizonList[i]);

   /* if not then add one */
   /* check for space in list */
   if (nrHorizonList == sizeHorizonList)
   {
      sizeHorizonList += LIST_INC;
      horizonList = (HORIZON **)realloc(horizonList, sizeof(HORIZON *)*sizeHorizonList);
   }

   h = (HORIZON *)malloc(sizeof(HORIZON));
   horizonList[nrHorizonList++] = h;
   strcat(strcpy(fileName, tablePath), tableName);


   /* hook up table to temp array t */
   t = ReadSoilTable(fileName, &nrRowsa);

   lutCont = (double *)malloc(sizeof(double)*NR_COL*(nrRowsa+2));
   for(i=0; i < nrRowsa; i++)
   {
      lutCont[IND(i,THETA_COL)] =  t[i*3+THETA_COL];
      lutCont[IND(i,H_COL)]     =  t[i*3+H_COL];
      lutCont[IND(i,K_COL)]     =  t[i*3+K_COL];
   }

   // check some stuff
   for(i=0; i < (nrRowsa-1); i++)
   {
      if (lutCont[IND(i+1,H_COL)] <= lutCont[IND(i,H_COL)])
         Error(QString("matrix head not decreasing in table %1 at h = %2.").arg(tableName).arg(lutCont[IND(i,H_COL)]));
      if (lutCont[IND(i+1,THETA_COL)] <= lutCont[IND(i,THETA_COL)])
         Error(QString("moisture content not decreasing in table %1 at theta = %2.").arg(tableName).arg(lutCont[IND(i,THETA_COL)]));
   }

   for(i=0; i < (nrRowsa-1); i++)
   {
      lutCont[IND(i,DMCH_COL)] = 0.5 *
            (lutCont[IND(i+1,H_COL)] + lutCont[IND(i,H_COL)]);
      // dif moisture cap dh

      lutCont[IND(i,DMCC_COL)] =
            (lutCont[IND(i+1,THETA_COL)] - lutCont[IND(i,THETA_COL)])/
            (lutCont[IND(i+1,H_COL)] - lutCont[IND(i,H_COL)]);
      // dif moisture cap dtheta/dh
      //qDebug() << i << lutCont[IND(i+1,H_COL)]<< lutCont[IND(i,THETA_COL)] << lutCont[IND(i,DMCH_COL)]<<lutCont[IND(i,DMCC_COL)];
   }
   lutCont[IND(nrRowsa-1,DMCH_COL)] = 0;
   lutCont[IND(nrRowsa-1,DMCC_COL)] = lutCont[IND(nrRowsa-2,DMCC_COL)] ;

   free(t);

   h->name = strcpy((char *)malloc(strlen(tableName)+1), tableName);
   h->lut = CreateLutFromContents(lutCont, true, nrRowsa, LUT_COLS);


   return(h);
}
//----------------------------------------------------------------------------------------------
