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

/*!
  \file LisUIMapnames.cpp
  \brief fill map list structure and handle map tree

 * fill the map list structure with names and basic descriptions
 * enable and disable braches according to users choices
 * get and set mapames when direct editing by the user, edit mapList
 */


#include "lisemqt.h"


//--------------------------------------------------------------------
// put the mapnames in the namelist structure, e.g. after screen update
// if to = true then maplist to namelist, else namelist to maplist
void lisemqt::fillNamelistMapnames(bool to)
{
   for (int j = mapstartnr; j < nrnamelist; j++)
   {
      if (!namelist[j].name.startsWith("[") && !namelist[j].name.isEmpty())
      {
         for (int k = 0; k < nrmaplist; k++)
         {
            if (mapList[k].name.toUpper() == namelist[j].name.toUpper())
            {
               if (to)
               {
                  namelist[j].value = mapList[k].value;
                  // update namelist (the run file) with the current map names
               }
               else
                  mapList[k].value = namelist[j].value;
               // update the maplist with the namelist data when the runfile is parsed
               break;
            }
         }
      }
   }
   //   for (int k = 0; k < nrmaplist; k++)
   //      qDebug() << mapList[k].name << mapList[k].value << mapList[k].groupnr << mapList[k].varnr;
}
//--------------------------------------------------------------------
// DEFmaps has default var names, filenames and descriptions, in LisUIDefaultNames.cpp
// first number 0/1/2 is flag if title treenode or subnode
// fill the mapList structure
void lisemqt::fillMapnames()
{
   int subbranch = 0, branch = -1, nr = -1 /*VJ bug fix */, dec = 0;
   QStringList SL;

   DefaultMapnames();
   //VJ 110113 take the default map list

   for (int i = 0; i < DEFmaps.count(); i++)\
   {

      if (DEFmaps[i].startsWith("0"))
      {
         branch++;
         subbranch = 0;
         dec = 0;
      }
      if (DEFmaps[i].startsWith("1"))
      {
         dec += 10;
         subbranch = 0;
      }
      if (DEFmaps[i].startsWith("2"))
      {
         nr++;

         SL = DEFmaps[i].split(";",QString::SkipEmptyParts);
         mapList[nr].groupnr=branch;
         mapList[nr].varnr=subbranch+dec;
         mapList[nr].value=SL[2];   //<= mapname
         mapList[nr].name=SL[4];    //<= id string to recognize variable
         mapList[nr].dir="";//not used for now

         subbranch++; //VJ 110326 moved to here, branch numbers were wrong
      }
   }
   nrmaplist = nr+1; //VJ bug fix
}
//--------------------------------------------------------------------
/** enables or disables a branch and expands or contracts it */
void lisemqt::checkMapNameModel(int parentrow, int selrow, bool setit)
{
	if (MapNameModel)
	{
		QModelIndex indexParent = MapNameModel->index(parentrow, 0);
		// select and expand chosen main and sub-branch
		if (selrow >= 0)
		{
			MapNameModel->setFlag(setit, parentrow);
			// set top level node
			if (selrow == 0)
				for (int k = 0; k < MapNameModel->rowCount(indexParent); k++)
					MapNameModel->setFlag(setit, k, indexParent);
			else
				if (selrow >= 10)
				{
					selrow -= 10;
					MapNameModel->setFlag(setit, selrow, indexParent);

					QModelIndex indexChild = MapNameModel->index(selrow, 0, indexParent);

					for (int k = 0; k < MapNameModel->rowCount(indexChild); k++)
						MapNameModel->setFlag(setit, k, indexChild);
				}

			if (setit)
			{
				treeView->expand(MapNameModel->index(parentrow,0));
				if (selrow >= 0)
					treeView->expand(MapNameModel->index(selrow, 0, indexParent));
			}
			else
			{
				if (selrow >= 0)
					treeView->collapse(MapNameModel->index(selrow, 0, indexParent));
				if (selrow == 0)
					treeView->collapse(MapNameModel->index(parentrow,0));
			}
		}
	}
}
//--------------------------------------------------------------------
/** initialize the map tree interface, this function is called twice: first to initialize the interface
also after each call of a runfile so that the runfile mapnames are loaded */
void lisemqt::initMapTree()
{
	if (MapNameModel)
	{
		delete MapNameModel;
		MapNameModel = NULL;
	}

   MapNameModel = new TreeModel(DEFmaps);
   // DEFmaps is used to construct the interface tree struture

	treeView->setModel(MapNameModel);

	treeView->setColumnWidth(0,196);
	treeView->setColumnWidth(1,196);
	treeView->setColumnWidth(2,400);
	treeView->setColumnWidth(3,0);
	treeView->setColumnWidth(4,0);
	treeView->setAlternatingRowColors(true);

   checkMapNameModel(RAINFALLMAPS,0, true);
   checkMapNameModel(CATCHMENTMAPS,0, true);
   checkMapNameModel(LANDUSEMAPS,0, true);
   checkMapNameModel(SURFACEMAPS,0, true);
   checkMapNameModel(EROSIONMAPS,0, true);
	// enable basic maps, tree nodes 0-4 = first 5 branches
   
	treeView->collapseAll();
}
//--------------------------------------------------------------------
/** edit mapname in reponse to edit keys (like F2)
 get user edited name from the tree and put it in mapList structure for saving
 */
void lisemqt::editMapname(QModelIndex topLeft, QModelIndex bottomRight )
{
   int groupnr = topLeft.parent().parent().row(); // assume 3 level structure
   int varnr = topLeft.row();

   // if not 3 level structure then use 2nd level
   if (topLeft.parent().parent().row() < 0)
      groupnr = topLeft.parent().row();

   if (groupnr == INFILTRATIONMAPS || groupnr == CHANNELSMAPS || groupnr == NUTRIENTSMAPS)
      varnr = (topLeft.parent().row()+1)*10 + topLeft.row();

   for (int k = 0; k < nrmaplist; k++)
   {
      if (mapList[k].groupnr == groupnr && mapList[k].varnr == varnr)
      {
         QVariant d = MapNameModel->data(topLeft,Qt::DisplayRole);;//MapNameModel->data(MapNameModel->index(j, k, indexParent),0);
         mapList[k].value = d.toString(); //<== put new map name
      }
   }
}
//--------------------------------------------------------------------
/** open a openfile window when double clicked on a treeView mapname.
 Look for the correct map in mapList and copy the new file name into mapList */
void lisemqt::openMapname(QModelIndex topLeft)
{
   if (topLeft.column() != 1 || (topLeft.row() == 0 && topLeft.parent().row() < 0))
      return;
   // if not a mapname return

   // else find out which mapname it is

   int groupnr = topLeft.parent().parent().row();
   int k, varnr = topLeft.row();
   // assume 3 level structure

   // if not 3 level structure then use 2nd level
   if (topLeft.parent().parent().row() < 0)
      groupnr = topLeft.parent().row();

   if (groupnr == INFILTRATIONMAPS || groupnr == CHANNELSMAPS || groupnr == NUTRIENTSMAPS)
      varnr = (topLeft.parent().row()+1)*10 + topLeft.row();
   // correct for 3 level structures

   // find the map in mapList
   for (k = 0; k < nrmaplist; k++)
   {
      if (mapList[k].groupnr == groupnr && mapList[k].varnr == varnr)
      {
         QVariant d = MapNameModel->data(topLeft,Qt::DisplayRole);
         break;
      }
   }

   QString path = QFileDialog::getOpenFileName(this,	QString("Select the map: %1;")
                                               .arg(mapList[k].value),E_MapDir->text(),QString("*.map *.csf;;*.*"));
   // open file dialog


   if (!path.isEmpty())// && QFileInfo(path).exists())
   {
      MAP *m = Mopen(QFileInfo(path).absoluteFilePath().toAscii().constData(), M_READ);
      if (m == NULL)
      {
         QMessageBox::critical(this, "openLISEM",
                               QString("File \"%1\" is not a PCRaster map.")
                               .arg(path));
         return;
      }

      mapList[k].value = QFileInfo(path).fileName();
      mapList[k].dir = QFileInfo(path).dir().path();
      // put the name and path into the mapList structure

      //qDebug() << "mapname edit" <<  mapList[k].name << mapList[k].id << k;

      QVariant d(mapList[k].value);
      MapNameModel->setData(topLeft, d, Qt::EditRole);
      // put the name into the treeView and MapNameModel
   }
}
//--------------------------------------------------------------------
