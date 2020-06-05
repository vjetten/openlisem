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
\file LisUItreemodel.cpp
\brief treemodel basic funtions. treemodel is the tree structure with input map names
 */

#include "lisemqt.h"

//-----------------------------------------------------------------------------------------------
// is called in FillMapTree() in LisUImapnames.cpp
TreeModel::TreeModel(const QStringList &data, QObject *parent)
   : QAbstractItemModel(parent)
{
   QVector<QVariant> header;

   //VJ 101223 simplified
   header << "Variable name";
   header << "Map name";
   header << "Description";
   header << " ";

   rootItem = new TreeItem(header);

   setupModelData(data, rootItem);
   // make the treestructure, data contains DEFnames

}
//-----------------------------------------------------------------------------------------------
TreeModel::~TreeModel()
{
   delete rootItem;
}
//-----------------------------------------------------------------------------------------------
int TreeModel::columnCount(const QModelIndex & /* parent */) const
{
   return rootItem->columnCount();
}
//-----------------------------------------------------------------------------------------------
QVariant TreeModel::data(const QModelIndex &index, int role) const
{

   if (!index.isValid())
      return QVariant();

   if (role != Qt::DisplayRole && role != Qt::EditRole && role)// != Qt::CheckStateRole)
      return QVariant();

   TreeItem *item = getItem(index);
   // int dt = item->getDatatype();

   // if (dt == 3 && role == Qt::CheckStateRole && index.column() == 1)
   //     return (item->data(1).toInt() != 0) ? Qt::Checked : Qt::Unchecked;

   return item->data(index.column());

}
//-----------------------------------------------------------------------------------------------
Qt::ItemFlags TreeModel::flags(const QModelIndex &index) const
{
   if (!index.isValid())
      return Qt::NoItemFlags;

   Qt::ItemFlags flags = Qt::ItemIsEnabled | Qt::ItemIsSelectable ;

   if (index.column() == 1)
   {
      flags |= Qt::ItemIsEditable;
   }
   TreeItem *item = getItem(index);
   if (!item->nodeEnabled)
        flags = Qt::NoItemFlags;

   return flags;

}
//-----------------------------------------------------------------------------------------------
TreeItem *TreeModel::getItem(const QModelIndex &index) const
{
   if (index.isValid()) {
      TreeItem *item = static_cast<TreeItem*>(index.internalPointer());
      if (item) return item;
   }
   return rootItem;
}
//-----------------------------------------------------------------------------------------------
QVariant TreeModel::headerData(int section, Qt::Orientation orientation,
                               int role) const
{
   if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
      return rootItem->data(section);

   return QVariant();
}
//-----------------------------------------------------------------------------------------------
QModelIndex TreeModel::index(int row, int column, const QModelIndex &parent) const
{

   if (parent.isValid() && parent.column() != 0)
      return QModelIndex();
   TreeItem *parentItem = getItem(parent);

   TreeItem *childItem = parentItem->child(row);
   if (childItem)
      return createIndex(row, column, childItem);
   else
      return QModelIndex();
}
//-----------------------------------------------------------------------------------------------
void TreeModel::setFlag(int f, int row, const QModelIndex &parent) const
{
   TreeItem *parentItem = getItem(parent);
   TreeItem *childItem = parentItem->child(row);
   if (childItem)
      childItem->nodeEnabled = f;
}
//-----------------------------------------------------------------------------------------------
bool TreeModel::getFlag(int row, const QModelIndex &parent) const
{
   TreeItem *parentItem = getItem(parent);
   TreeItem *childItem = parentItem->child(row);
   if (childItem)
      return (childItem->nodeEnabled);
   return (false);
}
//-----------------------------------------------------------------------------------------------
QModelIndex TreeModel::parent(const QModelIndex &index) const
{
   if (!index.isValid())
      return QModelIndex();

   TreeItem *childItem = getItem(index);
   TreeItem *parentItem = childItem->parent();

   if (parentItem == rootItem)
      return QModelIndex();

   return createIndex(parentItem->childNumber(), 0, parentItem);
}
//-----------------------------------------------------------------------------------------------
int TreeModel::rowCount(const QModelIndex &parent) const
{
   TreeItem *parentItem = getItem(parent);

   return parentItem->childCount();
}
//-----------------------------------------------------------------------------------------------
bool TreeModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
   if (role != Qt::EditRole )
      return false;

   TreeItem *item = getItem(index);
   bool result = item->setData(index.column(), value);
   if (result)
   {
      emit dataChanged(index, index);
   }
   return result;
}
//-----------------------------------------------------------------------------------------------
bool TreeModel::setHeaderData(int section, Qt::Orientation orientation,
                              const QVariant &value, int role)
{
   if (role != Qt::EditRole || orientation != Qt::Horizontal)
      return false;

   bool result = rootItem->setData(section, value);

   if (result)
      emit headerDataChanged(orientation, section, section);

   return result;
}
//-----------------------------------------------------------------------------------------------
// fills the tree structure with the QStringList with mapnames and descriptions
// this QStringlist is constructed in DefaultMapnames() in LisUImapnames.cpp
void TreeModel::setupModelData(const QStringList &lines, TreeItem *parent)
{
   QList <TreeItem*> parents;
   QList <int> indentations;
   parents << parent; // add rootitem to parents list
   indentations << 0;

   int number = 0;

   while (number < lines.count())
   {
      int position = 0;
    	// get the lineData=string and parse it to a qlist = columnData
      QString lineData = lines[number].trimmed();

      if (!lineData.isEmpty())
      {
            QStringList columnStrings = lineData.split(";", Qt::KeepEmptyParts);
         QList<QVariant> columnData;

         position = columnStrings[0].toInt();
         if (position == 0 || position == 1)
            columnData << columnStrings[1];
         // string is a title/category field when "position" is 0 or 1
         else
            for (int column = 1; column < columnStrings.count()-1; ++column)
               columnData << columnStrings[column];
         //fill column data, "position" > 1

         // determine level, where is a columnData linked, what is child, what is parent
         // indentations is an integer list
         if (position > indentations.last())
         {
            if (parents.last()->childCount() > 0)
            {
               parents << parents.last()->child(parents.last()->childCount()-1);
               indentations << position;
            }
         }
         else
         {
            while (position < indentations.last() && parents.count() > 0)
            {
               parents.pop_back();  //pop_back means last item in list is removed?
               indentations.pop_back();

            }
         }

         // Append a new item to the current parent's list of children.
         TreeItem *parent = parents.last();
         parent->insertChildren(parent->childCount(), 1, rootItem->columnCount());
         for (int column = 0; column < columnData.size(); ++column)
            parent->child(parent->childCount() - 1)->setData(column, columnData[column]);
      }
      number++;
   }

}
//-----------------------------------------------------------------------------------------------


/* just a check write function
//label->setText("saving ...");
QFile fout("hup.txt");
fout.open(QIODevice::ReadWrite);
for (int i = 0; i < MapNameModel->rowCount(); i++)
{
 QVariant d;
 int _r = MapNameModel->index(i,0).row();
 QString S = "["+QString::number(_r)+"-" + QString::number(0) + "] ";
 for (int j = 0; j < 2; j++)
 {
  d = MapNameModel->data(MapNameModel->index(i,j),0); // parent
  S = S + d.toString()+";";
 }
 S = S + "\n";
 QModelIndex indexParent = MapNameModel->index(i, 0);
 for (int j = 0; j < MapNameModel->rowCount(indexParent); j++)
 {
  int _r = MapNameModel->index(i,j).row();
  int _c = MapNameModel->index(i,j).column();
  S = S + "["+QString::number(_r)+"-" + QString::number(_c) + "] ";

  for (int k = 0; k < MapNameModel->columnCount(indexParent); k++)
  {
   int _rr = MapNameModel->index(j, k, indexParent).row();
   int _cc = MapNameModel->index(j, k, indexParent).column();
   d = MapNameModel->data(MapNameModel->index(j, k, indexParent),0);
   S = S + "["+QString::number(_rr)+"-" + QString::number(_cc) + "] ";
   S = S + d.toString()+";";
  }
  S = S + "\n";
 }
 S = S + "\n";
 QByteArray line(S.toLatin1());
 fout.write(line);
}
fout.close();
 */
