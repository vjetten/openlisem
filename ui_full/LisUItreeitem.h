/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

#ifndef TREEITEM_H
#define TREEITEM_H

#include <QtGui>

class TreeItem
{
public:
    TreeItem(const QVector<QVariant> &data, TreeItem *parent = 0);
    ~TreeItem();

    TreeItem *child(int number);
    int childCount() const;
    int columnCount() const;
    bool insertChildren(int position, int count, int columns);
    QVariant data(int column) const;
    TreeItem *parent();
    int childNumber() const;
    bool setData(int column, const QVariant &value);
    bool _flag;  // for setting flags to disable node
 //   int _datatype; //0,1 = node level; 2 = string; 3 = int; 4 = float
 //   void setDatatype(int datatype);
 //   int getDatatype();


private:
    QList<TreeItem*> childItems;
    QVector<QVariant> itemData;
    TreeItem *parentItem;
};


#endif
