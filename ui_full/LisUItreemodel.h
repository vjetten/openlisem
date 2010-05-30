
#ifndef TREEMODEL_H
#define TREEMODEL_H

#include <QtGui>

//#include <QAbstractItemModel>
//#include <QModelIndex>
//#include <QVariant>

#include "LisUItreeitem.h"

class TreeItem;

//! [0]
class TreeModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    TreeModel(const QStringList &headers, const QStringList &data,
              QObject *parent = 0);
    ~TreeModel();
//! [0]

//! [1]
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;

    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;

    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
//! [1]

//! [2]
    Qt::ItemFlags flags(const QModelIndex &index) const;
    bool setData(const QModelIndex &index, const QVariant &value,
                 int role = Qt::EditRole);
    bool setHeaderData(int section, Qt::Orientation orientation,
                       const QVariant &value, int role = Qt::EditRole);
    void setflag(int f, int row, const QModelIndex &parent = QModelIndex()) const;
    bool getflag(int row, const QModelIndex &parent = QModelIndex()) const;

private:
    void setupModelData(const QStringList &lines, TreeItem *parent);
    TreeItem *getItem(const QModelIndex &index) const;

    TreeItem *rootItem;
};
//! [2]



#endif
