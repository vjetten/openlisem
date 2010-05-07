/****************************************************************************
** Meta object code from reading C++ file 'ifacebasic.h'
**
** Created: Fri May 7 08:56:08 2010
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../ifacebasic.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ifacebasic.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ifacebasic[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      35,   11,   11,   11, 0x08,
      71,   11,   11,   11, 0x08,
      88,   80,   11,   11, 0x08,
     107,   80,   11,   11, 0x08,
     127,   11,   11,   11, 0x08,
     163,   11,   11,   11, 0x08,
     175,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ifacebasic[] = {
    "ifacebasic\0\0on_runButton_clicked()\0"
    "on_toolButton_ShowRunfile_clicked()\0"
    "Showit()\0results\0worldDone(QString)\0"
    "worldDebug(QString)\0"
    "on_toolButton_runfilename_clicked()\0"
    "StorePath()\0GetStorePath()\0"
};

const QMetaObject ifacebasic::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ifacebasic,
      qt_meta_data_ifacebasic, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ifacebasic::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ifacebasic::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ifacebasic::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ifacebasic))
        return static_cast<void*>(const_cast< ifacebasic*>(this));
    if (!strcmp(_clname, "Ui::ifacebasicClass"))
        return static_cast< Ui::ifacebasicClass*>(const_cast< ifacebasic*>(this));
    return QWidget::qt_metacast(_clname);
}

int ifacebasic::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: on_runButton_clicked(); break;
        case 1: on_toolButton_ShowRunfile_clicked(); break;
        case 2: Showit(); break;
        case 3: worldDone((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: worldDebug((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: on_toolButton_runfilename_clicked(); break;
        case 6: StorePath(); break;
        case 7: GetStorePath(); break;
        default: ;
        }
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
