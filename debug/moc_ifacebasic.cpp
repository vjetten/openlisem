/****************************************************************************
** Meta object code from reading C++ file 'ifacebasic.h'
**
** Created: Wed Mar 3 07:22:10 2010
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
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      40,   35,   11,   11, 0x08,
      60,   52,   11,   11, 0x08,
      79,   52,   11,   11, 0x08,
      99,   11,   11,   11, 0x08,
     125,   11,   11,   11, 0x08,
     151,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_ifacebasic[] = {
    "ifacebasic\0\0on_runButton_clicked()\0"
    "step\0Showit(int)\0results\0worldDone(QString)\0"
    "worldDebug(QString)\0on_checkChannel_clicked()\0"
    "on_checkErosion_clicked()\0"
    "on_toolButton_runfilename_clicked()\0"
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
        case 1: Showit((*reinterpret_cast< const int(*)>(_a[1]))); break;
        case 2: worldDone((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: worldDebug((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: on_checkChannel_clicked(); break;
        case 5: on_checkErosion_clicked(); break;
        case 6: on_toolButton_runfilename_clicked(); break;
        default: ;
        }
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
