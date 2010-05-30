TEMPLATE = app
TARGET = openLisem
QT += core \
    gui
HEADERS += ui/ifacebasic.h \
    ui_full/LisUItreeitem.h \
    ui_full/LisUItreemodel.h \
    ui_full/lisemqt.h \
    include/version.h \
    include/CsfMap.h \
    include/csf.h \
    include/csfattr.h \
    include/csfimpl.h \
    include/csftypes.h \
    include/error.h \
    include/global.h \
    include/mmath.h \
    include/model.h \
    include/stable.h
SOURCES += ui/ifacebasic.cpp \
    ui_full/LisUImapnames.cpp \
    ui_full/LisUItreeitem.cpp \
    ui_full/LisUItreemodel.cpp \
    ui_full/lisemqt.cpp \
    CsfMap.cpp \
    lisChannelflow.cpp \
    lisDataInit.cpp \
    lisErosion.cpp \
    lisInfiltration.cpp \
    lisKinematic.cpp \
    lisModel.cpp \
    lisOverlandflow.cpp \
    lisRainintc.cpp \
    lisReportfile.cpp \
    lisRunfile.cpp \
    lisSurfstor.cpp \
    main.cpp \
    mmath.cpp
FORMS += ui_full/lisemqt.ui \
    ui/ifacebasic.ui
LIBS += -L"libs/debug/" \
    -llibcsf
INCLUDEPATH += "include"
INCLUDEPATH += "ui_full"
RESOURCES += resources/openlisem.qrc
CONFIG += precompile_header
PRECOMPILED_HEADER = include/stable.h
