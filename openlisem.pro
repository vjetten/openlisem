TEMPLATE = app
TARGET = openLisem
QT += core \
    gui
HEADERS += ui_full/LisUItreeitem.h \
    ui_full/LisUItreemodel.h \
    ui_full/lisemqt.h \
    include/LisUIoutput.h \
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
SOURCES += ui_full/LisUItreecheck.cpp \
    ui_full/LisUIModel.cpp \
    ui_full/LisUIrunfile.cpp \
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
FORMS += ui_full/lisemqt.ui 
    #ui/ifacebasic.ui
CONFIG(debug, debug|release) {
LIBS += -L"libs" \
    -llibcsfd \
    -lqwtd5  
} else {
LIBS += -L"libs" \
    -llibcsf \
    -lqwt5 
}
INCLUDEPATH += "include"
INCLUDEPATH += "ui_full"
INCLUDEPATH += "d:\prgc\qwt-5.2.1\src"
RESOURCES += resources/openlisem.qrc
CONFIG += precompile_header
PRECOMPILED_HEADER = include/stable.h
CONFIG(debug, debug|release) {
    DESTDIR = debug
} else {
    DESTDIR = bin
}
