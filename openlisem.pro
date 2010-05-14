TEMPLATE = app
TARGET = openLisem
QT += core \
    gui
HEADERS += include/CsfMap.h \
    include/csf.h \
    include/csfattr.h \
    include/csfimpl.h \
    include/csftypes.h \
    include/error.h \
    include/global.h \
    include/ifacebasic.h \
    include/mmath.h \
    include/model.h \
    include/stable.h
SOURCES += CsfMap.cpp \
    ifacebasic.cpp \
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
FORMS += ui/ifacebasic.ui
LIBS += -L"libs/debug/" \
    -llibcsf
INCLUDEPATH += "include"
RESOURCES += resources/openlisem.qrc
#CONFIG += precompile_header
#PRECOMPILED_HEADER = include/stable.h
