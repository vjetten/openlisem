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
SOURCES += source/CsfMap.cpp \
    source/ifacebasic.cpp \
    source/lisChannelflow.cpp \
    source/lisDataInit.cpp \
    source/lisErosion.cpp \
    source/lisInfiltration.cpp \
    source/lisKinematic.cpp \
    source/lisModel.cpp \
    source/lisOverlandflow.cpp \
    source/lisRainintc.cpp \
    source/lisReportfile.cpp \
    source/lisRunfile.cpp \
    source/lisSurfstor.cpp \
    source/main.cpp \
    source/mmath.cpp
FORMS += ui/ifacebasic.ui
# LIBS += -L"d:\prgc\libcsfs\debug" \
# -llibcsfs
LIBS += -L"d:\prgc\libcsf\debug" \
    -llibcsf
INCLUDEPATH += "include"
RESOURCES += resources/openlisem.qrc
CONFIG += precompile_header
PRECOMPILED_HEADER = include/stable.h
