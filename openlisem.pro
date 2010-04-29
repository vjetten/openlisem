TEMPLATE = app
TARGET = ifacebasic
QT += core \
    gui
HEADERS += error.h \
    CsfMap.h \
    csf.h \
    csfattr.h \
    csfimpl.h \
    csftypes.h \
    mmath.h \
    model.h \
    global.h \
    ifacebasic.h
SOURCES += lisRunfile.cpp \
    lisReportfile.cpp \
    lisDataInit.cpp \
    lisModel.cpp \
    lisKinematic.cpp \
    lisChannelflow.cpp \
    lisSurfstor.cpp \
    lisOverlandflow.cpp \
    lisRainintc.cpp \
    lisInfiltration.cpp \
    lisErosion.cpp \
    CsfMap.cpp \
    mmath.cpp \
    main.cpp \
    ifacebasic.cpp
FORMS += ifacebasic.ui
LIBS += -L"d:\prgc\libcsfs\release" \
    -llibcsfs
RESOURCES += resources/openlisem.qrc
