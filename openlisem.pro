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
SOURCES += reportfile.cpp \
    runfile.cpp \
    CsfMap.cpp \
    DataInit.cpp \
    channelflow.cpp \
    erosion.cpp \
    infiltration.cpp \
    kinematic.cpp \
    mmath.cpp \
    overlandflow.cpp \
    rainintc.cpp \
    surfstor.cpp \
    model.cpp \
    main.cpp \
    ifacebasic.cpp
FORMS += ifacebasic.ui
LIBS += -L"d:\prgc\libcsfs\release" \
    -llibcsfs
RESOURCES += 
