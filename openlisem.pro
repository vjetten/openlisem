TEMPLATE = app
TARGET = openLisem
QT += core \
    gui
QWTDIR = c:/Qt/qwt
#CSFDIR = d:/prgc/libcsf/
#change the QWT directory to your own install
# CONFIG += console
# this echos qdebug to screen in dos mode
CONFIG += exceptions
HEADERS += \
    ui_full/LisUItreeitem.h \
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
    include/stable.h \
    include/swatre_p.h \
    include/swatreLookup.h \
    include/swatresoillut.h \
    include/swatre_g.h \
    include/TMmapVariables.h \
    ui_full/LisUImapplot.h
SOURCES += lisTotalsMB.cpp \
    ui_full/LisUItreecheck.cpp \
    ui_full/LisUIModel.cpp \
    ui_full/LisUIrunfile.cpp \
    ui_full/LisUImapnames.cpp \
    ui_full/LisUItreeitem.cpp \
    ui_full/LisUItreemodel.cpp \
    ui_full/LisUIDefaultNames.cpp \
    ui_full/lisemqt.cpp \
    ui_full/LisUIplot.cpp \
    CsfMap.cpp \
    lisChannelflow.cpp \
    lisTiledrainflow.cpp \
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
    lisSnowmelt.cpp \
    main.cpp \
    mmath.cpp \
    swatre/swatstep.cpp \
    swatre/swatinit.cpp \
    swatre/soillut.cpp \
    swatre/lutio.cpp \
    swatre/lookup.cpp \
    swatre/swatinp.cpp \
    LisKinematicSorted.cpp \
    ui_full/LisUImapplot.cpp \
    lisChannelflood.cpp \
    lisSWOF2D.cpp
FORMS += ui_full/lisemqt.ui
# debug version
CONFIG(debug, debug|release) {
  win32:win32-msvc2010{
    # VC version
    DEFINES += _CRT_SECURE_NO_WARNINGS
    LIBS += -L"debug/vc/static" -lcsfvcd -lqwtd
 #   LIBS += -L"debug/vc" -lqwtvcd  #-L"$${QWTDIR}/lib"
    DESTDIR = debug/vc
  }else{
    # mingw version
    LIBS += -L"debug" -llibcsfd
    LIBS += -L"$${QWTDIR}/lib" -lqwtd
    DESTDIR = debug
  }
  MOC_DIR = debug/moc
  OBJECTS_DIR= debug/objs
  UI_DIR= debug/ui
  RCC_DIR = debug/rc
}
#release version
else {
  win32:win32-msvc2010{
    DEFINES += _CRT_SECURE_NO_WARNINGS
    LIBS += -L"bin/vc/static" -lcsfvc -lqwt
    #LIBS += -L"bin/vc" -L"$${QWTDIR}/lib" -lqwtvc
    DESTDIR = bin/vc
  }else{
    # mingw version
    LIBS += -L"bin" -llibcsf
    LIBS += -L"$${QWTDIR}/lib" -lqwt
    DESTDIR = bin
  }
  MOC_DIR = bin/moc
  OBJECTS_DIR= bin/objs
  UI_DIR = bin/ui
  RCC_DIR = bin/rc
}
INCLUDEPATH += "include"
INCLUDEPATH += "ui_full"
INCLUDEPATH += $${QWTDIR}/src
RESOURCES += resources/openlisem.qrc
CONFIG += precompile_header
PRECOMPILED_HEADER = include/stable.h
RC_FILE = openlisemico.rc
#OTHER_FILES +=
