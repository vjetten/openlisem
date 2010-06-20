TEMPLATE = app
TARGET = openLisem
QWTDIR = "c:\Qt\qwt-5.2.1"
QT += core \
    gui
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
    include/lookup.h \
    include/lutio.h \
    include/misc.h \
    include/soillut.h \
    include/swat_inp.h \
    include/swatre_p.h \
    include/swatre_g.h
SOURCES += lisTotalsMB.cpp \
    ui_full/LisUItreecheck.cpp \
    ui_full/LisUIModel.cpp \
    ui_full/LisUIrunfile.cpp \
    ui_full/LisUImapnames.cpp \
    ui_full/LisUItreeitem.cpp \
    ui_full/LisUItreemodel.cpp \
    ui_full/LisUIDefaultNames.cpp \
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
    lisSnowmelt.cpp \
    main.cpp \
    mmath.cpp \
    swatre/swatstep.cpp \
    swatre/swatinit.cpp \
    swatre/swat_inp.cpp \
    swatre/soillut.cpp \
    swatre/qsortcmp.cpp \
    swatre/lutio.cpp \
    swatre/lookup.cpp
FORMS += ui_full/lisemqt.ui
CONFIG(debug, debug|release) { 
    LIBS += -L"debug" \
        -llibcsfd \
        -lqwtd5 #\
        #-lswatred
    DESTDIR = debug
}
else { 
    LIBS += -L"bin" \
        -llibcsf \
        -lqwt5 #\
        #-lswatre
    DESTDIR = bin
}
INCLUDEPATH += "include"
INCLUDEPATH += "../libswatreQt/include"
INCLUDEPATH += "ui_full"
INCLUDEPATH += $${QWTDIR}/src
RESOURCES += resources/openlisem.qrc
CONFIG += precompile_header
PRECOMPILED_HEADER = include/stable.h

OTHER_FILES +=
