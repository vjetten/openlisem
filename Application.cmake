cmake_minimum_required(VERSION 2.8.11)

# path qwt and gdal build directories on local machine
# following windows MSYS2.0 installation
IF(WIN32)
 SET(QWT_BUILD_DIR "c:/qt/qwtma")
 SET(GDAL_BUILD_DIR "c:/qt/msys64/mingw64")
ENDIF()

# linux ubuntu, qwt installation is a bit messy,
# should be in usr if you followed the instructions
IF(UNIX)
 SET(QWT_BUILD_DIR "/usr/local/qwt-6.1.4")
ENDIF()

INCLUDE(LisemCompiler)
INCLUDE(LisemExternal)

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

SET(APP_SOURCES
    main
    CsfMap
    CsfRGBMap
    error
    fixture
    io
    operation
    ui_full/LisUItreecheck
    ui_full/LisUIModel
    ui_full/LisUIrunfile
    ui_full/LisUImapnames
    ui_full/LisUItreeitem
    ui_full/LisUItreemodel
    ui_full/LisUIDefaultNames
    ui_full/lisemqt
    ui_full/LisUIplot
    ui_full/LisUImapplot
    swatre/swatstep
    swatre/swatinit
    swatre/soillut
    swatre/lutio
    swatre/lookup
    swatre/swatinp    
    lisChannelErosion
    lisChannelflood
    lisChannelflow
    lisDataInit
    lisErosion
    lisExtendedChannel
    lisFlowBarriers
    lisEvaporation
    lisInfiltration
    lisInterception
    lisKinematic
    lisModel
    lisOverlandflow
    lisPesticide
    lisPercolation
    lisRainfall
    lisDischargein
    lisReportfile
    lisReportmaps
    lisRunfile
    lisSnowmelt
    lisSurfstor
    lisSWOF2D
    lisSWOF2Daux
    lisSWOF2Dopen
    lisSWOF2DSediment
    lisSWOF2DChannel
    lisTiledrainflow
    lisTotalsMB
    include/array.h
    include/CsfMap.h
    include/CsfRGBMap.h
    include/error.h
    include/fixture.h
    include/global.h
    include/io.h
    include/LisUIoutput.h
    include/masked_raster.h
    include/mmath.h
    include/model.h
    include/operation.h
    include/option.h
    include/raster.h
    include/swatre_g.h
    include/swatre_p.h
    include/swatreLookup.h
    include/swatremisc.h
    include/swatresoillut.h
    include/TMmapVariables.h
    include/version.h
    include/pcrtypes.h
    include/csf.h
    include/csfattr.h
    include/csfimpl.h
    include/csftypes.h
    openlisemico.rc
)

SET(PCR_SOURCES
    PCRlib/_getcell.c
    PCRlib/_getrow.c
    PCRlib/_gsomece.c
    PCRlib/_putcell.c
    PCRlib/_rputrow.c
    PCRlib/angle.c
    PCRlib/attravai.c
    PCRlib/attrsize.c
    PCRlib/cellsize.c
    PCRlib/create2.c
    PCRlib/csfglob.c
    PCRlib/csfsup.c
    PCRlib/delattr.c
    PCRlib/dumconv.c
    PCRlib/endian.c
    PCRlib/filename.c
    PCRlib/gattrblk.c
    PCRlib/gattridx.c
    PCRlib/gcellrep.c
    PCRlib/gdattype.c
    PCRlib/getattr.c
    PCRlib/getx0.c
    PCRlib/gety0.c
    PCRlib/ggisfid.c
    PCRlib/gmaxval.c
    PCRlib/gminval.c
    PCRlib/gnrcols.c
    PCRlib/gnrrows.c
    PCRlib/gproj.c
    PCRlib/gputproj.c
    PCRlib/gvalscal.c
    PCRlib/gvartype.c
    PCRlib/gversion.c
    PCRlib/ismv.c
    PCRlib/kernlcsf.c
    PCRlib/legend.c
    PCRlib/mclose.c
    PCRlib/mopen.c
    PCRlib/moreattr.c
    PCRlib/mperror.c
    PCRlib/pgisfid.c
    PCRlib/pmaxval.c
    PCRlib/pminval.c
    PCRlib/putallmv.c
    PCRlib/putattr.c
    PCRlib/putsomec.c
    PCRlib/putx0.c
    PCRlib/puty0.c
    PCRlib/pvalscal.c
    PCRlib/rattrblk.c
    PCRlib/rcomp.c
    PCRlib/rcoords.c
    PCRlib/rdup2.c
    PCRlib/reseterr.c
    PCRlib/rextend.c
    PCRlib/rmalloc.c
    PCRlib/rrowcol.c
    PCRlib/ruseas.c
    PCRlib/setangle.c
    PCRlib/setmv.c
    PCRlib/setvtmv.c
    PCRlib/strconst.c
    PCRlib/strpad.c
    PCRlib/swapio.c
    PCRlib/trackmm.c
    PCRlib/vs2.c
    PCRlib/vsdef.c
    PCRlib/vsis.c
    PCRlib/vsvers.c
    PCRlib/wattrblk.c
)

QT5_WRAP_UI(UI_SOURCES ui_full/lisemqt.ui)

QT5_ADD_RESOURCES(RCC_SOURCES resources/openlisem.qrc)

# change exec name here
add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
    ${PCR_SOURCES}
)
# Use the Widgets module from Qt 5.
TARGET_LINK_LIBRARIES(Lisem Qt5::Widgets Qt5::Gui Qt5::Core  ${LISEM_EXTERNAL_LIBRARIES})
#Qt5::Charts
