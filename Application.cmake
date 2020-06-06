cmake_minimum_required(VERSION 2.8.11)

# path to pcraster and qwt build directories on local machine
# example windows
IF(WIN32)
 SET(LISEM_QWT_ROOT "c:/qt/msys64/mingw64")
 SET(PCRASTER_BUILD_DIR "c:/prgc/lisem_external/3rd_party_root-mwqt513/PCR")
 SET(GDAL_BUILD_DIR "c:/prgc/lisem_external/3rd_party_root-mwqt513/GDAL")
ENDIF()

# example linux
IF(UNIX)
 SET(LISEM_QWT_ROOT "/usr/local/qwt-6.1.4")
 SET(PCRASTER_BUILD_DIR "~/pcraster-4.2.1")
ENDIF()

INCLUDE(LisemCompiler)
INCLUDE(LisemExternal)

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

SET(APP_SOURCES
    CsfMap
    error
    fixture
    io
    lisTotalsMB
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
    CsfRGBMap
    lisChannelflow
    lisTiledrainflow
    lisDataInit
    lisErosion
    lisInfiltration
    lisKinematic
    lisModel
    lisOverlandflow
    lisRainintc
    lisReportfile
    lisRunfile
    lisSurfstor
    lisSnowmelt
  #  lisKinematic2d
    main
    operation
    swatre/swatstep
    swatre/swatinit
    swatre/soillut
    swatre/lutio
    swatre/lookup
    swatre/swatinp
    lisChannelflood
    lisSWOF2Dopen
    lisSWOF2Daux
    lisSWOF2DSediment
    lisInterception
    lisPesticide
    lisFlowBarriers
    lisUnifiedFlowThreadPool
    lisUnifiedFlowThread
    lisExtendedChannel
    include/version.h
    include/model.h
    include/TMmapVariables.h
    include/LisUIoutput.h
 #   include/lisUnifiedFlowThreadPool
 #   include/lisUnifiedFlowThread
    include/CsfMap.h
    include/CsfRGBMap.h
    include/array.h
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

    openlisemico.rc

)

QT5_WRAP_UI(UI_SOURCES ui_full/lisemqt.ui)

QT5_ADD_RESOURCES(RCC_SOURCES resources/openlisem.qrc)

# change exec name here
add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
)
# Use the Widgets module from Qt 5.
TARGET_LINK_LIBRARIES(Lisem Qt5::Widgets Qt5::Gui Qt5::Core ${LISEM_EXTERNAL_LIBRARIES})

#TARGET_LINK_LIBRARIES(Lisem ${QWT_LIBRARIES})
