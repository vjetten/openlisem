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
    lisInfiltration
    lisInterception
    lisKinematic
    lisModel
    lisOverlandflow
    lisPesticide
    lisRainintc
    lisReportfile
    lisRunfile
    lisSnowmelt
    lisSurfstor
    lisSWOF2Daux
    lisSWOF2Dopen
    lisSWOF2DSediment
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

