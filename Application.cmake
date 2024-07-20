cmake_minimum_required(VERSION 3.9)

#============ WIN ========================

IF(WIN32)
   # NOTE: a branch of QWT is used for double axis display:
   # https://sourceforge.net/p/qwt/code/HEAD/tree/branches/qwt-6.1-multiaxes/
    #SET(QWT_BUILD_DIR "c:/qt/qwt-6.1-ma")          # <= give your own folder names here
    #SET(QWT_BUILD_DIR "c:/qt/qwt-6.3.0")          # <= give your own folder names here
    SET(QWT_BUILD_DIR "C:/prgc/lisemgit/qwt/git")          # <= give your own folder names here
    SET(MINGW_BUILD_DIR "c:/qt/msys64/mingw64")     # <= give your own folder names here

    SET(GDAL_INCLUDE_DIRS "${MINGW_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${MINGW_BUILD_DIR}/lib/libgdal.dll.a")

    # Lisem uses a QWT branch with quadruple axes support
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/src")
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")
   # SET(QWT_LIBRARIES "${MINGW_BUILD_DIR}/lib/libqwt-qt6.dll.a")

    #SET(PCR_LIBRARIES "C:/prgc/PCR/libs/sources/libpcraster_raster_format.a")

    FIND_PATH(OMP_INCLUDE_DIRS
        NAMES omp.h
        PATHS "${MINGW_BUILD_DIR}/lib/gcc/x86_64-w64-mingw32/14.1.0/include"
    )

ENDIF()

#============= LINUX ======================

# linux ubuntu, qwt installation should be in usr if you followed the instructions, version nr may be different
IF(UNIX AND NOT CYGWIN)
    SET(QWT_BUILD_DIR "/usr/local/qwt-6.1.4")
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/")
ENDIF()

#============ INCLDUES ====================

INCLUDE_DIRECTORIES(
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${OMP_INCLUDE_DIRS}
    SYSTEM
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

#============ OMP ===========================

find_package(OpenMP REQUIRED)

#============ FLAGS =========================

INCLUDE(CheckCXXCompilerFlag)

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -Wno-var-tracking-assignments -std=c++11 ${OpenMP_CXX_FLAGS}")
#--param=max-vartrack-size=1500000
    #add_compile_options(-Wno-var-tracking-assignments)

    IF(UNIX)
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wl,-rpath=${ORIGIN}./lib")
       # extra flags for thread and looking for so libs in ./lib
    ENDIF()
ENDIF()

#============ sourcecode files ===============

SET(APP_SOURCES
    fixesandbugs.txt    
    main.cpp
    CsfMap.cpp
    CsfRGBMap.cpp
    error.cpp
    fixture.cpp
    io.cpp
    operation.cpp
    ui_full/LisUIDialogs.cpp
    ui_full/LisUIScreenshot.cpp
    ui_full/LisUItreecheck.cpp
    ui_full/LisUIModel.cpp
    ui_full/LisUIrunfile.cpp
    ui_full/LisUImapnames.cpp
    ui_full/LisUItreeitem.cpp
    ui_full/LisUItreemodel.cpp
    ui_full/LisUIDefaultNames.cpp
    ui_full/lisemqt.cpp
    ui_full/LisUIplot.cpp
    ui_full/LisUImapplot.cpp
    ui_full/LisUImapplot.h
    ui_full/Lismpeg.cpp
    ui_full/lisUIStyle.cpp
    ui_full/Lismpeg.h
    ui_full/lisemqt.h
    swatre/swatstep.cpp
    swatre/swatinit.cpp
    swatre/soillut.cpp
    swatre/lutio.cpp
    swatre/lookup.cpp
    swatre/swatinp.cpp    
    lisBoundary.cpp
    lisChannelErosion.cpp
    lisChannelflood.cpp
    lisChannelflow.cpp
    lisDataInit.cpp
    lisDataFunctions.cpp
    lisErosion.cpp
    lisExtendedChannel.cpp
    lisFlowBarriers.cpp
    lisEvaporation.cpp
    lisGWflow.cpp
    lisInfiltration.cpp
    lisInterception.cpp
    lisKinematic.cpp
    lisModel.cpp
    lisOverlandflow.cpp
    lisPesticide.cpp
    lisPercolation.cpp
    lisRainfall.cpp
    lisDischargein.cpp
    lisReportfile.cpp
    lisReportmaps.cpp
    lisRunfile.cpp
    lisSnowmelt.cpp
    lisSurfstor.cpp
    lisSoilmoisture.cpp
   # lisSWOF2D.cpp  OBSOLETE
    lisSWOF2Daux.cpp
    lisSWOF2Dopen.cpp
    lisSWOF2DSediment.cpp
    lisSWOF2DChannel.cpp
    lisTiledrainflow.cpp
    lisTotalsMB.cpp
    include/array.h
    include/CsfMap.h
    include/CsfRGBMap.h
    include/lerror.h
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
    include/VectormapVariables.h
    include/version.h
    #include/pcrtypes.h
    #include/csf.h
    #include/csfattr.h
    #include/csfimpl.h
    #include/csftypes.h
    openlisemico.rc
)

# Generate UI source files
qt_wrap_ui(UI_SOURCES ui_full/lisemqt.ui ui_full/lismpeg.ui)

# Generate resource source files
qt_add_resources(RCC_SOURCES resources/openlisem.qrc)

# Add executable target
add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
)

# Link the necessary libraries
if (Qt6_FOUND)
    target_link_libraries(Lisem
        Qt6::Widgets Qt6::Gui Qt6::Core
        ${GDAL_LIBRARIES} ${QWT_LIBRARIES}
        OpenMP::OpenMP_CXX
    )
else()
    target_link_libraries(Lisem
        Qt5::Widgets Qt5::Gui Qt5::Core
        ${GDAL_LIBRARIES} ${QWT_LIBRARIES}
        OpenMP::OpenMP_CXX
    )
endif()

#QT5_WRAP_UI(UI_SOURCES ui_full/lisemqt.ui ui_full/lismpeg.ui)
#
#QT5_ADD_RESOURCES(RCC_SOURCES resources/openlisem.qrc)
#
#
#target_link_libraries(Lisem
#    Qt5::Widgets Qt5::Gui Qt5::Core
#    ${GDAL_LIBRARIES} ${QWT_LIBRARIES}
#    OpenMP::OpenMP_CXX
#)
#${PCR_LIBRARIES}
