cmake_minimum_required(VERSION 3.9)

#============ WIN ========================

IF(WIN32)
   # NOTE: a branch of QWT is used for double axis display:
   # https://sourceforge.net/p/qwt/code/HEAD/tree/branches/qwt-6.1-multiaxes/
   SET(QWT_BUILD_DIR "C:/Qwt61ma/qwt-6.1-ma")       # <= give your own folder names here
   SET(MINGW_BUILD_DIR "C:/msys/mingw64")     		# <= give your own folder names here
       SET(GDAL_BUILD_DIR "C:/msys/mingw64")

    SET(GDAL_INCLUDE_DIRS "${MINGW_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${MINGW_BUILD_DIR}/lib/libgdal.dll.a")

    # QWT standard MSYS install
    #   SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/qwt-qt5")
    #   SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")

    # Lisem uses a QWT branch with quadruple axes support
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/src")
    SET(QWT_LIBRARIES "${CMAKE_CURRENT_SOURCE_DIR}/qwtlib/libqwt.dll.a")

  #  SET(OMP_INCLUDE_DIRS "${MINGW_BUILD_DIR}/lib/gcc/x86_64-w64-mingw32/11.3.0/include")
    FIND_PATH(OMP_INCLUDE_DIRS
        NAMES omp.h
        PATHS "${MINGW_BUILD_DIR}/lib/gcc/x86_64-w64-mingw32"
    )

ENDIF()

#============= LINUX ======================

# linux ubuntu, qwt installation should be in usr if you followed the instructions, version nr may be different
IF(UNIX AND NOT CYGWIN)
 SET(QWT_BUILD_DIR "/usr/local/qwt-6.4.0-svn")
 SET(GDAL_INCLUDE_DIRS "/usr/include/gdal")
 SET(GDAL_LIBRARIES "/usr/lib/libgdal.so")
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include")
    SET(OpenMP_CXX_INCLUDE_DIRS "/usr/lib/gcc/x86_64-linux-gnu/11/include")
ENDIF()

#============ INCLDUES ====================

INCLUDE_DIRECTORIES(
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${OpenMP_CXX_INCLUDE_DIRS}
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
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -std=c++11 ${OpenMP_CXX_FLAGS} -fopenmp")

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
    lisSWOF2D.cpp
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

add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
    ${PCR_SOURCES}
)

target_link_libraries(Lisem Qt5::Widgets Qt5::Gui Qt5::Core ${GDAL_LIBRARIES} ${QWT_LIBRARIES} OpenMP::OpenMP_CXX)

