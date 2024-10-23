cmake_minimum_required(VERSION 3.9)

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()

# Platform-specific configurations
IF(WIN32)
    # QWT configuration for double axis display, note a double axis branch of qwt is used
    SET(QWT_BUILD_DIR "C:/prgc/lisemgit/qwt/git")    # Adjust to your folder names
    SET(MINGW_BUILD_DIR "c:/qt/msys64/mingw64")     # Adjust to your folder names
    SET(GDAL_INCLUDE_DIRS "${MINGW_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${MINGW_BUILD_DIR}/lib/libgdal.dll.a")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/src")
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")

    FIND_PATH(OMP_INCLUDE_DIRS
        NAMES omp.h
        PATHS "${MINGW_BUILD_DIR}/lib/gcc/x86_64-w64-mingw32/14.1.0/include"
    )
ENDIF()

IF(UNIX AND NOT CYGWIN)
    SET(QWT_BUILD_DIR "/usr/local/qwt-6.1.4")
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/")
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${OMP_INCLUDE_DIRS}
    SYSTEM
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Enable automatic handling of MOC, UIC, and RCC
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTOUIC ON)
#set(CMAKE_AUTORCC ON)

# Optionally skip rule dependency checks to avoid timestamp issues
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set_property(DIRECTORY PROPERTY CMAKE_CONFIGURE_DEPENDS "")

# Cache the autogen output in a specific directory
set(CMAKE_AUTOGEN_BUILD_DIR "${CMAKE_BINARY_DIR}/autogen")
set_source_files_properties(ui_full/LisUIDialogs.cpp PROPERTIES AUTOMOC ON)
set_source_files_properties(ui_full/lisemqt.cpp PROPERTIES AUTOMOC ON)
set_source_files_properties(ui_full/LisUItreemodel.cpp PROPERTIES AUTOMOC ON)
set_source_files_properties(ui_full/Lismpeg.cpp PROPERTIES AUTOUIC ON)
set_source_files_properties(resources/openlisem.qrc PROPERTIES AUTORCC ON)

set_property(SOURCE ui_full/lisemqt.ui PROPERTY SKIP_AUTOUIC ON)

# Compiler flags
IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -Wno-var-tracking-assignments -std=c++11 ${OpenMP_CXX_FLAGS}")
    IF(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wl,-rpath=${ORIGIN}./lib")
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
    #include/model.h
    include/operation.h
    include/option.h
    include/raster.h
    include/swatre_p.h
    include/TMmapVariables.h
    include/VectormapVariables.h
    include/version.h
    openlisemico.rc
)

qt6_wrap_cpp(MOC_FILES
    include/model.h
    #include/lisemqt.h
    #include/lismpeg.h
    # Add all header files with Q_OBJECT here
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
    ${MOC_FILES}
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


