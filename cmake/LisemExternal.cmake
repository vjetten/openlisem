# Configure packages. ----------------------------------------------------------
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_STATIC_RUNTIME OFF)
ADD_DEFINITIONS(
    # Use dynamic libraries.
    -DBOOST_ALL_DYN_LINK
    # Prevent auto-linking.
    -DBOOST_ALL_NO_LIB
)


# Find packages. ---------------------------------------------------------------
FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
FIND_PACKAGE(GDAL REQUIRED)
FIND_PACKAGE(Qt4 4 REQUIRED QtCore QtGui)
FIND_PACKAGE(Qwt REQUIRED)
FIND_PACKAGE(PCRasterRasterFormat REQUIRED)
FIND_PACKAGE(Fern REQUIRED)
FIND_PACKAGE(Doxygen)


# Fixup GDAL variable.
IF(WIN32 AND MINGW)
    # On Windows, find_package always selects the release dll named
    # libgdal.dll. This is a symbolic links which trips the linker. Also,
    # for debug builds, we want to link against the debug library.
    # Before fixing this, verify that the variable contains a single string
    # and that the name is libgdal.dll. The CMake module may be fixed at
    # some point.
    LIST(LENGTH GDAL_LIBRARIES NR_GDAL_LIBRARIES)
    IF(NR_GDAL_LIBRARIES EQUAL 1)
        GET_FILENAME_COMPONENT(GDAL_LIBRARY_NAME ${GDAL_LIBRARIES} NAME)

        IF(GDAL_LIBRARY_NAME STREQUAL "libgdal.dll")
            SET(GDAL_VERSION "1.11.1")
            SET(OPTIMIZED_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            SET(DEBUG_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            STRING(REPLACE gdal gdald DEBUG_GDAL_LIBRARY ${DEBUG_GDAL_LIBRARY})
            SET(GDAL_LIBRARIES
                optimized ${OPTIMIZED_GDAL_LIBRARY}
                debug ${DEBUG_GDAL_LIBRARY})
        ENDIF()
    ENDIF()
ENDIF()


# Configure project. -----------------------------------------------------------
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
    SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
)

INCLUDE_DIRECTORIES(
    ${PCRASTER_RASTER_FORMAT_INCLUDE_DIRS}
    ${FERN_INCLUDE_DIRS}
)

SET(LISEM_EXTERNAL_LIBRARIES
    ${QWT_LIBRARIES}
    ${QT_LIBRARIES}
    ${GDAL_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_LIBRARIES}
    ${FERN_LIBRARIES}
)
