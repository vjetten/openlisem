IF(DEFINED ENV{LISEM_EXTERN})
    # Configure search path to find packages.
    SET(CMAKE_PREFIX_PATH
        ${CMAKE_PREFIX_PATH}
        $ENV{LISEM_EXTERN}
    )
ENDIF()


# Configure packages. ----------------------------------------------------------


# Find packages. ---------------------------------------------------------------
FIND_PACKAGE(Qt4 4.8 REQUIRED QtCore QtGui)
FIND_PACKAGE(Qwt 6 REQUIRED)
FIND_PACKAGE(Doxygen)

# <HACK>
SET(PCRASTER_RASTER_FORMAT_SOURCE_ROOT
    /home/kor/Development/projects/rasterformat/sources/pcraster_raster_format
)
SET(PCRASTER_RASTER_FORMAT_BINARY_ROOT
    /home/kor/Development/objects/gcc-4_x86-64/Debug/rasterformat/bin
)
# </HACK>

# Configure project. -----------------------------------------------------------
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
    SYSTEM
    ${QWT_INCLUDE_DIRS}
    ${PCRASTER_RASTER_FORMAT_SOURCE_ROOT}
)

SET(LISEM_EXTERNAL_LIBRARIES
    ${QWT_LIBRARIES}
    ${QT_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_BINARY_ROOT}/libpcraster_raster_format.a
)
