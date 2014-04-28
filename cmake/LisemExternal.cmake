IF(DEFINED ENV{LISEM_EXTERN})
    # Configure search path to find packages.
    SET(CMAKE_PREFIX_PATH
        ${CMAKE_PREFIX_PATH}
        $ENV{LISEM_EXTERN}
    )
ENDIF()


# Configure packages. ----------------------------------------------------------


# Find packages. ---------------------------------------------------------------
FIND_PACKAGE(Qt4 4 REQUIRED QtCore QtGui)
FIND_PACKAGE(Qwt REQUIRED)
FIND_PACKAGE(PCRasterRasterFormat REQUIRED)
FIND_PACKAGE(Doxygen)


# Configure project. -----------------------------------------------------------
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
    SYSTEM
    ${QWT_INCLUDE_DIRS}
    ${PCRASTER_RASTER_FORMAT_INCLUDE_DIRS}
)

SET(LISEM_EXTERNAL_LIBRARIES
    ${QWT_LIBRARIES}
    ${QT_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_LIBRARIES}
)
