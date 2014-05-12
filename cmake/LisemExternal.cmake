# Configure packages. ----------------------------------------------------------


# Find packages. ---------------------------------------------------------------
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(Qt4 4 REQUIRED QtCore QtGui)
FIND_PACKAGE(Qwt REQUIRED)
FIND_PACKAGE(PCRasterRasterFormat REQUIRED)
FIND_PACKAGE(Fern REQUIRED)
FIND_PACKAGE(Doxygen)


# Configure project. -----------------------------------------------------------
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
    SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
)

INCLUDE_DIRECTORIES(
    ${PCRASTER_RASTER_FORMAT_INCLUDE_DIRS}
    ${FERN_INCLUDE_DIRS}
)

SET(LISEM_EXTERNAL_LIBRARIES
    ${QWT_LIBRARIES}
    ${QT_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_LIBRARIES}
    ${FERN_LIBRARIES}
)
