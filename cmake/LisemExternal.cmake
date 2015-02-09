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
