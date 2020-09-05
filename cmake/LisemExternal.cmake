
# Find packages. ---------------------------------------------------------------
# this calls fndq
FIND_PACKAGE(GDAL REQUIRED)
FIND_PACKAGE(Qwt REQUIRED)
FIND_PACKAGE(PCRasterRasterFormat REQUIRED)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OpenMP::OpenMP_CXX)
endif()

# Configure project. -----------------------------------------------------------

INCLUDE_DIRECTORIES(
    SYSTEM
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${PCRASTER_RASTER_FORMAT_INCLUDE_DIRS}
    ${OMP_INCLUDE_DIRS}
)


SET(LISEM_EXTERNAL_LIBRARIES
    ${GDAL_LIBRARIES}
    ${QWT_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_LIBRARIES}
)
