IF(UNIX AND NOT CYGWIN)
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

    SET(PCRASTER_RASTER_FORMAT_LIBRARIES "${PCRASTER_BUILD_DIR}/build/bin/libpcraster_raster_format.a")
    SET(PCRASTER_RASTER_FORMAT_INCLUDE_DIRS "${PCRASTER_BUILD_DIR}/source/rasterformat/sources/pcraster_raster_format")
    SET(QWT_LIBRARIES "${LISEM_QWT_ROOT}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${LISEM_QWT_ROOT}/include/")
ENDIF()

IF(WIN32)
    SET(PCRASTER_RASTER_FORMAT_LIBRARIES "${PCRASTER_BUILD_DIR}/lib/libpcraster_raster_format.a")
    SET(PCRASTER_RASTER_FORMAT_INCLUDE_DIRS "${PCRASTER_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${GDAL_BUILD_DIR}/lib/libgdal.dll")
    SET(GDAL_INCLUDE_DIRS "${GDAL_BUILD_DIR}/include")
    SET(QWT_INCLUDE_DIRS "${LISEM_QWT_ROOT}/include/qwt")
    SET(QWT_LIBRARIES "${LISEM_QWT_ROOT}/lib/libqwt.dll.a")
ENDIF()


INCLUDE(CheckCXXCompilerFlag)

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -std=c++11 -fopenmp")
    IF(UNIX)
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wl,-rpath=${ORIGIN}./lib")
       # extra flags for thread and looking for so libs in ./lib
    ENDIF()
ENDIF()

